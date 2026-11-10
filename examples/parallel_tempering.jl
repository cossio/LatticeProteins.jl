# Reproducible full-model experiment (structure A, beta_max=3000).
# Run from the repository root:
# julia --project=. --threads=4 examples/parallel_tempering.jl OUTPUT [RUNS [BURNIN [NSEQS [THIN [all|shakhnovich]]]]]
# Defaults: 2 independent ladders, 5000 burn-in rounds, 500 samples, thin=10.
# Each round has 10 mutation proposals per replica. CSVs include cold sequences,
# pnat/replica-ID traces, and per-edge acceptance; summary.toml records parameters.
# Controls are independent beta=3000 Metropolis chains (`metropolis!`), which make the
# same local moves as the replicas of a ladder, with the same total proposal budget as
# each ladder (one control chain per temperature). They isolate the effect of exchanges.
# The default is all 103,346 structures; the optional final argument selects the
# 10,000-structure subset. Target indices are translated to the selected set.

using LatticeProteins: sample_msa, ParallelTempering, metropolis!, pnat, aaseq, structure_index,
    STRUCTURES, shakhnovich_structures, STRUCTURE_A, L
using Random: Xoshiro
using DelimitedFiles: writedlm
using TOML: print as print_toml

function save_csv(path, header, rows)
    return open(path, "w") do io
        println(io, join(header, ','))
        writedlm(io, rows, ',')
    end
end

function controls(initial, seed, burnin, nseqs, thin, local_steps, s, structures)
    rng = Xoshiro(seed)
    ps = Matrix{Float64}(undef, size(initial, 2), nseqs)
    for i in axes(initial, 2)
        seq = initial[:, i]
        metropolis!(seq, s; β = 3000, nsteps = burnin * local_steps, structures, rng)
        for sample in 1:nseqs
            metropolis!(seq, s; β = 3000, nsteps = thin * local_steps, structures, rng)
            ps[i, sample] = pnat(s, seq; structures)
        end
        println("control seed=$seed chain=$i final_pnat=$(ps[i, end])")
        flush(stdout)
    end
    return ps
end

function main(args)
    isempty(args) && error("provide an output directory")
    out = abspath(args[1])
    ispath(out) && error("output path already exists: $out")
    runs = length(args) >= 2 ? parse(Int, args[2]) : 2
    burnin = length(args) >= 3 ? parse(Int, args[3]) : 5000
    nseqs = length(args) >= 4 ? parse(Int, args[4]) : 500
    thin = length(args) >= 5 ? parse(Int, args[5]) : 10
    which = length(args) >= 6 ? args[6] : "all"
    which in ("all", "shakhnovich") || error("structure set must be all or shakhnovich")
    structures = which == "all" ? STRUCTURES : shakhnovich_structures()
    s = structure_index(structures, STRUCTURE_A) # index of structure A in the chosen set
    runs > 0 && burnin >= 0 && nseqs > 0 && thin > 0 || error("invalid run lengths")
    local_steps = 10
    betas = [0.0; exp.(range(log(0.1), log(3000.0); length = 31))]
    betas[end] = 3000.0
    mkdir(out)
    initial = [rand(Xoshiro(1000 + run), 1:20, L, length(betas)) for run in 1:runs]
    summaries = Vector{Dict{String, Any}}(undef, 2runs)
    Threads.@threads for job in 1:2runs
        run = (job + 1) ÷ 2
        seed = 2000 + run
        if isodd(job)
            # Independent jobs already run in parallel; keep each ladder serial.
            method = ParallelTempering(; betas, burnin, thin, local_steps, threaded = false)
            elapsed = @elapsed r = sample_msa(
                s, method; nseqs, structures, initial_seqs = initial[run], rng = Xoshiro(seed)
            )
            rows = [
                (
                        sample * thin, exp(r.logpnat[end, sample]), r.replica_ids[end, sample],
                        aaseq(r.samples[:, sample]),
                    ) for sample in 1:nseqs
            ]
            save_csv(
                joinpath(out, "pt_$(run)_cold.csv"),
                ["production_round", "pnat", "replica_id", "sequence"], rows
            )
            save_csv(
                joinpath(out, "pt_$(run)_trace.csv"),
                ["production_round", "temperature_index", "beta", "pnat", "replica_id"],
                [
                    (sample * thin, i, betas[i], exp(r.logpnat[i, sample]), r.replica_ids[i, sample])
                        for sample in 1:nseqs for i in eachindex(betas)
                ]
            )
            d = r.diagnostics
            rates = d.swap_accepts ./ d.swap_attempts
            save_csv(
                joinpath(out, "pt_$(run)_swaps.csv"),
                ["beta_hot", "beta_cold", "attempts", "accepts", "acceptance"],
                [
                    (betas[i], betas[i + 1], d.swap_attempts[i], d.swap_accepts[i], rates[i])
                        for i in eachindex(rates)
                ]
            )
            summaries[job] = Dict(
                "method" => "parallel_tempering", "run" => run,
                "seed" => seed, "seconds" => elapsed,
                "pnat_min" => minimum(exp.(r.logpnat[end, :])),
                "pnat_max" => maximum(exp.(r.logpnat[end, :])),
                "pnat_mean" => sum(exp.(r.logpnat[end, :])) / nseqs,
                "cold_replica_ids" => sort(unique(r.replica_ids[end, :])),
                "round_trips" => d.round_trips,
                "burnin_round_trips" => r.burnin_diagnostics.round_trips,
                "swap_acceptance_min" => minimum(rates),
                "swap_acceptance_max" => maximum(rates)
            )
            println("PT run=$run mean_pnat=$(summaries[job]["pnat_mean"]) round_trips=$(sum(d.round_trips)) seconds=$elapsed")
        else
            elapsed = @elapsed ps = controls(initial[run], seed, burnin, nseqs, thin, local_steps, s, structures)
            save_csv(
                joinpath(out, "control_$(run)_trace.csv"),
                ["production_round", "chain", "pnat"],
                [(sample * thin, i, ps[i, sample]) for sample in 1:nseqs for i in axes(ps, 1)]
            )
            summaries[job] = Dict(
                "method" => "metropolis_control", "run" => run,
                "seed" => seed, "seconds" => elapsed, "pnat_min" => minimum(ps),
                "pnat_max" => maximum(ps), "pnat_mean" => sum(ps) / length(ps),
                "final_pnat" => ps[:, end]
            )
        end
        flush(stdout)
    end
    open(joinpath(out, "summary.toml"), "w") do io
        print_toml(
            io, Dict(
                "structure" => s, "structure_set" => which,
                "nstructures" => length(structures), "betas" => betas,
                "runs" => runs, "burnin_rounds" => burnin, "nseqs" => nseqs,
                "thin_rounds" => thin, "local_steps" => local_steps,
                "proposals_per_method_per_run" => length(betas) * (burnin + nseqs * thin) * local_steps,
                "julia_version" => string(VERSION), "results" => summaries
            )
        )
    end
    return println("Results written to $out")
end

main(ARGS)
