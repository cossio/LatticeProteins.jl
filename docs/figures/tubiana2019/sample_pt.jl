# Reproduce the full-structure alignment with independent replica-exchange ladders.
# Usage: julia -t 4 --project=. sample_pt.jl OUT [RUNS NSEQS BETA BURNIN THIN LOCAL_STEPS CACHE]
# Defaults: 4 ladders, 9000 saved sequences each, beta_max=1000, 4000 burn-in
# rounds, saving every 2 rounds, 5 mutations per replica per round.
# Keep the unfiltered trajectories: convergence must be checked before filtering.
using LatticeProteins: LatticeProteins, parallel_tempering, STRUCTURES, CONTACT_MAP_A,
    L, aaseq, log_pnat, Hugo_MSA_path
using Random: Xoshiro
using DelimitedFiles: writedlm
using Statistics: mean
using TOML
using SHA: sha256

function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ','))
        writedlm(io, rows, ',')
    end
end

function main(args)
    isempty(args) && error("provide a new output directory")
    out = abspath(args[1])
    ispath(out) && error("output directory already exists: $out")
    runs = length(args) >= 2 ? parse(Int, args[2]) : 4
    nseqs = length(args) >= 3 ? parse(Int, args[3]) : 9000
    beta = length(args) >= 4 ? parse(Float64, args[4]) : 1000.0
    burnin = length(args) >= 5 ? parse(Int, args[5]) : 4000
    thin = length(args) >= 6 ? parse(Int, args[6]) : 2
    local_steps = length(args) >= 7 ? parse(Int, args[7]) : 5
    energy_cache = length(args) >= 8 ? Symbol(args[8]) : :auto
    energy_cache in (:auto, :float, :integer) || error("CACHE must be auto, float, or integer")
    runs >= 4 && nseqs >= 4 && beta > 0.1 && burnin >= 0 && thin > 0 && local_steps > 0 ||
        error("need >=4 independent runs, >=4 samples, beta > 0.1, and valid run lengths")
    betas = [0.0; exp.(range(log(0.1), log(beta); length=25))]
    betas[end] = beta
    seeds = collect(2026103001:2026103000+runs)
    mkpath(out)
    parameters = Dict("structure_set" => "all", "target" => CONTACT_MAP_A,
        "nstructures" => size(STRUCTURES, 2), "betas" => betas, "runs" => runs,
        "samples_per_run" => nseqs, "burnin_rounds" => burnin, "thin_rounds" => thin,
        "local_steps" => local_steps, "seeds" => seeds, "julia_version" => string(VERSION),
        "energy_cache_requested" => string(energy_cache),
        "rng_mode" => "shared", "threaded_replicas" => false,
        "reference_alignment" => Hugo_MSA_path(:A),
        "script_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "sampler_sha256" => bytes2hex(sha256(read(joinpath(@__DIR__, "../../../src/parallel_tempering.jl")))),
        "proposals_per_replica" => (burnin+nseqs*thin)*local_steps,
        "status" => "Unfiltered samples; run diagnose_pt.py before interpreting the figure.")
    open(joinpath(out, "parameters.toml"), "w") do io
        TOML.print(io, parameters)
    end
    # Initialize the shared, immutable contact index before starting threads.
    LatticeProteins.site_contacts(STRUCTURES)
    println("$(runs) ladders, $(length(betas)) replicas each, beta_max=$beta, " *
        "$((burnin+nseqs*thin)*local_steps) proposals per replica")
    flush(stdout)
    Threads.@threads for run in 1:runs
        elapsed = @elapsed r = parallel_tempering(CONTACT_MAP_A; betas, burnin,
            nseqs, thin, local_steps, energy_cache, structures=STRUCTURES, rng=Xoshiro(seeds[run]),
            threaded=false, rng_mode=:shared)
        p = exp.(r.log_pnat[end, :])
        # Independent, direct energy sums verify the cached probabilities.
        checks = unique(round.(Int, range(1, nseqs; length=min(nseqs, 50))))
        errors = [abs(p[i] - exp(log_pnat(CONTACT_MAP_A, r.samples[:, i]))) for i in checks]
        maximum(errors) < 1e-10 || error("cached probability check failed in run $run")
        write_csv(joinpath(out, "run_$(run)_cold.csv"),
            ["sample", "pnat", "replica_id", "sequence"],
            [(i, p[i], r.replica_ids[end,i], aaseq(r.samples[:,i])) for i in 1:nseqs])
        write_csv(joinpath(out, "run_$(run)_trace.csv"),
            ["sample", "beta", "pnat", "replica_id"],
            [(i, betas[j], exp(r.log_pnat[j,i]), r.replica_ids[j,i])
             for i in 1:nseqs for j in eachindex(betas)])
        write_csv(joinpath(out, "run_$(run)_swaps.csv"),
            ["beta_hot", "beta_cold", "attempts", "accepts", "rate"],
            [(betas[j], betas[j+1], r.diagnostics.swap_attempts[j],
              r.diagnostics.swap_accepts[j], r.diagnostics.swap_accepts[j]/r.diagnostics.swap_attempts[j])
             for j in 1:length(betas)-1])
        write_csv(joinpath(out, "run_$(run)_final.csv"), ["beta", "pnat", "sequence"],
            [(betas[j], exp(r.final_log_pnat[j]), aaseq(r.final_seqs[:,j])) for j in eachindex(betas)])
        open(joinpath(out, "run_$(run).msa"), "w") do io
            for i in 1:nseqs
                println(io, beta, ' ', p[i], ' ', aaseq(r.samples[:,i]))
            end
        end
        summary = Dict("run"=>run, "seed"=>seeds[run], "seconds"=>elapsed,
            "energy_cache"=>string(r.energy_cache),
            "mean_pnat"=>mean(p), "minimum_pnat"=>minimum(p), "maximum_pnat"=>maximum(p),
            "fraction_above_099"=>mean(p .> .99), "fraction_above_0995"=>mean(p .> .995),
            "first_half_mean"=>mean(p[1:nseqs÷2]), "last_half_mean"=>mean(p[nseqs÷2+1:end]),
            "burnin_round_trips"=>r.burnin_diagnostics.round_trips,
            "production_round_trips"=>r.diagnostics.round_trips,
            "local_acceptance"=>r.diagnostics.local_accepts ./ r.diagnostics.local_attempts,
            "direct_probability_checks"=>length(checks), "maximum_probability_error"=>maximum(errors))
        open(joinpath(out, "run_$(run)_summary.toml"), "w") do io
            TOML.print(io, summary)
        end
        println("run=$run mean_pnat=$(mean(p)) fraction_above_099=$(mean(p .> .99)) " *
            "round_trips=$(sum(r.diagnostics.round_trips)) seconds=$elapsed")
        flush(stdout)
    end
    open(joinpath(out, "alignment.msa"), "w") do io
        for run in 1:runs
            write(io, read(joinpath(out, "run_$(run).msa")))
        end
    end
    println("Wrote $out; convergence has not yet been assessed.")
end

main(ARGS)
