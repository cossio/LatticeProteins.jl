# Continue the four independent ladders until the unchanged diagnostic gates pass.
# julia -t 4 --project=. continue_pt.jl SOURCE OUT [BATCH_SAMPLES=2000] [CHECK_EVERY=5]
# Set TUBIANA_PYTHON to the Python executable with requirements.txt installed.
# Checkpoints preserve the RNG stream and temperature states. Re-running the same
# command resumes after the last committed chunk, without dropping saved draws.
using LatticeProteins: LatticeProteins, parallel_tempering, STRUCTURES, CONTACT_MAP_A,
    potts, aaseq, log_pnat
using Random: Xoshiro
using DelimitedFiles: readdlm, writedlm
using Serialization: serialize, deserialize
using TOML
using SHA: sha256
using Statistics: mean
const LP = LatticeProteins

function csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ',')); writedlm(io, rows, ',')
    end
end
function checkpoint(path, state)
    temporary = path * ".tmp"
    serialize(temporary, state)
    mv(temporary, path; force=true)
end
function initial_state(source, run, parameters)
    final = readdlm(joinpath(source, "run_$(run)_final.csv"), ',', Any; skipstart=1)
    @assert Float64.(final[:,1]) == parameters["betas"]
    ids = zeros(Int, length(parameters["betas"]))
    for line in Iterators.drop(eachline(joinpath(source, "run_$(run)_trace.csv")), 1)
        fields = split(line, ',')
        if parse(Int, fields[1]) == parameters["samples_per_run"]
            i = findfirst(==(parse(Float64,fields[2])), parameters["betas"])
            ids[i] = parse(Int, fields[4])
        end
    end
    @assert sort(ids) == collect(eachindex(ids))
    stages = zeros(Int, length(ids)); stages[ids[1]] = 1
    seed = 2026103100 + run
    Dict{Symbol,Any}(:seqs=>hcat([potts(String(s)) for s in final[:,3]]...),
        :rng=>Xoshiro(seed), :ids=>ids, :stages=>stages, :trips=>zeros(Int,length(ids)),
        :batch=>0, :added_samples=>0, :seed=>seed, :seconds=>0.0,
        :max_error=>0.0, :checks=>0,
        :local_attempts=>zeros(Int,length(ids)), :local_accepts=>zeros(Int,length(ids)),
        :swap_attempts=>zeros(Int,length(ids)-1), :swap_accepts=>zeros(Int,length(ids)-1))
end
function advance!(state, parameters, out, run, batch_samples)
    betas = parameters["betas"]
    elapsed = @elapsed result = parallel_tempering(CONTACT_MAP_A; betas,
        burnin=0, nseqs=batch_samples, thin=parameters["thin_rounds"],
        local_steps=parameters["local_steps"], rng=state[:rng], initial_seqs=state[:seqs],
        structures=STRUCTURES, energy_cache=:integer)
    global_ids = state[:ids][result.replica_ids]
    # With an even ladder and thin=2, both endpoint edges are attempted on odd
    # rounds. The saved even-round snapshots therefore capture every endpoint
    # change, including across checkpoint boundaries.
    for ids in eachcol(global_ids)
        hot, cold = ids[1], ids[end]
        state[:stages][hot] == 2 && (state[:trips][hot] += 1)
        state[:stages][hot] = 1
        state[:stages][cold] == 1 && (state[:stages][cold] = 2)
    end
    offset = parameters["initial_samples_per_run"] + state[:added_samples]
    p = exp.(result.log_pnat[end,:])
    indices = unique(round.(Int, range(1,batch_samples; length=min(50,batch_samples))))
    errors = [abs(p[i]-exp(log_pnat(CONTACT_MAP_A,result.samples[:,i]))) for i in indices]
    maximum(errors) < 1e-10 || error("direct probability check failed")
    batch = state[:batch]+1
    prefix = joinpath(out,"chunks","run_$(run)_batch_$(batch)")
    csv(prefix*"_cold.csv", ["sample","pnat","replica_id","sequence"],
        [(offset+i,p[i],global_ids[end,i],aaseq(result.samples[:,i])) for i in 1:batch_samples])
    csv(prefix*"_trace.csv", ["sample","beta","pnat","replica_id"],
        [(offset+i,betas[j],exp(result.log_pnat[j,i]),global_ids[j,i])
         for i in 1:batch_samples for j in eachindex(betas)])
    state[:seqs] = result.final_seqs
    state[:ids] = global_ids[:,end]
    state[:batch] = batch
    state[:added_samples] += batch_samples
    state[:seconds] += elapsed
    state[:max_error] = max(state[:max_error],maximum(errors))
    state[:checks] += length(errors)
    for field in (:local_attempts,:local_accepts,:swap_attempts,:swap_accepts)
        state[field] .+= getproperty(result.diagnostics,field)
    end
    checkpoint(joinpath(out,"checkpoint_$(run).jls"),state)
    println("run=$run batch=$batch total_samples=$(offset+batch_samples) continuation_round_trips=$(sum(state[:trips])) seconds=$(round(elapsed;digits=1))")
    flush(stdout)
end
function aggregate(source,out,parameters,states)
    for (run,state) in enumerate(states)
        coldpath = joinpath(out,"run_$(run)_cold.csv")
        open(coldpath,"w") do io
            write(io,read(joinpath(source,"run_$(run)_cold.csv")))
            for batch in 1:state[:batch]
                path = joinpath(out,"chunks","run_$(run)_batch_$(batch)_cold.csv")
                for line in Iterators.drop(eachline(path),1)
                    println(io,line)
                end
            end
        end
        p = Float64[]
        open(joinpath(out,"run_$(run).msa"),"w") do io
            for line in Iterators.drop(eachline(coldpath),1)
                fields=split(line,','); value=parse(Float64,fields[2]); push!(p,value)
                println(io,last(parameters["betas"]),' ',value,' ',fields[4])
            end
        end
        summary = TOML.parsefile(joinpath(source,"run_$(run)_summary.toml"))
        summary["production_round_trips"] .+= state[:trips]
        summary["continuation_round_trips"] = state[:trips]
        summary["seconds"] += state[:seconds]
        summary["maximum_probability_error"] = max(summary["maximum_probability_error"],state[:max_error])
        summary["direct_probability_checks"] += state[:checks]
        summary["mean_pnat"] = mean(p)
        summary["minimum_pnat"] = minimum(p); summary["maximum_pnat"] = maximum(p)
        summary["first_half_mean"] = mean(p[1:length(p)÷2])
        summary["last_half_mean"] = mean(p[length(p)÷2+1:end])
        summary["fraction_above_099"] = mean(p .> .99)
        summary["fraction_above_0995"] = mean(p .> .995)
        summary["continuation_local_acceptance"] = state[:local_accepts] ./ state[:local_attempts]
        summary["continuation_seed"] = state[:seed]
        open(io->TOML.print(io,summary),joinpath(out,"run_$(run)_summary.toml"),"w")
        csv(joinpath(out,"run_$(run)_final.csv"),["beta","pnat","sequence"],
            [(parameters["betas"][j],exp(log_pnat(CONTACT_MAP_A,state[:seqs][:,j])),aaseq(state[:seqs][:,j]))
             for j in eachindex(parameters["betas"])])
    end
    open(joinpath(out,"alignment.msa"),"w") do io
        for run in eachindex(states);write(io,read(joinpath(out,"run_$(run).msa")));end
    end
    parameters["samples_per_run"] = parameters["initial_samples_per_run"]+states[1][:added_samples]
    parameters["proposals_per_replica"] = parameters["initial_proposals_per_replica"] +
        states[1][:added_samples]*parameters["thin_rounds"]*parameters["local_steps"]
    open(io->TOML.print(io,parameters),joinpath(out,"parameters.toml"),"w")
end
function main(args)
    length(args)>=2 || error("provide SOURCE OUT [BATCH_SAMPLES] [CHECK_EVERY]")
    source,out=abspath.(args[1:2])
    source!=out || error("source and output must differ")
    batch_samples=length(args)>=3 ? parse(Int,args[3]) : 2000
    check_every=length(args)>=4 ? parse(Int,args[4]) : 5
    batch_samples>0 && check_every>0 || error("invalid batch/check interval")
    parameters=TOML.parsefile(joinpath(source,"parameters.toml"))
    parameters["thin_rounds"]==2 && iseven(length(parameters["betas"])) ||
        error("endpoint accounting requires an even ladder and thin=2")
    batch_samples*parameters["thin_rounds"] % 100 == 0 || error("chunks must end on 100-round boundaries")
    root=normpath(joinpath(@__DIR__,"../../.."))
    sampler_hash=bytes2hex(sha256(read(joinpath(root,"src/parallel_tempering.jl"))))
    if isdir(out)
        saved=TOML.parsefile(joinpath(out,"parameters.toml"))
        saved["continuation_sampler_sha256"]==sampler_hash || error("sampler changed since checkpoint")
        saved["source_directory"]==source || error("checkpoint source differs")
        parameters=saved
        states=[deserialize(joinpath(out,"checkpoint_$(r).jls")) for r in 1:parameters["runs"]]
    else
        mkpath(joinpath(out,"chunks"))
        states=[initial_state(source,r,parameters) for r in 1:parameters["runs"]]
        parameters["initial_samples_per_run"]=parameters["samples_per_run"]
        parameters["initial_proposals_per_replica"]=parameters["proposals_per_replica"]
        parameters["source_directory"]=source
        parameters["continuation_sampler_sha256"]=sampler_hash
        parameters["continuation_script_sha256"]=bytes2hex(sha256(read(@__FILE__)))
        parameters["continuation_energy_cache"]="integer"
        parameters["continuation_seeds"]=[state[:seed] for state in states]
        parameters["round_trip_accounting"]="Sum of complete trips in original and continued production phases; checkpoint boundaries preserve continued-phase endpoint history."
        for (r,state) in enumerate(states);checkpoint(joinpath(out,"checkpoint_$(r).jls"),state);end
        open(io->TOML.print(io,parameters),joinpath(out,"parameters.toml"),"w")
    end
    LP.site_contacts(STRUCTURES)
    python=get(ENV,"TUBIANA_PYTHON","python")
    while true
        batch=minimum(state[:batch] for state in states)+1
        Threads.@threads for r in eachindex(states)
            states[r][:batch]>=batch && continue
            advance!(states[r],parameters,out,r,batch_samples)
        end
        if batch % check_every == 0
            aggregate(source,out,parameters,states)
            script=joinpath(@__DIR__,"diagnose_pt.py")
            run(`$python $script $out --threshold 0.99`)
            cp(joinpath(out,"diagnostics.json"),joinpath(out,"diagnostics_batch_$(batch).json");force=true)
            # Python reads JSON and returns success only on the unchanged gates.
            code="import json,sys;sys.exit(0 if json.load(open(sys.argv[1]))['diagnostic_gates_passed'] else 1)"
            if success(`$python -c $code $(joinpath(out,"diagnostics.json"))`)
                println("CONVERGENCE GATES PASSED after $(parameters["samples_per_run"]) samples per ladder.")
                break
            end
            println("Continuing: convergence gates not met.");flush(stdout)
        end
    end
end
main(ARGS)
