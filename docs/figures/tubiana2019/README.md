# Fig. 4 of Tubiana, Cocco & Monasson (2019)

Scripts behind the page "Reproducing Fig. 4 of Tubiana et al. (2019)" of the documentation.
They are not run by the documentation build; the figures in `docs/src/assets/` are
pre-rendered.

The full-structure figure uses the saved parallel-tempering alignment in
`docs/src/assets/tubiana2019_full_pt.msa.gz`: four independent ladders with 79,000
draws each, 316,000 in total. Filtering at `P_nat > 0.99` retains 172,717 sequences.
All 534 monitored observables pass the convergence criteria (worst R-hat 1.00958,
minimum bulk ESS 443.9). Parameters and diagnostics are saved alongside the alignment
in `tubiana2019_full_pt_diagnostics.json`.

Requirements: the package (run the Julia scripts with `--project=<package directory>`),
and Python >=3.11 with the packages in `requirements.txt`:

```bash
python -m pip install -r docs/figures/tubiana2019/requirements.txt
```

To redraw the documented full-structure figure from those samples, run from the
repository root:

```bash
gzip -dc docs/src/assets/tubiana2019_full_pt.msa.gz > full_pt.msa
julia --project=. docs/figures/tubiana2019/analyze_msa.jl full_pt.msa stats_all 0.99
python docs/figures/tubiana2019/plot_fig4.py stats_all docs/src/assets/tubiana2019_fig4_all.jpg "Full-model parallel tempering: monitored convergence checks passed"
```

To generate new alignments:

```bash
# alignment against the subset of 10,000 structures (the model of the paper), 4 chains x 9,000
julia -t 4 --project=. docs/figures/tubiana2019/sample_msa.jl shakhnovich 4 9000 msa_subset.txt
# statistics of the reference-subset sample
julia --project=. docs/figures/tubiana2019/analyze_msa.jl msa_subset.txt stats_subset
python docs/figures/tubiana2019/plot_fig4.py stats_subset fig4_subset.png "subset of 10,000 structures"

# Full set: four independent parallel-tempering ladders, 26 replicas each,
# beta_max=1000, 4000 burn-in rounds, initially 9000 draws per ladder, thin=2,
# five local mutations per replica per round. Use NEW output directories.
# The documented run began with the original floating-point kernel.
julia -t 4 --project=. docs/figures/tubiana2019/sample_pt.jl full_pt_initial 4 9000 1000 4000 2 5 float
# Continue with exact integer caching until all monitored gates pass.
# Set TUBIANA_PYTHON if the required Python executable is not named python.
julia -t 4 --project=. docs/figures/tubiana2019/continue_pt.jl full_pt_initial full_pt 2000 5
# Inspect diagnostics.json/observables.csv and diagnostics.png before interpreting.
julia --project=. docs/figures/tubiana2019/analyze_msa.jl full_pt/analysis_alignment.msa stats_all 0.99
python docs/figures/tubiana2019/plot_fig4.py stats_all fig4_all.png "parallel tempering against all 103,346 structures"

# agreement between two alignments
python docs/figures/tubiana2019/compare_stats.py stats_subset stats_all comparison.png
```

`sample_pt.jl OUT [RUNS NSEQS BETA BURNIN THIN LOCAL_STEPS CACHE]` records all inputs,
independent RNG seeds, code hashes, cold sequences, temperature traces, exchange
acceptance, and round trips. It checks 50 saved probabilities per run by recomputing
the complete energy sum. Increasing `NSEQS` or changing the thinning interval does
not by itself establish equilibrium. `CACHE` defaults to `auto`; use `float` to
reproduce the initial phase above or `integer` to require the accelerated kernel.

`continue_pt.jl SOURCE OUT [BATCH_SAMPLES=2000] [CHECK_EVERY=5]` preserves all initial
draws and extends the same four ladders with the same temperatures and proposal
schedule, with no additional burn-in or post-hoc discards. It checks convergence
after every five batches by default and stops only when the gates below pass.
It can take a long time: the saved sample count greatly exceeds the effective
sample size for slow sequence observables.

The initial driver records final states but not RNG states, so the continuation
starts with documented fresh independent RNG seeds. Thereafter, atomic per-ladder
checkpoints preserve both RNGs and temperature states. Repeating the same command
resumes the last completed chunks. Keep the Julia version, sampler source, source
directory, and batch settings fixed when resuming. Chunks end on 100-round
boundaries so cache refreshes and alternating swap phases match an uninterrupted
call. Endpoint accounting requires the even ladder and `thin=2` used here.
The output retains cold draws, temperature traces, direct probability checks, and
each convergence report. Original and continuation phases have separate code
hashes and seeds in `parameters.toml`.

Both figure drivers explicitly use `threaded=false, rng_mode=:shared` inside each
ladder and parallelize the four ladders externally. This preserves the single-stream
protocol used to generate the saved alignment, including across chunk boundaries.
The package's default `rng_mode=:per_temperature` uses the newer independent
temperature streams; disabling threads alone does not select the original protocol.
The archived report retains the original run's code hashes.

`diagnose_pt.py` checks rank-normalized split R-hat and bulk effective sample size
using ArviZ. It monitors folding probability, threshold membership, five leading
PCs, the reference PC1, and residue indicators with frequencies between 1% and 99%.
It also checks threshold-weighted indicators/PCs on the original time grid, and
reports Monte Carlo errors for conditional residue frequencies. The default gates
are R-hat <1.01, bulk ESS >=400, and a production round trip in every ladder.
Passing these gates is evidence for the monitored observables, not proof of global
equilibrium or a guarantee about rare residues. The explicit `--discard` option
can inspect a trajectory's tail; full-phase round-trip counts cannot certify that
tail, so it does not automatically receive a passing status.

The sampler's correctness is tested separately against direct energy sums, swap
detailed balance, an exactly enumerable three-residue model, and uninterrupted
versus checkpointed trajectories in
`test/parallel_tempering.jl`. Check the Python diagnostics with:

```bash
python -m unittest discover -s docs/figures/tubiana2019 -p 'test_diagnose_pt.py'
```

`analyze_msa.jl` records its filter in `analysis.toml`. `plot_fig4.py` reads that
metadata and rejects conflicting threshold labels. Panel (c) uses a seeded random
selection across the retained alignment rather than the first sequences of one chain.

The reference alignment of Jacquin et al. (2016) can be analyzed in the same way, its path
is `LatticeProteins.Hugo_MSA_path(:A)`.
