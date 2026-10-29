# Fig. 4 of Tubiana, Cocco & Monasson (2019)

Scripts behind the page "Reproducing Fig. 4 of Tubiana et al. (2019)" of the documentation.
They are not run by the documentation build; the figures in `docs/src/assets/` are
pre-rendered.

Requirements: the package (run the Julia scripts with `--project=<package directory>`),
Python 3 with `numpy`, `pandas`, `matplotlib` and `logomaker` for the plots.

```bash
# alignment against the subset of 10,000 structures (the model of the paper), 4 chains x 9,000
julia -t 4 --project=. docs/figures/tubiana2019/sample_msa.jl shakhnovich 4 9000 msa_subset.txt
# alignment against all 103,346 compact structures, 20 chains x 100 (optional trailing
# arguments: β, default 1000, and burn-in steps, default 10,000)
julia -t 4 --project=. docs/figures/tubiana2019/sample_msa.jl all 20 100 msa_all.txt
# statistics of panels (c)-(g), then the figures
julia --project=. docs/figures/tubiana2019/analyze_msa.jl msa_subset.txt stats_subset
julia --project=. docs/figures/tubiana2019/analyze_msa.jl msa_all.txt stats_all
python docs/figures/tubiana2019/plot_fig4.py stats_subset fig4_subset.png "subset of 10,000 structures"
python docs/figures/tubiana2019/plot_fig4.py stats_all fig4_all.png "all 103,346 compact structures"
# agreement between two alignments
python docs/figures/tubiana2019/compare_stats.py stats_subset stats_all comparison.png
```

The reference alignment of Jacquin et al. (2016) can be analyzed in the same way, its path
is `LatticeProteins.Hugo_MSA_path(:A)`.
