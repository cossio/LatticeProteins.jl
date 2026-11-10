"""
    sample_msa(s, method; nseqs = 100, rng = default_rng())

Sample a multiple sequence alignment (MSA) of `nseqs` sequences that fold into structure
`s`: sequences drawn with probability proportional to `pnat(s, seq)^β`, where the inverse
temperature `β` sets the strength of the selection for folding.

`method` chooses the sampler, and holds its parameters:

- [`Metropolis`](@ref)`(; β, nsteps)`: independent Metropolis chains, one per sequence.
  Simple, and adequate at moderate selection.
- [`ParallelTempering`](@ref)`(; betas, ...)`: replica exchange over a ladder of inverse
  temperatures. For strong selection, where a single chain does not equilibrate.

Returns a named tuple. Its field `samples` is the alignment, a `Matrix{Int}` of size
`(27, nseqs)` whose columns are integer-coded sequences (see [`aaseq`](@ref)). The other
fields depend on the method, and are described with it.

```julia
msa = sample_msa(STRUCTURE_A, Metropolis(β = 100); nseqs = 20).samples
aaseq(msa)
```
"""
function sample_msa end
