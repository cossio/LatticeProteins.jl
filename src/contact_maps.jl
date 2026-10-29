"""
    L

Sequence length for lattice proteins on the 3×3×3 cubic lattice (27).
"""
const L = 27 # sequence length

"""
    load_contact_maps()

Load the pre-computed contact maps of the representative subset of 10,000 compact
structures of Heo, Maslov & Shakhnovich (2011). Returns an `Array{Int, 3}` of size
`(28, 10000, 2)`, where `cm[:, s, :]` gives the 28 pairs of contacting residue indices of
structure `s`. See [`SHAKHNOVICH_STRUCTURES`](@ref).
"""
function load_contact_maps()
    cm = readdlm(joinpath(artifact"ContactMaps", "contact_maps_10000.dat"), Int)[:, 2:end]
    return reshape(cm, 28, 10_000, 2) # array[:, s, :] gives contacts of structure 's'
end

"""
    SHAKHNOVICH_STRUCTURES

Contact maps of the representative subset of 10,000 compact structures selected by the
Shakhnovich group (Heo, Maslov & Shakhnovich, 2011) out of the 103,346 of
[`STRUCTURES`](@ref), and used as the set of competing structures by Jacquin et al. (2016)
and Tubiana et al. (2019). An `Array{Int, 3}` of size `(28, 10000, 2)`; use
`SHAKHNOVICH_STRUCTURES[:, s, :]` to get the contact pairs of structure `s`.
See [`load_contact_maps`](@ref).

Pass `structures = SHAKHNOVICH_STRUCTURES` to [`pnat`](@ref), [`metropolis!`](@ref) and
related functions to compete only against this set, as in those papers and in the
[`Hugo_MSA_pnat`](@ref) values. Structure indices refer to the set in use; use
[`structure_index`](@ref) to translate between the two sets. For example,
`structure_index(SHAKHNOVICH_STRUCTURES, CONTACT_MAP_A)` gives the index (625) of reference
structure A in this set.

# References

- Heo, M., Maslov, S., and Shakhnovich, E. "Topology of protein interaction network shapes
  protein abundances and strengths of their functional and nonspecific interactions."
  *Proceedings of the National Academy of Sciences* 108.10 (2011): 4258-4263.
- Jacquin, H., et al. "Benchmarking inverse statistical approaches for protein structure
  and design with exactly solvable models." *PLoS Computational Biology* 12.5 (2016):
  e1004889.
"""
const SHAKHNOVICH_STRUCTURES = load_contact_maps()
