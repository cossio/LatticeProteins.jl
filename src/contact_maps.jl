function load_contact_maps()
    cm = readdlm(joinpath(artifact"ContactMaps", "contact_maps_10000.dat"), Int)[:, 2:end]
    return reshape(cm, 28, 10_000, 2) # array[:, s, :] gives contacts of structure 's'
end

const CONTACT_MAPS = load_contact_maps()
const N_STRUCTURES = size(CONTACT_MAPS, 2)
