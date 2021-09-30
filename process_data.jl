using RDKitMinimalLib
using GeometricFlux
using JSON
using GeometricFlux

r"*^"

function get_df(filename)
    raw_string = open(filename) do io
        read(io, String)
    end
    raw_lines = split(raw_string, '\n')
    n_atoms = Int(raw_lines[1])


    props_meta = ["tag+i", "A", "B", "C", "mu", "alpha", "eps_homo",
                  "eps_lumo", "eps_gap", "R2", "zpve", "U0", "U",
                  "H", "G", "Cv"]
    raw_props = split(raw_lines[2], '\t')
    props_dict = Dict(zip(props_meta, raw_props))
    raw_atom_data = raw_lines[3:3+n_atoms]
    smiles = split(raw_lines[4+n_atoms], '\t')

    mol_data = RDKitMinimalLib.get_mol(smiles[1]) |> RDKitMinimalLib.get_json

    return Dict(
                props_dict...
           )
end
