using RDKitMinimalLib
using GeometricFlux
using JSON
using GeometricFlux
using DataFrames
using LightGraphs

function process_atom_data(lines)
    matx = zeros(Float64, 4, length(lines))
    @show lines
    for (i,l) in lines
        @show i,l
        coord_strs = strip(l[2:end]) |> x -> split(x, '\t')
        coord_strs = coord_strs[coord_strs .!= ""]
        coords = map(x -> parse(Float64, x), coord_strs)
        @assert length(coords) == 4
        matx[:,i] = coords  
    end
    return matx
end

function get_df(filename)
    raw_lines = open(filename) do io
        readlines(io, keep=false)
    end
    n_atoms = parse(Int, raw_lines[1])


    props_meta = ["tag+i", "A", "B", "C", "mu", "alpha", "eps_homo",
                  "eps_lumo", "eps_gap", "R2", "zpve", "U0", "U",
                  "H", "G", "Cv"]
    raw_props = split(raw_lines[2], '\t')
    props_dict = Dict(zip(props_meta, raw_props))
    raw_atom_data = raw_lines[3:3+n_atoms]
    smiles = split(raw_lines[4+n_atoms], '\t')

    mol_data = RDKitMinimalLib.get_mol(smiles[1]) |> RDKitMinimalLib.get_json |> JSON.parse
    moljson = mol_data["molecules"][1]

    mol_graph = SimpleGraph(n_atoms)

    num_bonds = length(moljson["bonds"])
    nf = process_atom_data(raw_atom_data)
    ef = zeros(1, num_bonds)
    for (i,bond) in moljson["bonds"]
        add_edge!(mol_graph, bond["atoms"][1] + 1,  bond["atoms"][2] + 1)
        ef[:,i] = get(bond, "bo", 1)
    end

    mol_fg = FeaturedGraph(mol_graph, nf=nf, ef=ef)

    @show data

    return Dict(
                props_dict...,
                "molecule" => mol_fg
           )
end
