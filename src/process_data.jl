using RDKitMinimalLib
using GeometricFlux
using JSON
using GeometricFlux
using DataFrames
using LightGraphs
using Serialization

function process_atom_data(lines)
    matx = zeros(Float64, 4, length(lines))
    for (i,l) in enumerate(lines)
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
    raw_atom_data = raw_lines[3:3+n_atoms-1]
    raw_numbers = raw_lines[3+n_atoms]
    smiles = split(raw_lines[4+n_atoms], '\t')

    mol_data = RDKitMinimalLib.get_mol(smiles[1]) |> RDKitMinimalLib.get_json |> JSON.parse
    moljson = mol_data["molecules"][1]

    mol_graph = SimpleGraph(n_atoms)

    num_bonds = length(moljson["bonds"])
    num_non_hydrogen = length(moljson["atoms"])

    for (i,bond) in enumerate(moljson["bonds"])
        add_edge!(mol_graph, bond["atoms"][1] + 1,  bond["atoms"][2] + 1)
    end

    Hcount = 0
    for (i,atom) in enumerate(moljson["atoms"])
        n_implicit_H = get(atom, "impHs", 0)
        for _ in 1:n_implicit_H
            add_edge!(mol_graph, i, num_non_hydrogen + 1 + Hcount)
            Hcount += 1
            num_bonds += 1
        end
    end
    @assert Hcount == n_atoms - num_non_hydrogen

    nf = process_atom_data(raw_atom_data)
    ef = zeros(1, num_bonds)
    bond_counter = 1
    for bond in moljson["bonds"]
        ef[bond_counter] = get(bond, "bo", 1)
        bond_counter += 1
    end

    for atom in moljson["atoms"]
        n_implicit_H = get(atom, "impHs", 0)
        for _ in 1:n_implicit_H
            ef[bond_counter] = 1
            bond_counter += 1
        end
    end
    mol_fg = FeaturedGraph(mol_graph, nf=nf, ef=ef)


    return Dict(
                props_dict...,
                "molecule" => mol_fg,
                "rawline" => raw_numbers,
           )
end


datarows = []

DATADIR = "chemdata/"
for datafile in readdir(DATADIR)
    try
       push!(datarows, get_df(DATADIR * datafile))
    catch
        println("Error in file $(datafile)")
    end
end

serialize("data.blob", datarows)
