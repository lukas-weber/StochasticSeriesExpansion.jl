struct Site
    dim::Int
end

struct Bond{NSites}
    type::Int
    sites::NTuple{NSites,Int}
end

struct SSEData{NSites}
    vertex_data::Vector{VertexData{NSites}}

    sites::Vector{Site}
    bonds::Vector{Bond{NSites}}

    energy_offset::Float64
end

function SSEData(
    vertex_data::AbstractVector{VertexData{NSites}},
    sites::Vector{Site},
    bonds::Vector{Bond{NSites}},
) where {NSites}
    energy_offset = sum(b -> vertex_data[b.type].energy_offset, bonds)

    return SSEData(vertex_data, sites, bonds, energy_offset)
end

get_vertex_data(data::SSEData, bond_idx::Integer) =
    data.vertex_data[data.bonds[bond_idx].type]
