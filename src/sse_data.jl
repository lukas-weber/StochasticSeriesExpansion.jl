struct SSESite
    dim::Int
end

"""
    SSEBond(type, (i, j, ...))

Defines a bond for the SSE algorithm. The `type` specifies which `VertexData` object to use. `i`, `j` are two or more site indices that are connected by a bond.
This number is variable but has to be consistent across the model. If the model has bonds with different number of sites, dummy sites can be used to express lower order bonds."""
struct SSEBond{NSites}
    type::Int
    sites::NTuple{NSites,Int}
end

struct SSEData{NSites}
    vertex_data::Vector{VertexData{NSites}}

    sites::Vector{SSESite}
    bonds::Vector{SSEBond{NSites}}

    energy_offset::Float64
end

"""
    SSEData(vertex_data::AbstractVector{<:VertexData}, bonds::AbstractVector{<:SSEBond})

This object holds everything StochasticSeriesExpansion needs to know to simulate a model using the abstract loop algorithm.

The array `vertex_data` contains one instance of [`VertexData`](@ref) for each distinct *type* of bond. The `bonds` define the graph of bonds along with the corresponding bond types.

## Public fields
* `bonds`: the bond information passed on construction
"""
function SSEData(
    vertex_data::AbstractVector{VertexData{NSites}},
    sites::AbstractVector{SSESite},
    bonds::AbstractVector{SSEBond{NSites}},
) where {NSites}
    energy_offset = sum(b -> vertex_data[b.type].energy_offset, bonds)

    return SSEData(vertex_data, sites, bonds, energy_offset)
end

"""
    get_vertex_data(data::SSEData, bond_idx) -> VertexData

The `VertexData` object for a given SSE bond index `bond_idx`.
"""
@inline get_vertex_data(data::SSEData, bond_idx::Integer) =
    @inbounds data.vertex_data[data.bonds[bond_idx].type]
