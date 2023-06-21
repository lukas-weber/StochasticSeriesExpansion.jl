
const OperCodeUInt = UInt64
const vertex_code_maxbits = 8 * 3 + 1

abstract type Identity end

"""
    VertexCode

represents a ketbra or *vertex* or *nonbranching operator* in the abstract loop formulation of the stochastic series expansion. 
"""
struct VertexCode
    code::OperCodeUInt
end

VertexCode(::Nothing) = VertexCode((1 << vertex_code_maxbits) + 1)

function VertexCode(diagonal::Bool, vertex_idx::Integer)
    v = VertexCode(diagonal | (vertex_idx << 1))
    return v
end

isdiagonal(v::VertexCode) = Bool(v.code & 1)
isinvalid(v::VertexCode) = v.code >= 1 << vertex_code_maxbits
function get_vertex_idx(v::VertexCode)
    return v.code >> 1
end

Base.show(io::IO, v::VertexCode) =
    print(io, isinvalid(v) ? "VertexCode(invalid)" : "VertexCode($(get_vertex_idx(v)))")

"""
    OperCode

Represents an operator in the operator string.
"""
struct OperCode
    code::OperCodeUInt
end

OperCode(::Type{Identity}) = OperCode(0)

function OperCode(bond::Integer, vertex::VertexCode)
    v = vertex.code

    return OperCode(1 | v << 1 | bond << (1 + vertex_code_maxbits))
end

"""
    get_bond(op::OperCode) -> Int

Get the SSE bond index of a specific operator. This is the same index as used by the `bonds` field of [`SSEData`](@ref).
"""
@inline get_bond(o::OperCode)::Int = o.code >> (1 + vertex_code_maxbits)

"""
    get_vertex(op::OperCode) -> VertexCode

Get the [`VertexCode`](@ref) of the corresponding opercode. This can be used in conjunction with [`VertexData`](@ref) to get the leg states (or matrix elements) of that operator.
"""
@inline get_vertex(o::OperCode) =
    VertexCode((o.code & ((1 << vertex_code_maxbits) - 1)) >> 1)

@inline isidentity(o::OperCode) = o.code == 0

"""
    isdiagonal(op::OperCode) -> Bool

Returns true if operator is diagonal.
"""
@inline isdiagonal(o::OperCode) = isdiagonal(get_vertex(o))
