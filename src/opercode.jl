
const OperCodeUInt = UInt64
const vertex_code_maxbits = 8 * 3 + 1

abstract type Identity end

struct VertexCode
    code::OperCodeUInt
end

function VertexCode(diagonal::Bool, vertex_idx::Integer)
    v = VertexCode(diagonal | (vertex_idx << 1))
    @debug begin
        @assert !isinvalid(v)
    end
    return v
end

@inline isdiagonal(v::VertexCode) = v.code & 1
@inline get_vertex_idx(v::VertexCode) = v.code >> 1
@inline isinvalid(v::VertexCode) = v.code >= 1 << vertex_code_maxbits

struct OperCode
    code::OperCodeUInt
end

OperCode(::Type{Identity}) = OperCode(0)

function OperCode(bond::Integer, vertex::VertexCode)
    v = vertex.code

    @debug begin
        @assert !isinvalid(v)
        @assert bond < (1 << (8 * sizeof(OperCodeUInt) - vertex_code_maxbits - 1))
    end
    return OperCode(1 | v << 1 | bond << (1 + vertex_code_maxbits))
end

@inline get_bond(o::OperCode) = o.code >> (1 + vertex_code_maxbits)
@inline get_vertex(o::OperCode) =
    VertexCode((o.code & ((1 << vertex_code_maxbits) - 1)) >> 1)

@inline isidentity(o::OperCode) = o.code == 0
@inline isdiagonal(o::OperCode) = isdiagonal(get_vertex(o))
