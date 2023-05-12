using StaticArrays

macro stub(func::Expr)
    return :(
        $func = error(
            "$(func.args[2].args[2]) interface not implemented for type $(typeof(func.args[2].args[1]))",
        )
    )
end

macro stubT(func::Expr)
    return :(
        $func = error(
            "$(only(func.args[2].args[2].args[2].args)) interface not implemented for type $(func.args[2].args[1])",
        )
    )
end

function split_idx(dims::NTuple{D,<:Integer}, idx::Integer) where {D}
    idx -= 1

    r = MVector{D,Int}(undef)
    for (i, d) in enumerate(dims)
        r[i] = idx % d + 1
        idx ÷= d
    end

    return Tuple(r)
end

function join_idx(dims::NTuple{D,<:Integer}, idxs::NTuple{D,<:Integer}) where {D}
    r = 0
    for (idx, d) in zip(Iterators.reverse(idxs), Iterators.reverse(dims))
        r *= d
        r += idx - 1
    end

    return r + 1
end
