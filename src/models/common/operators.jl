using LinearAlgebra

"""Returns the tuple (S+,Sz) of spin operators for `S=(dimension-1)/2`"""
function spin_operators(::Type{T}, dimension::Integer) where {T<:AbstractFloat}
    S = T(dimension - 1) / 2

    sz = zeros(T, dimension, dimension)
    splus = zeros(T, dimension, dimension)

    for i = 1:dimension
        m = S - i + 1
        sz[i, i] = m

        if i < dimension
            splus[i, i+1] = sqrt((S - m + 1) * (S + m))
        end
    end

    return (splus, sz)
end

spin_operators(dimension::Integer) = spin_operators(Float64, dimension)

"""Returns the bosonic creation operator of `dimension` which corresponds to a boson number cutoff of `n < dimension`"""
function bosonic_creation_operator(::Type{T}, dimension) where {T<:AbstractFloat}
    return diagm(1 => sqrt.(T.(1:dimension)))
end

bosonic_creation_operator(dimension::Integer) =
    bosonic_creation_operator(Float64, dimension)
