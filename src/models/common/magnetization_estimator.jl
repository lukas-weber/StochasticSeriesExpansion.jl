@doc raw"""
    MagnetizationEstimator{
        OrderingVector,
        StaggerUC,
        Model,
        Prefix} <: AbstractOpstringEstimator

Generic operator string estimator that can be used for all models that have

* a field `lattice::Lattice`
* implement [`magnetization_state`](@ref).

optionally, [`magnetization_lattice_site_idx`](@ref).

It computes the
* `:Mag`, magnetization ``\langle m\rangle``
* `:AbsMag`, absolute magnetization ``\langle |m|\rangle``
* `:Mag2`, squared magnetization ``\langle m^2\rangle``
* `:Mag4`, quartic magnetization ``\langle m^4\rangle``
* `:MagChi`, susceptibility ``N \int_0^\beta d\tau \langle m(\tau) m\rangle``
* `:BinderRatio`, Binder ratio ``\langle m^2\rangle^2/\langle m^4\rangle``

Here, ``m = \frac{1}{N} \sum_i m_i``, the ``m_i`` are given by [`magnetization_state`](@ref) and ``N`` is returned by [`normalization_site_count`](@ref).

## Type parameters

* `OrderingVector`: Fourier component to compute in units of π. For example ``(1,1)`` corresponds to ``(π,π)``
* If `StaggerUC` is true, the sublattice sign of the unitcell is additionally taken into account.
* `Model`: model type to apply this estimator to.
* `Prefix`: (Symbol) this prefix is added to all observable names. Consider using [`magnetization_estimator_standard_prefix`](@ref).
"""
Base.@kwdef mutable struct MagnetizationEstimator{
    OrderingVector,
    StaggerUC,
    Model<:AbstractModel,
    Prefix,
} <: AbstractOpstringEstimator
    model::Model

    n::Float64

    tmpmag::Float64
    mag::Float64
    absmag::Float64
    mag2::Float64
    mag4::Float64
end

"""
    all_magnetization_estimators(model::AbstractModel, dimension)

Generates a list of all valid `MagnetizationEstimator` types for a given model. Useful if you want to measure them all.
"""
function all_magnetization_estimators(model::Type{<:AbstractModel}, dimension::Integer)
    return [
        MagnetizationEstimator{
            q,
            stagger_uc,
            model,
            magnetization_estimator_standard_prefix(q, stagger_uc),
        } for stagger_uc in (false, true),
        q in Iterators.product(((false, true) for _ = 1:dimension)...)
    ]
end

"""
    magnetization_state(model::AbstractModel, site_idx, state_idx)

This interface needs to be implemented by any model that wants to use `MagnetizationEstimator`. 
"""
function magnetization_state end

"""
    magnetization_lattice_site_idx(model::AbstractModel, sse_site_idx) -> Integer

In models where additional degrees of freedom exist, this function maps sse site indices to physical lattice site indices.
"""
magnetization_lattice_site_idx(::AbstractModel, sse_site_idx::Integer) = sse_site_idx

function staggered_sign(
    est::MagnetizationEstimator{OrderingVector,StaggerUC,Model,Prefix},
    site_idx::Integer,
)::Int where {OrderingVector,StaggerUC,Model,Prefix}
    return staggered_sign(est.model.lattice, OrderingVector, StaggerUC, site_idx)
end

function init(
    ::Type{MagnetizationEstimator{OrderingVector,StaggerUC,Model,Prefix}},
    model::AbstractModel,
    state::AbstractVector{StateIndex},
) where {OrderingVector,StaggerUC,Model,Prefix}
    tmpmag = @inline @fastmath sum(
        site ->
            staggered_sign(model.lattice, OrderingVector, StaggerUC, site) *
            magnetization_state(model, site, state[site]),
        1:site_count(model.lattice),
    )

    mag = tmpmag
    absmag = abs(tmpmag)
    mag2 = tmpmag^2
    mag4 = tmpmag^4

    # use typeof(model) here because Model may be an abstract type
    return MagnetizationEstimator{OrderingVector,StaggerUC,typeof(model),Prefix}(;
        model,
        n = 1,
        tmpmag,
        mag,
        absmag,
        mag2,
        mag4,
    )
end

@inline function measure(
    est::MagnetizationEstimator,
    op::OperCode,
    ::AbstractVector{StateIndex},
    sse_data::SSEData,
)
    @inbounds if !isdiagonal(op)
        bond = sse_data.bonds[get_bond(op)]
        vd = get_vertex_data(sse_data, get_bond(op))
        leg_state = get_leg_state(vd, get_vertex(op))

        nsites = length(vd.dims)
        for l = 1:nsites
            site = magnetization_lattice_site_idx(est.model, bond.sites[l])
            if site !== nothing
                est.tmpmag +=
                    staggered_sign(est, site) * (
                        magnetization_state(est.model, site, leg_state[nsites+l]) -
                        magnetization_state(est.model, site, leg_state[l])
                    )
            end
        end
    end

    est.mag += est.tmpmag
    est.absmag += abs(est.tmpmag)
    tmpmag2 = est.tmpmag * est.tmpmag
    est.mag2 += tmpmag2
    est.mag4 += tmpmag2 * tmpmag2

    est.n += 1

    return nothing
end

"""
    magnetization_estimator_standard_prefix(q::Tuple{Bool...}, stagger_uc::Bool)

Returns the conventional prefix for a magnetization estimator with ordering vector `q` and `stagger_uc` set.
"""
function magnetization_estimator_standard_prefix(q::Tuple{Vararg{Bool}}, stagger_uc::Bool)
    if !any(q) && !stagger_uc
        return Symbol()
    end
    direction_name(i) = mod('X' - 'A' + i - 1, 26) + 'A'
    return Symbol(
        "Stag" * join("$(direction_name(i))" for i in findall(q)) * "uc"^stagger_uc,
    )
end

"""
    magnetization_estimator_obs_symbols(prefix::Symbol)
    magnetization_estimator_obs_symbols(est::MagnetizationEstimator)

Generates the symbols used for observable names by the `MagnetizationEstimator` based on `prefix`. Two named tuples are returned for the normal and sign-multiplied observables respectively.
"""
magnetization_estimator_obs_symbols(est::MagnetizationEstimator) =
    magnetization_estimator_obs_symbols(Val(get_prefix(typeof(est))))
magnetization_estimator_obs_symbols(prefix::Symbol) =
    magnetization_estimator_obs_symbols(Val(prefix))
@generated function magnetization_estimator_obs_symbols(::Val{prefix}) where {prefix}
    obsnames = ["Mag", "AbsMag", "Mag2", "Mag4", "MagChi", "BinderRatio"]
    obs = [Symbol.((lowercase(name), String(prefix) * name)) for name in obsnames]
    obssigned =
        [Symbol.((lowercase(name), "Sign" * String(prefix) * name)) for name in obsnames]

    return quote
        return ($([:($name = $(QuoteNode(prefixname))) for (name, prefixname) in obs]...),),
        ($([:($name = $(QuoteNode(prefixname))) for (name, prefixname) in obssigned]...),)
    end
end


function result(
    est::MagnetizationEstimator,
    ctx::MCContext,
    T::AbstractFloat,
    sign::AbstractFloat,
)
    _, signsymbols = magnetization_estimator_obs_symbols(est)
    norm = 1 / normalization_site_count(est.model)
    est.mag *= norm
    est.absmag *= norm
    est.mag2 *= norm^2
    est.mag4 *= norm^4

    measure!(ctx, signsymbols.mag, sign * est.mag / est.n)
    measure!(ctx, signsymbols.absmag, sign * est.absmag / est.n)
    measure!(ctx, signsymbols.mag2, sign * est.mag2 / est.n)
    measure!(ctx, signsymbols.mag4, sign * est.mag4 / est.n)

    chi =
        1 / T / (est.n + 1) / est.n *
        (est.mag^2 + est.mag2) *
        normalization_site_count(est.model)
    measure!(ctx, signsymbols.magchi, sign * chi)

    return nothing
end

get_prefix(
    ::Type{MagnetizationEstimator{OrderingVector,StaggerUC,Model,Prefix}},
) where {OrderingVector,StaggerUC,Model,Prefix} = Prefix

function register_evaluables(
    est::Type{<:MagnetizationEstimator},
    eval::Carlo.Evaluator,
    ::AbstractDict,
)
    symbols, signsymbols = magnetization_estimator_obs_symbols(get_prefix(est))

    for obs in (:mag, :absmag, :mag2, :mag4, :magchi)
        evaluate!(unsign, eval, symbols[obs], (signsymbols[obs], :Sign))
    end
    evaluate!(
        eval,
        symbols.binderratio,
        (signsymbols.mag4, signsymbols.mag2, :Sign),
    ) do smag4, smag2, sign
        if smag2 == smag4 == 0
            return zero(smag2)
        end
        return smag2^2 / smag4 / sign
    end

    return nothing
end
