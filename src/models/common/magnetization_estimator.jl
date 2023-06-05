"""Generic OpstringEstimator that can be used for all models that have (i) a field `lattice::Lattice` and (ii) implement [`magnetization_state(::Model, state_idx)`](@ref).

`OrderingVector` is the Fourier component to compute in units of π. For example ``(1,1)`` corresponds to ``(π,π)``. If `StaggerUC` is true, the sublattice sign of the unitcell is additionally taken into account."""
mutable struct MagnetizationEstimator{
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

# Model extension interface:
@stub magnetization_state(model::AbstractModel, site_idx::Integer, state_idx::Integer)

"""In models where additional degrees of freedom exist, this function maps sse site indices to physical lattice site indices"""
magnetization_lattice_site_idx(::AbstractModel, sse_site_idx::Integer) = sse_site_idx

function MagnetizationEstimator{OrderingVector,StaggerUC,Model,Prefix}(
    model::Model,
) where {OrderingVector,StaggerUC,Model,Prefix}
    return MagnetizationEstimator{OrderingVector,StaggerUC,Model,Prefix}(
        model,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    )
end

function stag_sign(
    est::MagnetizationEstimator{OrderingVector,StaggerUC},
    site_idx::Integer,
) where {OrderingVector,StaggerUC}
    iuc, ix = split_idx(est.model.lattice, site_idx)

    sign = StaggerUC ? est.model.lattice.unitcell.sites[iuc].sublattice_sign : 1
    sign *= (-1)^sum(OrderingVector .* ix)

    return sign
end

function init!(est::MagnetizationEstimator, state::AbstractVector)
    est.tmpmag = sum(
        site ->
            stag_sign(est, site) * magnetization_state(est.model, site, state[site]),
        1:site_count(est.model.lattice),
    )

    est.mag = est.tmpmag
    est.absmag = abs(est.tmpmag)
    est.mag2 = est.tmpmag^2
    est.mag4 = est.tmpmag^4

    return nothing
end

function measure(
    est::MagnetizationEstimator,
    op::OperCode,
    ::AbstractVector{<:UInt8},
    sse_data::SSEData,
)
    if !isdiagonal(op)
        bond = sse_data.bonds[get_bond(op)]
        vd = get_vertex_data(sse_data, get_bond(op))
        leg_state = get_leg_state(vd, get_vertex(op))

        nsites = length(vd.dims)
        for l = 1:nsites
            site = magnetization_lattice_site_idx(est.model, bond.sites[l])
            if site !== nothing
                est.tmpmag +=
                    stag_sign(est, site) * (
                        magnetization_state(est.model, site, leg_state[nsites+l]) -
                        magnetization_state(est.model, site, leg_state[l])
                    )
            end
        end
    end

    est.mag += est.tmpmag
    est.absmag += abs(est.tmpmag)
    est.mag2 += abs(est.tmpmag^2)
    est.mag4 += abs(est.tmpmag^4)

    est.n += 1

    return nothing
end

function get_prefix(::Type{MagnetizationEstimator{G,S,M,Prefix}}) where {G,S,M,Prefix}
    return Prefix
end

function result(
    est::MagnetizationEstimator,
    ctx::MCContext,
    T::AbstractFloat,
    sign::AbstractFloat,
)
    prefix = "Sign$(get_prefix(typeof(est)))"

    norm = 1 / normalization_site_count(est.model)
    est.mag *= norm
    est.absmag *= norm
    est.mag2 *= norm^2
    est.mag4 *= norm^4

    measure!(ctx, Symbol(prefix * "Mag"), sign * est.mag / est.n)
    measure!(ctx, Symbol(prefix * "AbsMag"), sign * est.absmag / est.n)
    measure!(ctx, Symbol(prefix * "Mag2"), sign * est.mag2 / est.n)
    measure!(ctx, Symbol(prefix * "Mag4"), sign * est.mag4 / est.n)

    chi =
        1 / T / (est.n + 1) / est.n *
        (est.mag^2 + est.mag2) *
        normalization_site_count(est.model)
    measure!(ctx, Symbol(prefix * "MagChi"), sign * chi)

    return nothing
end

function register_evaluables(
    est::Type{<:MagnetizationEstimator},
    eval::LoadLeveller.Evaluator,
)
    prefix = get_prefix(est)

    evaluate!(unsign, eval, Symbol("$(prefix)Mag"), [Symbol("Sign$(prefix)Mag"), :Sign])
    evaluate!(
        unsign,
        eval,
        Symbol("$(prefix)AbsMag"),
        [Symbol("Sign$(prefix)AbsMag"), :Sign],
    )
    evaluate!(unsign, eval, Symbol("$(prefix)Mag2"), [Symbol("Sign$(prefix)Mag2"), :Sign])
    evaluate!(unsign, eval, Symbol("$(prefix)Mag4"), [Symbol("Sign$(prefix)Mag4"), :Sign])
    evaluate!(
        unsign,
        eval,
        Symbol("$(prefix)MagChi"),
        [Symbol("Sign$(prefix)MagChi"), :Sign],
    )

    evaluate!(
        eval,
        Symbol("$(prefix)BinderRatio"),
        [Symbol("Sign$(prefix)Mag4"), Symbol("Sign$(prefix)Mag2"), :Sign],
    ) do smag4, smag2, sign
        if smag2 == smag4 == 0
            return zero(smag2)
        end
        return smag2^2 / smag4 / sign
    end

    return nothing
end
