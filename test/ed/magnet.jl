function make_operator(func, magnet::S.MagnetModel)
    dims = [site.spin_states for site in magnet.site_params]

    lifter = Lifter(dims)
    op = spzeros(prod(dims), prod(dims))
    op = func(op, lifter)
    return op
end

function hamiltonian(magnet::S.MagnetModel)
    return make_operator(magnet) do H, lifter
        for (bond, params) in zip(magnet.lattice.bonds, magnet.bond_params)
            H +=
                params.J * heisen_bond(lifter, bond.i, bond.j) +
                params.J * params.d * spin(lifter, bond.i, 3) * spin(lifter, bond.j, 3) +
                params.hz[1] * spin(lifter, bond.i, 3) +
                params.hz[2] * spin(lifter, bond.j, 3) +
                params.Dx[1] * spin(lifter, bond.i, 1)^2 +
                params.Dx[2] * spin(lifter, bond.j, 1)^2 +
                params.Dz[1] * spin(lifter, bond.i, 3)^2 +
                params.Dz[2] * spin(lifter, bond.j, 3)^2
        end
        return H
    end
end

function calc_magnetization!(
    obs::AbstractDict{Symbol,<:Any},
    model::S.AbstractModel,
    ens::Ensemble,
    ::Type{S.MagnetizationEstimator{OrderingVector,StaggerUC,Model,Prefix,Tag}},
) where {OrderingVector,StaggerUC,Model,Prefix,Tag}
    calc_magnetization!(
        obs,
        model,
        ens;
        ordering_vector = OrderingVector,
        stagger_uc = StaggerUC,
        prefix = Prefix,
        tag = Tag,
    )
end

function calc_magnetization!(
    obs::AbstractDict{Symbol,<:Any},
    model::S.AbstractModel,
    ens::Ensemble;
    ordering_vector::Tuple,
    stagger_uc::Bool,
    prefix::Symbol,
    tag,
)

    M =
        make_operator(model) do M, lifter
            for i in eachindex(model.site_params)
                M +=
                    S.staggered_sign(model.lattice, ordering_vector, stagger_uc, i) *
                    lift(
                        lifter,
                        i,
                        Diagonal([
                            S.magnetization_state(model, Val(tag), i, s) for
                            s = 1:lifter.site_dims[i]
                        ]),
                    )
            end
            return M
        end ./ S.normalization_site_count(model)

    symbols, _ = S.magnetization_estimator_obs_symbols(prefix)

    obs[symbols.mag] = mean(ens, M)
    m2 = mean(ens, M^2)
    m4 = mean(ens, M^4)
    obs[symbols.mag2] = m2
    obs[symbols.mag4] = m4
    obs[symbols.absmag] = mean(ens, abs.(M))

    obs[symbols.binderratio] = m2 .^ 2 ./ m4
    obs[symbols.magchi] =
        integrated_correlator(ens, M, M) * S.normalization_site_count(model)

    return nothing
end

function calc_observables!(
    obs::AbstractDict{Symbol,<:Any},
    model::S.AbstractModel,
    ens::Ensemble,
)
    for est in S.get_opstring_estimators(model)
        calc_magnetization!(obs, model, ens, est)
    end

    return nothing
end
