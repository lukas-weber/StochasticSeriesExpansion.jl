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
    magnet::S.MagnetModel,
    ens::Ensemble,
    ::Type{S.MagnetizationEstimator{OrderingVector,StaggerUC,Model,Prefix}},
) where {OrderingVector,StaggerUC,Model,Prefix}
    calc_magnetization!(
        obs,
        magnet,
        ens;
        ordering_vector = OrderingVector,
        stagger_uc = StaggerUC,
        prefix = Prefix,
    )
end

function calc_magnetization!(
    obs::AbstractDict{Symbol,<:Any},
    magnet::S.MagnetModel,
    ens::Ensemble;
    ordering_vector::Tuple,
    stagger_uc::Bool,
    prefix::Symbol,
)

    M =
        make_operator(magnet) do M, lifter
            for i in eachindex(magnet.site_params)
                M +=
                    S.staggered_sign(magnet.lattice, ordering_vector, stagger_uc, i) *
                    spin(lifter, i, 3)
            end
            return M
        end ./ S.normalization_site_count(magnet)

    symbols, _ = S.magnetization_estimator_obs_symbols(prefix)

    obs[symbols.mag] = mean(ens, M)
    m2 = mean(ens, M^2)
    m4 = mean(ens, M^4)
    obs[symbols.mag2] = m2
    obs[symbols.mag4] = m4
    obs[symbols.absmag] = mean(ens, abs.(M))

    obs[symbols.binderratio] = m2 .^ 2 ./ m4
    obs[symbols.magchi] =
        integrated_correlator(ens, M, M) * S.normalization_site_count(magnet)

    return nothing
end

function calc_observables!(
    obs::AbstractDict{Symbol,<:Any},
    magnet::S.MagnetModel,
    ens::Ensemble,
)
    for est in S.get_opstring_estimators(magnet)
        calc_magnetization!(obs, magnet, ens, est)
    end

    return nothing
end
