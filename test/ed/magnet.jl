function make_operator(func, magnet::S.Models.Magnet)
    dims = [Int(site.spin_mag * 2 + 1) for site in magnet.site_params]

    lifter = Lifter(dims)
    op = spzeros(prod(dims), prod(dims))
    op = func(op, lifter)
    return op
end

function hamiltonian(magnet::S.Models.Magnet)
    return make_operator(magnet) do H, lifter
        for (bond, params) in zip(magnet.lattice.bonds, magnet.bond_params)
            H +=
                params.J * heisen_bond(lifter, bond.i, bond.j) +
                params.d * spin(lifter, bond.i, 3) * spin(lifter, bond.j, 3) +
                params.hz * (spin(lifter, bond.i, 3) + spin(lifter, bond.j, 3))
            +params.Dx * (spin(lifter, bond.i, 1)^2 + spin(lifter, bond.j, 1)^2) +
            params.Dz * (spin(lifter, bond.i, 1)^2 + spin(lifter, bond.j, 3)^2)
        end
        return H
    end
end


function calc_observables!(
    obs::AbstractDict{Symbol,<:Any},
    magnet::S.Models.Magnet,
    ens::Ensemble,
)
    M = make_operator(magnet) do M, lifter
        for i in eachindex(magnet.site_params)
            M += spin(lifter, i, 3)
        end
        return M
    end./S.normalization_site_count(magnet)

    obs[:Mag] = mean(ens, M)
    obs[:Mag2] = mean(ens, M^2)
    obs[:Mag4] = mean(ens, M^4)
    obs[:AbsMag] = mean(ens, abs.(M))

    
    obs[:BinderRatio] = ((m4, m2)-> m4==m2==0 ? zero(m2) : m4/m2^2).(obs[:Mag4], obs[:Mag2])
end
