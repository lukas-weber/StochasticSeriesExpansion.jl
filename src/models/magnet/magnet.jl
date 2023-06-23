module Models

using LinearAlgebra

import ..StochasticSeriesExpansion as S

struct MagnetBondParams{F<:AbstractFloat}
    J::F
    d::F

    Dx::F
    Dz::F
    hz::F
end

struct MagnetSiteParams
    spin_states::Int8
end


"""Describes the parameters of each bond

``H = \\sum_{ij} J_ij S_i\\cdotS_j + d_ij S_i^z S_j^z``
"""
struct Magnet <: S.AbstractModel
    lattice::S.Lattice

    bond_params::Vector{MagnetBondParams}
    site_params::Vector{MagnetSiteParams}

    opstring_estimators::Vector{DataType}
end

struct ParameterMap
    map::Union{Nothing,Dict}
end


function map(parameter_map::ParameterMap, path...)::Symbol
    if parameter_map.map === nothing
        return path[end]
    end

    res = parameter_map.map
    for p in path
        try
            res = res[p]
        catch e
            if e isa KeyError
                return path[end]
            end
        end
    end

    return res
end

function Magnet(params::AbstractDict{Symbol,<:Any})
    Lx = params[:Lx]
    Ly = get(params, :Ly, Lx)
    lat = S.Lattice(params[:unitcell], (Lx, Ly))

    @assert Lx * Ly > 0
    @assert length(lat.bonds) > 0

    parameter_map = ParameterMap(get(params, :parameter_map, nothing))

    pm(path...) = map(parameter_map, path...)

    function split_site(param::Symbol, bond::S.LatticeBond, default)
        (iuc, _) = S.split_idx(lat, bond.i)
        (juc, _) = S.split_idx(lat, bond.j)

        first =
            get(params, pm(:sites, iuc, param), default) / lat.uc.sites[iuc].coordination
        second =
            get(params, pm(:sites, juc, param), default) / lat.uc.sites[juc].coordination

        return (first + second) / 2
    end

    full_bond_params = [
        MagnetBondParams(
            params[pm(:bonds, bond.type, :J)],
            get(params, pm(:bonds, bond.type, :d), 0.0),
            split_site(:Dx, bond, 0.0),
            split_site(:Dz, bond, 0.0),
            split_site(:hz, bond, 0.0),
        ) for bond in lat.bonds
    ]

    full_site_params = repeat(
        [
            MagnetSiteParams(get(params, pm(:sites, i, :S), 1 // 2) * 2 + 1) for
            i in eachindex(lat.uc.sites)
        ],
        prod(lat.Ls),
    )

    opstring_ests = gen_opstring_estimators(lat, params)
    return Magnet(lat, full_bond_params, full_site_params, opstring_ests)
end

function S.magnetization_state(mag::Magnet, site_idx::Integer, state_idx::Integer)
    return @fastmath (mag.site_params[site_idx].spin_states - 1) * 0.5 - state_idx + 1
end

S.normalization_site_count(mag::Magnet) = S.site_count(mag.lattice)

function generate_vertex_data(mag::Magnet, uc_bond, bond::MagnetBondParams)
    dimi = mag.site_params[uc_bond.iuc].spin_states
    dimj = mag.site_params[uc_bond.juc].spin_states

    splusi, szi = S.spin_operators(Float64, dimi)
    splusj, szj = S.spin_operators(Float64, dimj)
    idi = diagm(ones(Float64, dimi))
    idj = diagm(ones(Float64, dimj))

    H =
        bond.J * (
            0.5 * (kron(splusi', splusj) + kron(splusi, splusj')) +
            (1 + bond.d) * kron(szi, szj)
        ) +
        bond.hz * (kron(idi, szj) + kron(szi, idj)) +
        bond.Dx / 4 * (kron((splusi + splusi')^2, idj) + kron(idi, (splusj + splusj')^2)) +
        bond.Dz * (kron(szi^2, idj) + kron(idi, szj^2))

    energy_offset_factor = 0.25
    # use the deterministic solution for S == 1//2
    if dimi == dimj == 2 && bond.hz == 0
        energy_offset_factor = 0.0
    end

    return S.VertexData((dimi, dimj), H; energy_offset_factor = energy_offset_factor)
end

S.leg_count(::Type{Magnet}) = 4

function S.generate_sse_data(mag::Magnet)
    vertex_data = [
        generate_vertex_data(mag, uc_bond, bond) for
        (uc_bond, bond) in zip(mag.lattice.uc.bonds, mag.bond_params)
    ]

    bonds = [S.SSEBond(bond.type, (bond.i, bond.j)) for bond in mag.lattice.bonds]

    return S.SSEData(vertex_data, bonds)
end

function gen_opstring_estimators(lattice::S.Lattice, params::AbstractDict)
    ests = DataType[]
    for est in get(params, :measure, [])
        if est == :magnetization
            q = ntuple(i -> false, S.dimension(lattice))
            push!(ests, S.MagnetizationEstimator{q,false,Magnet,Symbol()})
        elseif est == :staggered_magnetization
            neel = S.neel_vector(lattice.uc)
            if neel === nothing
                error(
                    "selected :staggered_magnetization measurement, but lattice does not have a Neel vector.",
                )
            end
            q, stagger_uc = neel
            push!(ests, S.MagnetizationEstimator{q,stagger_uc,Magnet,:Stag})
        elseif est <: S.MagnetizationEstimator
            push!(ests, est)
        else
            error("Unrecognized measure option '$est'")
        end
    end

    return ests
end

S.get_opstring_estimators(mag::Magnet) = mag.opstring_estimators

end
