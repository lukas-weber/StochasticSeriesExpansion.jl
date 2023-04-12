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
    spin_mag::Rational
end


"""Describes the parameters of each bond

``H = \\sum_{ij} J_ij S_i\\cdotS_j + d_ij S_i^z S_j^z``
"""
struct Magnet{D,F} <: S.AbstractModel
    lat::S.Lattice{D,F}

    bond_params::Vector{MagnetBondParams{F}}
    site_params::Vector{MagnetSiteParams}
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
        if !haskey(res, p)
            return path[end]
        end
        res = res[p]
    end

    return res
end

function Magnet(params::AbstractDict{Symbol,<:Any})
    Lx = params[:Lx]
    Ly = get(params, :Ly, Lx)
    lat = S.Lattice(params[:unitcell], (Lx, Ly))

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
            MagnetSiteParams(get(params, pm(:sites, i, :S), 1 // 2)) for
            i in eachindex(lat.uc.sites)
        ],
        prod(lat.Ls),
    )

    return Magnet(lat, full_bond_params, full_site_params)
end

function magnetization_state(
    mag::Magnet{D,F},
    site_idx::Integer,
    state_idx::Integer,
)::F where {D,F}
    return -mag.site_params[site_idx].spin_mag + state_idx - 1
end

function generate_vertex_data(mag::Magnet{D,F}, uc_bond, bond::MagnetBondParams) where {D,F}
    dimi = Int(mag.site_params[uc_bond.iuc].spin_mag * 2 + 1)
    dimj = Int(mag.site_params[uc_bond.juc].spin_mag * 2 + 1)

    splusi, szi = S.spin_operators(F, dimi)
    splusj, szj = S.spin_operators(F, dimj)
    idi = diagm(ones(F, dimi))
    idj = diagm(ones(F, dimj))

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
        (uc_bond, bond) in zip(mag.lat.uc.bonds, mag.bond_params)
    ]

    sites = [S.Site(Int(s.spin_mag * 2 + 1)) for s in mag.site_params]
    bonds = [S.Bond(bond.type, (bond.i, bond.j)) for bond in mag.lat.bonds]

    return S.SSEData(vertex_data, sites, bonds)
end

end
