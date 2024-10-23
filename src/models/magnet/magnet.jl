struct MagnetBondParams{F<:Real}
    J::F
    d::F

    Dx::Tuple{F,F}
    Dz::Tuple{F,F}
    hz::Tuple{F,F}
end

function MagnetBondParams(
    J::Real,
    d::Real,
    Dx::Tuple{Real,Real},
    Dz::Tuple{Real,Real},
    hz::Tuple{Real,Real},
)
    T = promote_type(typeof.((J, d, Dx..., Dz..., hz...))...)
    return MagnetBondParams{T}(J, d, Dx, Dz, hz)
end

struct MagnetSiteParams
    spin_states::Int8
end


"""Describes the parameters of each bond

``H = \\sum_{ij} J_ij S_i\\cdotS_j + d_ij S_i^z S_j^z``
"""
struct MagnetModel <: AbstractModel
    lattice::Lattice

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

function MagnetModel(params::AbstractDict{Symbol,<:Any})
    Lx = params[:Lx]
    Ly = get(params, :Ly, Lx)
    lat = Lattice(params[:unitcell], (Lx, Ly))

    @assert Lx * Ly > 0
    @assert length(lat.bonds) > 0

    parameter_map = ParameterMap(get(params, :parameter_map, nothing))

    pm(path...) = map(parameter_map, path...)

    function split_site(param::Symbol, bond::LatticeBond, default)
        (iuc, _) = split_idx(lat, bond.i)
        (juc, _) = split_idx(lat, bond.j)

        first =
            get(params, pm(:sites, iuc, param), default) / lat.uc.sites[iuc].coordination
        second =
            get(params, pm(:sites, juc, param), default) / lat.uc.sites[juc].coordination

        return first, second
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
    return MagnetModel(lat, full_bond_params, full_site_params, opstring_ests)
end

function magnetization_state(mag::MagnetModel, site_idx::Integer, state_idx::Integer)
    return @fastmath (mag.site_params[site_idx].spin_states - 1) * 0.5 - state_idx + 1
end

normalization_site_count(mag::MagnetModel) = site_count(mag.lattice)

function generate_vertex_data(mag::MagnetModel, uc_bond, bond::MagnetBondParams)
    dimi = mag.site_params[uc_bond.iuc].spin_states
    dimj = mag.site_params[uc_bond.juc].spin_states

    splusi, szi = spin_operators(Float64, dimi)
    splusj, szj = spin_operators(Float64, dimj)
    idi = Matrix{Float64}(I, dimi, dimi)
    idj = Matrix{Float64}(I, dimj, dimj)

    H =
        bond.J * (
            0.5 * (kron(splusi', splusj) + kron(splusi, splusj')) +
            (1 + bond.d) * kron(szi, szj)
        ) +
        bond.hz[1] * kron(idi, szj) +
        bond.hz[2] * kron(szi, idj) +
        bond.Dx[1] / 4 * kron((splusi + splusi')^2, idj) +
        bond.Dx[2] / 4 * kron(idi, (splusj + splusj')^2) +
        bond.Dz[1] * kron(szi^2, idj) +
        bond.Dz[2] * kron(idi, szj^2)

    energy_offset_factor = 0.25
    # use the deterministic solution for S == 1//2
    if dimi == dimj == 2 && bond.hz == 0 && bond.Dz == 0 && bond.Dx == 0
        energy_offset_factor = 0.0
    end

    return VertexData((dimi, dimj), H; energy_offset_factor = energy_offset_factor)
end

leg_count(::Type{MagnetModel}) = 4

function generate_sse_data(mag::MagnetModel)
    vertex_data = [
        generate_vertex_data(mag, uc_bond, bond) for
        (uc_bond, bond) in zip(mag.lattice.uc.bonds, mag.bond_params)
    ]

    bonds = [SSEBond(bond.type, (bond.i, bond.j)) for bond in mag.lattice.bonds]

    return SSEData(vertex_data, bonds)
end

function gen_opstring_estimators(lattice::Lattice, params::AbstractDict)
    ests = DataType[]
    for est in get(params, :measure, [])
        if est == :magnetization
            q = ntuple(i -> false, dimension(lattice))
            push!(ests, MagnetizationEstimator{q,false,MagnetModel,Symbol()})
        elseif est == :staggered_magnetization
            neel = neel_vector(lattice.uc)
            if neel === nothing
                error(
                    "selected :staggered_magnetization measurement, but lattice does not have a Neel vector.",
                )
            end
            q, stagger_uc = neel
            push!(ests, MagnetizationEstimator{q,stagger_uc,MagnetModel,:Stag})
        elseif est <: AbstractOpstringEstimator
            push!(ests, est)
        else
            error("Unrecognized measure option '$est'")
        end
    end

    return ests
end

get_opstring_estimators(mag::MagnetModel) = mag.opstring_estimators
