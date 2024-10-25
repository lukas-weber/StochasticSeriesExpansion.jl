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


@doc raw"""
    MagnetModel <: AbstractModel

Describes an arbitrary-spin quantum magnet with the Hamiltonian
```math
H = \sum_{⟨i,j⟩} J_{ij} S_i \cdot S_j + d_{ij} S_i^z S_j^z + \sum_i h_i S^z_i + D^z_i (S^z_i)^2 + D^x_i (S^x_i)^2
```

## Task parameters

- `lattice`: sets the [`Lattice`](@ref)

- `S`: spin magnitude (default: `1//2`)
- `J`: exchange coupling ``J_{ij}``
- `d`: exchange anisotrozy ``d_{ij}`` (default: `0`)
- `h`: magnetic field in ``z``-direction ``h_i`` (default: `0`)
- `D_z`: single-ion anisotropy in ``z``-direction ``D^z_i`` (default: `0`)
- `D_x`: single-ion anisotropy in ``x``-direction ``D^x_i`` (default: `0`)

By default, these are the same for each site or bond. However, they (including `:S` )can be adjusted to be different across the unit cell using
- `parameter_map`: [parameter map](@ref parameter_maps) that assigns different parameters to different bonds/sites of the unit cell.

Further parameters:
- `measure`: control what observables should be measured. Vector containing a combination of
    * `:magnetization`: measures uniform magnetization
    * `:staggered_magnetization`: measures bipartite staggered (antiferromagnetic) magnetization
    * any `Type{<:AbstractOpstringEstimator}` for an implementation of the [operator string estimator interface](@ref abstract_opstring_estimator)
"""
struct MagnetModel <: AbstractModel
    lattice::Lattice

    bond_params::Vector{MagnetBondParams}
    site_params::Vector{MagnetSiteParams}

    opstring_estimators::Vector{DataType}
end

struct ParameterMap{Map}
    map::Map
end


function map(parameter_map::ParameterMap, name::Symbol, index::Int)::Symbol
    if parameter_map.map === nothing
        return name
    end

    res = parameter_map.map
    if !haskey(res, name) || !(index in eachindex(res[name]))
        return name
    end

    return res[name][index]
end

function MagnetModel(params::AbstractDict{Symbol,<:Any})
    lattice = Lattice(params[:lattice])

    @assert length(lattice.bonds) > 0

    parameter_map = ParameterMap(get(params, :parameter_map, nothing))

    pm(path...) = map(parameter_map, path...)

    function split_site(param::Symbol, bond::LatticeBond, default)
        (iuc, _) = split_idx(lattice, bond.i)
        (juc, _) = split_idx(lattice, bond.j)

        first = get(params, pm(param, iuc), default) / lattice.uc.sites[iuc].coordination
        second = get(params, pm(param, juc), default) / lattice.uc.sites[juc].coordination

        return first, second
    end

    full_bond_params = [
        MagnetBondParams(
            params[pm(:J, bond.type)],
            get(params, pm(:d, bond.type), 0.0),
            split_site(:Dx, bond, 0.0),
            split_site(:Dz, bond, 0.0),
            split_site(:hz, bond, 0.0),
        ) for bond in lattice.bonds
    ]

    full_site_params = repeat(
        [
            MagnetSiteParams(get(params, pm(:S, i), 1 // 2) * 2 + 1) for
            i in eachindex(lattice.uc.sites)
        ],
        prod(lattice.Ls),
    )

    opstring_ests = gen_opstring_estimators(lattice, params)
    return MagnetModel(lattice, full_bond_params, full_site_params, opstring_ests)
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
