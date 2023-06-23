using StaticArrays

"""
    UCBond{D}

Defines a unit cell bond between site `i` and `j`. Construct it using

    UCBond(iuc, jd, juc)

`iuc` is the intra-unit-cell index of `i`. `jd = (jdx, jdy, ...)` are the `D` (nonnegative) offsets of the unit cell containing `j` with respect to the one containing `i`. `juc` is the intra-unit-cell index of `j`. For examples, see [`UnitCells`](@ref).

These parameters are accessible as public fields as well.
"""
struct UCBond{D}
    iuc::Int
    jd::NTuple{D,Int}
    juc::Int
end

"""
    UCSite{D}

Defines a site within the unit cell. Construct it using

    UCSite(pos::Tuple)

`pos` contains the position of the site inside the unit cell in the lattice vector basis.

## Public fields:

* `pos`: see above
* `coordination`: coordination number of the site
"""
struct UCSite{D}
    pos::SVector{D,Float64}
    sublattice_sign::Int
    coordination::Int
end

UCSite(pos::Tuple) = UCSite(SVector(pos), 0, 0)

"""
    UnitCell{D}

The unit cell of a `D`-dimensional lattice.

    UnitCell(lattice_vectors::SMatrix{D,D, Float64}, sites::AbstractVector{UCSite{D}}, bonds::AbstractVector{UCBond{D}})

The columns of `lattice_vector` are the lattice vectors. `sites` and `bonds` define the `bonds` within the unit cell.
"""
struct UnitCell{D}
    lattice_vectors::SMatrix{D,D,Float64} # as columns

    sites::Vector{UCSite{D}}
    bonds::Vector{UCBond{D}}
end

function calculate_uc_signs(bonds::AbstractVector{<:UCBond}, num_sites::Integer)
    signs = zeros(Int, num_sites)

    signs[1] = 1
    tries = 0
    while any(signs .== 0) || tries > length(bonds)^2
        for b in bonds
            if any(b.jd .!= 0)
                continue
            end

            if signs[b.iuc] != 0 && signs[b.juc] == 0
                signs[b.juc] = -signs[b.iuc]
            elseif signs[b.juc] != 0 && signs[b.iuc] == 0
                signs[b.iuc] = -signs[b.juc]
            elseif signs[b.iuc] == signs[b.juc]
                signs .= 1 # lattice not bipartite
                break
            end
            tries += 1
        end
    end

    if any(signs .== 0)
        signs .= 1
    end

    return signs
end

function calculate_uc_coordinations(bonds::AbstractVector{<:UCBond}, num_sites::Integer)
    coordinations = zeros(Int, num_sites)

    for b in bonds
        coordinations[b.iuc] += 1
        coordinations[b.juc] += 1
    end

    return coordinations
end

function UnitCell(
    lattice_vectors::AbstractMatrix,
    sites::Vector{<:UCSite{D}},
    bonds::Vector{<:UCBond{D}},
) where {D}

    signs = calculate_uc_signs(bonds, length(sites))
    coordinations = calculate_uc_coordinations(bonds, length(sites))

    return UnitCell(
        SMatrix{D,D}(lattice_vectors),
        [
            UCSite(site.pos, sign, coordination) for
            (site, sign, coordination) in zip(sites, signs, coordinations)
        ],
        bonds,
    )
end

function neel_vector(uc::UnitCell{D})::Union{Nothing,Tuple{NTuple{D,Bool},Bool}} where {D}
    for stagger_uc in (false, true)
        for q in Iterators.product(((false, true) for _ = 1:D)...)
            if all(uc.bonds) do bond
                signi = uc.sites[bond.iuc].sublattice_sign^stagger_uc
                signj =
                    uc.sites[bond.juc].sublattice_sign^stagger_uc * (-1)^(sum(bond.jd .* q))

                return signi != signj
            end
                return q, stagger_uc
            end
        end
    end

    return nothing
end

"""
    LatticeBond

Information about a lattice bond.

## Public fields

* `type`: translation symmetry equivalence class of the bond (does not take other symmetries into account!)
* `i`, `j`: sites connected by the bond

!!! note
In contrast to [`SSEData`](@ref), `Lattice` only supports 2-site bonds at this point. 
"""
struct LatticeBond
    type::Int
    i::Int
    j::Int
end

"""
    LatticeSite{D}

Information about a lattice site, locating its position inside the lattice.

## Public fields

* `iuc`: index within the unit cell
* `ix`: `D`-tuple of indices of the unit cell along the lattice dimensions
"""
struct LatticeSite{D}
    iuc::Int
    ix::NTuple{D,Int}
end

"""
    Lattice{D}

This is a model-agnostic implementation of a `D`-dimensional lattice. It allows generating translation symmetric bond graph based on
minimal information derived from a unit cell object.

## Public fields

* `uc`: underlying `UnitCell` object
* `Ls`: tuple of lattice lengths
* `bonds`: array of [`LatticeBond`](@ref).
* `site`: array of [`LatticeSite`](@ref), useful for quickly converting linear site indices to spatial ones.

"""
struct Lattice{D}
    uc::UnitCell{D}
    Ls::NTuple{D,Int}

    bonds::Vector{LatticeBond}
    sites::Vector{LatticeSite{D}}
end

dimension(lat::Lattice{D}) where {D} = D
dimension(unitcell::UnitCell{D}) where {D} = D

function Lattice(uc::UnitCell{D}, Ls::NTuple{D,<:Integer}) where {D}
    dims = (length(uc.sites), Ls...)
    bonds = LatticeBond[]
    sites = LatticeSite[]
    for r in Iterators.product([1:L for L in Ls]...)
        for (bond_type, b) in enumerate(uc.bonds)
            @assert 1 <= b.iuc <= length(uc.sites)
            @assert 1 <= b.juc <= length(uc.sites)

            i = join_idx(dims, (b.iuc, r...))
            j = join_idx(dims, (b.juc, ((r .+ b.jd .- 1) .% Ls .+ 1)...))

            push!(bonds, LatticeBond(bond_type, i, j))
        end
        for (iuc, uc_site) in enumerate(uc.sites)
            push!(sites, LatticeSite(iuc, r))
        end
    end

    return Lattice{D}(uc, Ls, bonds, sites)
end

function split_idx(l::Lattice, site_idx::Integer)
    r = split_idx((length(l.uc.sites), l.Ls...), site_idx)
    return r[1], r[2:end]
end

site_count(l::Lattice)::Int = length(l.uc.sites) * prod(l.Ls)

function site_pos(l::Lattice, site_idx::Integer)
    (iuc, r) = split_idx(l, site_idx)
    return l.uc.lattice_vectors * (r .+ l.uc.sites[iuc].pos)
end

function staggered_sign(
    l::Lattice{D},
    ordering_vector::NTuple{D,Bool},
    stagger_uc::Bool,
    site_idx::Integer,
)::Int where {D}
    sign = stagger_uc ? l.uc.sites[l.sites[site_idx].iuc].sublattice_sign : 1
    sign *= 1 - 2 * isodd(sum(ordering_vector .* l.sites[site_idx].ix))

    return sign
end

"""
The `UnitCells` module provides some predefined unitcells that can be used to construct lattices. It is also a helpful resource to see how to define your own unit cells.
"""
module UnitCells
import ..UnitCell, ..UCSite, ..UCBond

const square = UnitCell(
    [1.0 0.0; 0.0 1.0],
    [UCSite((0.0, 0.0))],
    [UCBond(1, (0, 1), 1), UCBond(1, (1, 0), 1)],
)

const columnar_dimer = UnitCell(
    [1.0 0.0; 0.0 2.0],
    [UCSite((0.0, 0.0)), UCSite((0.0, 0.5))],
    [
        UCBond(1, (0, 0), 2),
        UCBond(1, (1, 0), 1),
        UCBond(2, (1, 0), 2),
        UCBond(2, (0, 1), 1),
    ],
)

const honeycomb = UnitCell(
    [[sqrt(3) / 2, -0.5] [sqrt(3 / 2), 0.5]],
    [UCSite((0.0, 0.0)), UCSite((1 / 3, 1 / 3))],
    [UCBond(1, (0, 0), 2), UCBond(2, (0, 1), 1), UCBond(2, (1, 0), 1)],
)

const triangle = UnitCell(
    [[1.0, 0.0] [-0.5, sqrt(3 / 2)]],
    [UCSite((0.0, 0.0))],
    [UCBond(1, (0, 1), 1), UCBond(1, (1, 0), 1), UCBond(1, (1, 1), 1)],
)
end
