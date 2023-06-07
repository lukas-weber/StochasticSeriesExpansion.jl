using StaticArrays

struct UCBond{D}
    iuc::Int
    jd::NTuple{D,Int} # position of second unit cell
    juc::Int
end

struct UCSite{D}
    pos::SVector{D,Float64}
    sublattice_sign::Int
    coordination::Int # filled automatically
end

UCSite(pos::Tuple) = UCSite(SVector(pos), 0, 0)

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

struct LatticeBond
    type::Int
    i::Int
    j::Int
end

struct LatticeSite{D}
    iuc::Int
    ix::NTuple{D,Int}
end

convert(::Type{Bond{2}}, lb::LatticeBond) = Bond(lb.type, (lb.i, lb.j))

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
