using StaticArrays

struct UCBond{D}
    iuc::Int
    jd::NTuple{D,Int} # position of second unit cell
    juc::Int
end

struct UCSite{D,F<:Real}
    pos::SVector{D,F}
    sublattice_sign::Int
    coordination::Int # filled automatically
end

UCSite(pos::Tuple) = UCSite(SVector(pos), 0, 0)

struct UnitCell{D,F}
    lattice_vectors::SMatrix{D,D,F} # as columns

    sites::Vector{UCSite{D,F}}
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
    sites::Vector{<:UCSite{D,F}},
    bonds::Vector{<:UCBond{D}},
) where {D,F}

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

struct LatticeBond
    type::Int
    i::Int
    j::Int
end

struct Lattice{D,F}
    uc::UnitCell{D,F}
    Ls::NTuple{D,Int}

    bonds::Vector{LatticeBond}
end

function Lattice(uc::UnitCell{D,F}, Ls::NTuple{D,<:Integer}) where {D,F}
    dims = (length(uc.sites), Ls...)
    bonds = LatticeBond[]
    for r in Iterators.product([1:L for L in Ls]...)
        for (bond_type, b) in enumerate(uc.bonds)
            @assert 1 <= b.iuc <= length(uc.sites)
            @assert 1 <= b.juc <= length(uc.sites)

            i = join_idx(dims, (b.iuc, r...))
            j = join_idx(dims, (b.juc, ((r .+ b.jd .- 1) .% Ls .+ 1)...))

            push!(bonds, LatticeBond(bond_type, i, j))
        end
    end

    return Lattice{D,F}(uc, Ls, bonds)
end

function split_idx(l::Lattice, site_idx::Integer)
    r = split_idx((length(l.uc.sites), l.Ls...), site_idx)
    return r[1], r[2:end]
end

site_count(l::Lattice) = length(l.uc.sites) * prod(l.Ls)

function site_pos(l::Lattice, site_idx::Integer)
    (iuc, r) = split_idx(l, site_idx)
    return l.uc.lattice_vectors * (r .+ l.uc.sites[iuc].pos)
end

function sublattice_sign(l::Lattice, site_idx::Integer)
    return l.uc.sites[site_idx%length(l.uc.sites)].sublattice_sign
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

end
