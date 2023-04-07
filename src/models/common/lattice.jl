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

UCSite(pos::Tuple, sublattice_sign::Integer = 1) = UCSite(SVector(pos), sublattice_sign, 1)

struct UnitCell{D,F}
    lattice_vectors::SMatrix{D,D,F} # as columns

    sites::Vector{UCSite{D,F}}
    bonds::Vector{UCBond{D}}
end

UnitCell(
    lattice_vectors::AbstractMatrix,
    sites::Vector{<:UCSite{D,F}},
    bonds::Vector{<:UCBond{D}},
) where {D,F} = UnitCell(SMatrix{D,D}(lattice_vectors), sites, bonds)

struct Lattice{D,F}
    uc::UnitCell{D,F}
    Ls::NTuple{D,Int}

    bonds::Vector{NamedTuple{(:type, :i, :j),Tuple{Int,Int,Int}}}
end

function split_idx(l::Lattice{D,F}, site_idx::Integer) where {D,F}
    iuc = site_idx % length(l.uc.sites)
    site_idx ÷= length(l.uc.sites)

    r = MVector{D,Int}()
    for (i, L) in enumerate(l.Ls)
        r[i] = site_idx % L
        site_idx ÷= L
    end

    return iuc, SVector(r)
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
    [UCBond(0, (0, 1), 0), UCBond(0, (1, 0), 0)],
)
const columnar_dimer = UnitCell(
    [1.0 0.0; 0.0 2.0],
    [UCSite((0.0, 0.0)), UCSite((0.0, 0.5))],
    [
        UCBond(0, (0, 0), 1),
        UCBond(0, (1, 0), 0),
        UCBond(1, (1, 0), 1),
        UCBond(1, (0, 1), 0),
    ],
)

const honeycomb = UnitCell(
    [[sqrt(3) / 2, -0.5] [sqrt(3 / 2), 0.5]],
    [UCSite((0.0, 0.0)), UCSite((1 / 3, 1 / 3))],
    [UCBond(0, (0, 0), 1), UCBond(1, (0, 1), 0), UCBond(1, (1, 0), 0)],
)

end
