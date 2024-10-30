struct ClusterBasis{N,T}
    quantum_numbers::Vector{NTuple{N,T}}
    transformation::Matrix{T}
end

struct ClusterModel{InnerModel,N,T} <: AbstractModel
    inner_model::InnerModel
    basis::NTuple{N,ClusterBasis}
    cluster_ids::Vector{T}
end

function ClusterModel(params::AbstractDict{Symbol,<:Any})
    inner_model = params[:inner_model](params)
    bases = params[:cluster_bases]
    parameter_map = ParameterMap(get(params, :parameter_map, nothing))

    cluster_ids = [
        get(params, get_parameter(parameter_map, :cluster_id, siteidx), 1) for
        siteidx in eachindex(inner_model.lattice.uc.sites)
    ]

    nclusters = length(unique(cluster_ids))
    if nclusters != length(bases)
        error(
            "Number of cluster bases ($(length(bases))) does not match number of distinct cluster ids ($nclusters).",
        )
    end

    return ClusterModel(inner_model, bases, cluster_ids)
end

leg_count(::Type{<:ClusterModel{InnerModel}}) where {InnerModel} = leg_count(InnerModel)

normalization_site_count(model::ClusterModel) = normalization_site_count(model.inner_model)

"""
    lift_twobody_operator(op::AbstractMatrix, site_dims::Tuple, sites::Tuple)

Lift an operator `op` that has already been lifted to a two-body Hilbert space (e.g. via `kron`) to a higher-order Hilbert space.

`site_dims` are the local Hilbert space dimensions and `sites` is a doublet specifying to where in the total Hilbert space  the two sites of `op` should be mapped.
This function follows the column-major convention of `kron`, where the leftmost dimension in `site_dims` is actually the slowest moving index. So

```jldoctest
julia> A = rand(4,4)
julia> B = rand(2,2)
julia> StochasticSeriesExpansion.lift_twobody_operator(kron(A,B), (3,4,2), (2,3)) == kron(I(3), A,B,)
true
"""
function lift_twobody_operator(op::AbstractMatrix, site_dims::Tuple, sites::Tuple)
    # reverse everything to fit column-major kron convention
    site_dims = reverse(site_dims)
    sites = reverse(length(site_dims) .+ 1 .- sites)

    op4 = reshape(
        op,
        site_dims[sites[1]],
        site_dims[sites[2]],
        site_dims[sites[1]],
        site_dims[sites[2]],
    )

    res = zeros(eltype(op4), site_dims..., site_dims...)

    for I in CartesianIndices(Base.OneTo.(site_dims))
        Ituple = Tuple(I)
        leftidxs = getindex.(Ref(Ituple), sites)

        for J in CartesianIndices(axes(op4)[3:4])
            rightidxs = Tuple(J)

            Ktuple = Tuple(I)
            for (site, rightidx) in zip(sites, rightidxs)
                Ktuple = Base.setindex(Ktuple, rightidx, site)
            end
            K = CartesianIndex(Ktuple)

            res[CartesianIndex(Ituple), K] = op4[leftidxs..., rightidxs...]
        end
    end

    return reshape(res, prod(site_dims), prod(site_dims))
end

function build_cluster_hamiltonians(cluster_ids, uc_bonds, bonds, site_params)
    clusters = [
        findall(i -> cluster_ids[i] == id, eachindex(cluster_ids)) for
        id in unique(cluster_ids)
    ]
    cluster_dims = map(
        siteidxs ->
            Int.(tuple((site_params[siteidx].spin_states for siteidx in siteidxs)...)),
        clusters,
    )

    cluster_ordering = zeros(Int, length(cluster_ids))
    for cluster in clusters
        for (i, siteidx) in enumerate(cluster)
            cluster_ordering[siteidx] = i
        end
    end

    intracluster_hamiltonians = zeros.(prod.(cluster_dims), prod.(cluster_dims))
    intercluster_hamiltonians = Dict{
        @NamedTuple{iuc::Int, juc::Int, jd::Tuple{Int,Int}},
        eltype(intracluster_hamiltonians),
    }()

    for (uc_bond, bond) in zip(uc_bonds, bonds)
        (dimi, dimj), H, energy_offset_factor = generate_bond_hamiltonian(
            uc_bond,
            bond,
            (site_params[uc_bond.iuc], site_params[uc_bond.juc]),
        )

        if cluster_ids[uc_bond.iuc] == cluster_ids[uc_bond.juc] && all(iszero, uc_bond.jd)
            if uc_bond.iuc == uc_bond.juc
                error(
                    "found bond in Magnet connecting site to itself... not supported by ClusterModel",
                )
            end
            cluster_id = cluster_ids[uc_bond.iuc]
            intracluster_hamiltonians[cluster_id] += lift_twobody_operator(
                H,
                cluster_dims[cluster_id],
                (cluster_ordering[uc_bond.iuc], cluster_ordering[uc_bond.juc]),
            )
        else
            iuc_cluster = cluster_ids[uc_bond.iuc]
            juc_cluster = cluster_ids[uc_bond.juc]
            cluster_bond = (iuc = iuc_cluster, juc = juc_cluster, jd = uc_bond.jd)

            cluster_bond_dims = (cluster_dims[iuc_cluster]..., cluster_dims[juc_cluster]...)
            dim = prod(cluster_bond_dims)
            Hij = get!(
                () -> zeros(eltype(H), dim, dim),
                intercluster_hamiltonians,
                cluster_bond,
            )

            Hij .+= lift_twobody_operator(
                H,
                cluster_bond_dims,
                (
                    cluster_ordering[uc_bond.iuc],
                    length(cluster_dims[iuc_cluster]) + cluster_ordering[uc_bond.juc],
                ),
            )
        end
    end

    return intracluster_hamiltonians, intercluster_hamiltonians
end

function absorb_intracluster_hamiltonians(
    intracluster_hamiltonians,
    intercluster_hamiltonians,
)
    cluster_coordinations = [
        sum((bond.iuc == i) + (bond.juc == i) for bond in keys(intercluster_hamiltonians)) for i in eachindex(intracluster_hamiltonians)
    ]

    return Dict(
        bond =>
            intercluster_hamiltonian +
            kron(
                intracluster_hamiltonians[bond.iuc],
                I(size(intracluster_hamiltonians[bond.juc], 1)) /
                cluster_coordinations[bond.iuc],
            ) +
            kron(
                I(size(intracluster_hamiltonians[bond.iuc], 1)) /
                cluster_coordinations[bond.juc],
                intracluster_hamiltonians[bond.juc],
            ) for (bond, intercluster_hamiltonian) in intercluster_hamiltonians
    )
end

function get_opstring_estimators(model::ClusterModel)
    return []
end

generate_sse_data(cluster::ClusterModel) =
    generate_cluster_sse_data(cluster, cluster.inner_model)

function generate_cluster_sse_data(cluster::ClusterModel, mag::MagnetModel)
    intracluster_hamiltonians, intercluster_hamiltonians = build_cluster_hamiltonians(
        cluster.cluster_ids,
        mag.lattice.uc.bonds,
        mag.bond_params,
        mag.site_params,
    )

    num_clusters = length(intracluster_hamiltonians)
    bonds = unique(
        skipmissing(
            map(mag.lattice.bonds) do bond
                uc_bond = mag.lattice.uc.bonds[bond.type]

                iuc = cluster.cluster_ids[uc_bond.iuc]
                juc = cluster.cluster_ids[uc_bond.juc]

                i =
                    num_clusters * (fld1(bond.i, length(mag.lattice.uc.sites)) - 1) + iuc
                j =
                    num_clusters * (fld1(bond.j, length(mag.lattice.uc.sites)) - 1) + juc

                if iuc == juc && all(iszero, uc_bond.jd)
                    return missing
                end

                bond_type = findfirst(
                    b -> b.iuc == iuc && b.juc == juc && b.jd == uc_bond.jd,
                    collect(keys(intercluster_hamiltonians)),
                )

                @assert !isnothing(bond_type)

                return SSEBond(bond_type, (i, j))
            end,
        ),
    )

    intercluster_hamiltonians = absorb_intracluster_hamiltonians(
        intracluster_hamiltonians,
        intercluster_hamiltonians,
    )

    energy_offset_factor = 0.25
    vertex_data = [
        let U = kron(
                cluster.basis[bond.iuc].transformation,
                cluster.basis[bond.juc].transformation,
            )
            VertexData(
                (
                    size(intracluster_hamiltonians[bond.iuc], 1),
                    size(intracluster_hamiltonians[bond.juc], 1),
                ),
                U' * H * U;
                energy_offset_factor,
            )
        end for (bond, H) in intercluster_hamiltonians
    ]

    return SSEData(vertex_data, bonds)
end

module ClusterBases
import ..ClusterBasis

const dimer = ClusterBasis(
    [(0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (1.0, -1.0)],
    [
        0 1 0 0
        1/√2 0 1/√2 0
        -1/√2 0 1/√2 0
        0 0 0 1
    ],
)
end
