
hamiltonian(cluster::S.ClusterModel) = hamiltonian(cluster.inner_model)

function calc_observables!(
    obs::AbstractDict{Symbol,<:Any},
    magnet::S.ClusterModel,
    ens::Ensemble,
)
    # for est in S.get_opstring_estimators(magnet)
    #     calc_magnetization!(obs, magnet, ens, est)
    # end

    return nothing
end
