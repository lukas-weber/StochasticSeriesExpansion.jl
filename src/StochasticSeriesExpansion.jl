module StochasticSeriesExpansion
export AbstractModel,
    SSEData, AbstractOpstringEstimator, Bond, VertexData, get_vertex_data, get_leg_state

export MagnetizationEstimator,
    all_magests, magnetization_state, magnetization_lattice_site_idx, magest_standard_prefix

using Carlo

include("util.jl")
include("worms.jl")
include("opercode.jl")
include("vertex_data.jl")
include("sse_data.jl")
include("vertex_list.jl")
include("abstract_model.jl")
include("opstring_estimator.jl")
include("sse.jl")

include("models/common/operators.jl")
include("models/common/lattice.jl")
include("models/common/magnetization_estimator.jl")
include("models/magnet/magnet.jl")

end # module StochasticSeriesExpansion
