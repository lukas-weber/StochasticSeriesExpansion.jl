using Test
using LinearAlgebra
import StochasticSeriesExpansion as S

include("models/test_cluster.jl")
include("models/test_magnet.jl")

include("test_jobs.jl")
include("test_vertex_list.jl")
include("test_worms.jl")
include("test_vertex_data.jl")
include("test_sse.jl")
include("test_sse_data.jl")
include("test_opercode.jl")
include("test_util.jl")
include("test_common_operators.jl")
include("test_lattice.jl")
include("test_ed_compare.jl")
