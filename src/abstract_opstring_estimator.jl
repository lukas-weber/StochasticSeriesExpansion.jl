using Carlo

"""
This interface allows defining observable estimators that act on each operator of the the SSE operator string. The measurement code will call the functions of this interface like the following. However, for efficiency reasons, multiple estimators are automatically fused into a single loop.

```julia 
est = init(YourEstimator, model, state)
for op in operators
    if isidentity(op)
        continue
    end

    if !isdiagonal(op)
        # update state
    end

    if n < num_operators
        measure(est, op, state, sse_data)
    end
    n += 1
end

result(est, mccontext, T, sign)
```

However, in practice, StochasticSeriesExpansion will interlace different estimators into the same loop for efficiency.

A reference implementation is in the model-generic [`MagnetizationEstimator`](@ref), which can be used for most magnetization-like observables.
"""
abstract type AbstractOpstringEstimator end

"""
    init(::Type{YourOpstringEstimator},
         model::AbstractModel,
         state::AbstractVector{<:StateIdx}) -> YourOpstringEstimator

Constructs an opstring estimator based on a `model` and an initial `state`. The `state` is a vector of integers labelling
the computational site basis states. The estimator needs to interpret them in terms of physical quantities.
"""
function init end

"""
    measure(est, op::OperCode, state::AbstractVector{<:StateIdx}, sse_data::SSEData)

Perform the in-string measurement of estimator `est` on each operator `op` in the operator string. The `state` at the current position and the `sse_data` object are passed for additional reference.
"""
function measure end

"""
    result(est, ctx::Carlo.MCContext, T::AbstractFloat, sign::AbstractFloat)

Finalize the measurement by saving the results to the Carlo `MCContext`, e.g. by calling
`measure!(ctx, :Magnetization, est.mag)`. For more information, consult the Carlo documentation.

For some observables, knowing the temperature `T` is necessary. In the case of a signful simulation, `sign != 1` should be taken into account.
"""
function result end

"""
    register_evaluables(::Type{<:YourOpstringEstimator}, eval::Carlo.Evaluator, params::AbstractDict)

Operator string estimators `est` can define their own evaluables using this function, which passes a
`Carlo.Evaluator` and the task parameters. The state of the estimator is unavailable here since this runs in the postprocessing step.
"""
function register_evaluables end
