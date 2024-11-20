# [Parameters](@id parameters)

The parameters that are given in a StochasticSeriesExpansion.jl job script can be divided into three levels. The first level are parameters read by Carlo.jl, such as `sweeps`, `thermalization` and `binsize`. These are documented in the [Carlo.jl documentation](https://lukas-weber.github.io/Carlo.jl/stable/).
The second level are generic stochastic series expansion (SSE) parameters, which are documented on this page. The third level are model-specific parameters that are documented together with the respective model.

The two required generic SSE level parameters are

- `T`: the temperature of the simulation.
- `model`: the type of the [`AbstractModel`](@ref) implementation representing the model to be simulated.

Apart from this, there are several optional parameters that control some internal aspects of the SSE simulation. They have reasonable defaults and in the majority of cases do not need to be touched.

- `init_opstring_cutoff`: the initial length of the padded operator string. Will be grown dynamically at runtime. If you see it growing a lot after thermalization, you can give it a higher intial value.
- `diagonal_warmup_sweeps`: sometimes loop updates can perform very badly on dilute operator strings that appear in the first steps of the simulation. For this reason, a number of diagonal sweeps are performed at initialization before the first loop update. This parameter sets that number.
- `init_num_worms`: the number of worms or loops that should be formed per sweep is adjusted automatically by a simple controller loop. This is the starting value. Should not be less than one.
- `num_worms_attenuation_factor`: factor that sets the sensitivity of the worm number control loop. Higher values equilibrate quicker, but are more susceptible to noise.
- `target_worm_length_fraction`: the goal of the controller is to perform enough worm updates so that the combined lengths of the loops they form are `target_worm_length_fraction * num_operators`.
