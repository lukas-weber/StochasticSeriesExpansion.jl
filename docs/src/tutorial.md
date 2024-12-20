# [Tutorial](@id tutorial)

In this tutorial, we will calculate the magnetic susceptibility of the spin-1 honeycomb antiferromagnet BaNi₂V₂O₈, reproducing a calculation in [PRB **104**, 065502 (2021)](https://doi.org/10.1103/physrevb.104.064402), where an anisotropic Heisenberg model,

```math
    H = J \sum_{⟨i,j⟩} \mathbf{S}_i \cdot \mathbf{S}_j + D^z \sum_i (S^z_i)^2
```
 was fit to neutron-scattering data. Above, the first sum counts all nearest neighbor bonds (and we neglect next-nearest and next-next-nearest neighbor interactions). However, we do take into account the easy-plane anisotropy ``D^z`` to model the split between in-plane and out of plane magnetic susceptibility.

## Installation

As a first step, we need to install StochasticSeriesExpansion.jl and [Carlo.jl](https://lukas-weber.github.io/Carlo.jl/stable/), which can be done from the Julia REPL.

    julia>]
    pkg> add Carlo
    pkg> add StochasticSeriesExpansion

## Writing a job script

To specify the parameters we want to simulate, we want to write a job script. The following sets up the predefined [`MagnetModel`](@ref) for our Hamiltonian on a honeycomb lattice and defines tasks for a range of system sizes and temperatures to simulate.
```@eval
using Markdown
path = "../../examples/bani2v2o8.jl"
code = open(path, "r") do file
    return read(file, String)
end

Markdown.parse("```julia\n$code\n```")
```
The parameters of this job script control different components of StochasticSeriesExpansion.jl. First are the parameters for Carlo.jl, namely the number of `sweeps` (Monte Carlo timesteps), the number of `thermalization` steps and the internal `binsize` that is averaged before saving data to disk (use it to save disk space, but keep it small compared to `sweeps`).

Then, `T` and `model` set the temperature and the model (Hamiltonian) that should be simulated by StochasticSeriesExpansion.jl. See [Parameters](@ref parameters) for a complete list of parameters at this level. The remainder of the parameters are model parameters, interpreted in this case by the [`MagnetModel`](@ref).

`J` and `Dz` which correspond to the parameters in the Hamiltonian, there is the spin magnitude `S` and the vector `measure` which decides what kinds of observables should be measured. `:magnetization` is a shortcut for observables related to the magnetization in ``S^z`` direction, including the susceptibility. The meaning of all the model-specific parameters is listed in the documentation of [`MagnetModel`](@ref).

The `task` function takes a snapshot of the current properties of the [`TaskMaker`](https://lukas-weber.github.io/Carlo.jl/stable/jobtools.html#Carlo.JobTools.TaskMaker) and turns it into one parameter set to be simulated.

In the [`JobInfo`](https://lukas-weber.github.io/Carlo.jl/stable/jobtools.html#Carlo.JobTools.JobInfo) structure we specify

* the job directory (`splitext(@__FILE__)[1]` just evaluates to the location where the job script is saved)
* that we want to run a simulation with StochasticSeriesExpansion.jl
* the maximum runtime before the simulation stops
* the interval between checkpoints
* the parameter sets (tasks) that should be simulated, created using `TaskMaker`

## Running the simulation
Next we should run it using
```bash
mpirun -n $NCORES julia bani2v2o8.jl run
```
While you wait, you can run `julia bani2v2o8.jl --help` for an explanation of the Carlo command line interface.

In general, the computational cost is linear in the number of spins and the inverse temperature. For the system defined here it is small enough to run on a powerful laptop. For larger systems, running on a high-performance computing cluster may be necessary.

Once the simulation is done, we will find a new directory `bani2v2o8.data` and a file `bani2v2o8.results.json`. The former contains the raw data that was recorded during the simulation and checkpoints, the latter has the postprocessed averages and errorbars of the observables.

!!! note "Checkpointing"
    Carlo.jl automatically creates checkpoints for simulations. This means if it is aborted, it will restart from where it left off. For a finished simulation it is also possible to increase the `sweeps` and continue it to gather more statistics.
    
    This means that to restart a simulation from scratch, it is necessary to first delete the job data,
    ```bash
    julia $JOBSCRIPT delete
    ```
    or to run with the restart flag
    ```bash
    julia $JOBSCRIPT run --restart
    ```

## Evaluating the results

The `Carlo.ResultTools` submodule has a tool to directly import the postprocessed data into a Julia dataframe. The following snippet is an example of that.

```@example
using Plots
using DataFrames
using LaTeXStrings
using Carlo.ResultTools

df = DataFrame(ResultTools.dataframe("bani2v2o8.results.json"))

df.L = [lattice["size"][1] for lattice in df.lattice]

plot(df.T, df.MagChi, group=df.L,
    xlabel = L"Temperature $T/J$",
    ylabel = L"Magnetic susceptibility $χ^z J$",
    legend_title = "L"
)
```
The dataframe contains all parameters from the job file as well as the observables that were measured. `MagChi` is the uniform magnetic susceptibility in ``S^z`` direction, perpendicular to the easy-plane anisotropy created by ``D_z``. As we can see, the basic shape of the susceptibility (Fig. 3(b) in the paper) is reproduced, but it is not converged to the thermodynamic limit at ``L=10``. We should therefore run larger system sizes until we can assure convergence (the data in the paper is at ``L=42``).

To get the in-plane susceptibility ``\chi^{x/y}``, we can simply rotate the anisotropy into the ``x`` direction (using the `Dx` parameter) and run a simulation measuring ``\chi^z`` again.
