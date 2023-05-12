using LoadLeveller
using LoadLeveller.JobTools

function generate_test_jobs(
    jobdir::AbstractString,
    sweeps::Integer,
    thermalization::Integer,
)
    jobs = Dict{String,JobInfo}()

    function add_test_job(name::AbstractString, ::Type{Model}, tm::TaskMaker) where {Model}
        jobs[name] = JobInfo(
            "$jobdir/$name",
            S.mc(Model);
            run_time = "15:00",
            checkpoint_time = "15:00",
            tasks = make_tasks(tm),
        )

        return nothing
    end

    tm = TaskMaker()
    tm.sweeps = sweeps
    tm.thermalization = thermalization
    tm.binsize = 1000

    Ts = range(0.02, 4.00, length = 7)

    tm.unitcell = S.UnitCells.square
    tm.Lx = 2
    tm.Ly = 4

    tm.measure = [:magnetization]

    tm.J = 1.23
    tm.hz = -0.4

    tm.ed_run = 1

    for T in Ts
        task(tm; T = T)
    end

    add_test_job("magnet_square", S.Models.Magnet, tm)

    return jobs
end
