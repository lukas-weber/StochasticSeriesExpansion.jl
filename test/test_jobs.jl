using LoadLeveller.JobTools


function generate_test_jobs(jobdir::AbstractString, sweeps::Integer, thermalization::Integer) where {Model <: S.AbstractModel}
    jobs = Tuple{::Type{<:S.AbstractModel}, JobInfo}[]

    function test_job(name::AbstractString, tm::TaskMaker)
        return JobInfo("$jobdir/$name", run_time = "15:00", checkpoint_time = "15:00", tasks = make_tasks(tm))
    end

    tm = TaskMaker()
    tm.sweep = sweeps
    tm.thermalization = thermalization
    tm.binsize = 100

    Ts = range(0.02,1.0,steps=7)

    tm.unitcell = S.UnitCells.square
    tm.Lx = 4
    tm.J = 1.24
    tm.hz = 0.1

    for T in Ts
        task(tm, T)
    

    push!(jobs, (S.Models.Magnet, test_job("magnet_square")))


    return jobs
end
    
