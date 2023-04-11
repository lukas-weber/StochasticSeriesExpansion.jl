using LoadLeveller
using LoadLeveller.JobTools
using Random

@testset "SSE MC" begin
    sweeps = 1000
    params = Dict(
        :unitcell => S.UnitCells.honeycomb,
        :Lx => 4,
        :Ly => 4,
        :J => 1.4,
        :S => 1,
        :binsize => 1,
        :sweeps => sweeps,
        :thermalization => 100,
    )

    ctx = MCContext{Random.Xoshiro}(params)

    mc = S.MC{S.Models.Magnet}(params)

    LoadLeveller.init!(mc, ctx, params)

    for s = 1:sweeps
        LoadLeveller.sweep!(mc, ctx)
    end

    @test all(0 .< mc.state .< map(x -> x.dim, mc.sse_data.sites))
end
