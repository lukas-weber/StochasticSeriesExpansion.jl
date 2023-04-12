using LoadLeveller
using LoadLeveller.JobTools
using Random

@testset "SSE MC" begin
    sweeps = 1000
    params = Dict(
        :T => 0.1,
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

    state = copy(mc.state)

    for op in mc.operators
        sites = mc.sse_data.bonds[get_bond(op)].sites
        leg_states =
            @view get_vertex_data(sse_data, get_bond(op)).leg_states[:, get_vertex(op)]
        @test state[sites] .== leg_states[1:length(sites)]
        state[sites] .= leg_states[length(sites):end]
    end

    @test all(0 .< mc.state .< map(x -> x.dim, mc.sse_data.sites))
end
