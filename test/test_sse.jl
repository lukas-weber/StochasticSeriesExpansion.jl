using LoadLeveller
using LoadLeveller.JobTools
using Random

function isconsistent(
    operators::AbstractVector{<:S.OperCode},
    state0::AbstractVector,
    sse_data::S.SSEData,
)
    state = copy(state0)
    for op in operators
        if S.isidentity(op)
            continue
        end
        sites = sse_data.bonds[S.get_bond(op)].sites
        leg_states =
            S.get_leg_state(S.get_vertex_data(sse_data, S.get_bond(op)), S.get_vertex(op))
        if state[collect(sites)] != leg_states[1:length(sites)]
            return false
        end
        state[collect(sites)] .= leg_states[length(sites)+1:end]
    end

    return true
end

@testset "worm_traverse!" begin
    params = Dict(
        :T => 0.1,
        :unitcell => S.UnitCells.square,
        :Lx => 4,
        :Ly => 4,
        :J => 1.0,
        :S => 1 // 2,
    )

    sse_data = S.generate_sse_data(S.Models.Magnet(params))

    vl = S.VertexList{S.leg_count(S.Models.Magnet) ÷ 2}(length(sse_data.sites))
    operatorss = [
        [S.OperCode(1, S.get_vertex_data(sse_data, 1).diagonal_vertices[3])],
        [S.OperCode(1, S.VertexCode(true, 1)), S.OperCode(1, S.VertexCode(true, 1))],
    ]
    for operators in operatorss
        S.make_vertex_list!(vl, operators, sse_data.bonds)
        S.worm_traverse!((1, 1, 1), operators, vl.vertices, sse_data, Random.Xoshiro())

        state0 = [
            vf[1] < 0 ? 1 :
            S.get_leg_state(
                S.get_vertex_data(sse_data, S.get_bond(operators[vf[2]])),
                S.get_vertex(operators[vf[2]]),
            )[vf[1]] for vf in eachcol(vl.v_first)
        ]
        @test isconsistent(operators, state0, sse_data)
    end
end


@testset "SSE MC" begin
    sweeps = 1000
    params = Dict(
        :T => 0.1,
        :unitcell => S.UnitCells.honeycomb,
        :Lx => 4,
        :Ly => 4,
        :J => 1.4,
        :S => 1 // 2,
        :binsize => 1,
        :sweeps => sweeps,
        :thermalization => 100,
    )

    ctx = MCContext{Random.Xoshiro}(params)

    mc = S.mc(S.Models.Magnet)(params)

    LoadLeveller.init!(mc, ctx, params)

    for s = 1:sweeps
        if !isconsistent(mc.operators, mc.state, mc.sse_data)
            break
        end
        LoadLeveller.sweep!(mc, ctx)
    end

    @test isconsistent(mc.operators, mc.state, mc.sse_data)
    @test all(0 .< mc.state .<= map(x -> x.dim, mc.sse_data.sites))
end
