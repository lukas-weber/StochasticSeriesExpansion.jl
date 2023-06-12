function test_detailed_balance(vd::S.VertexData)
    for (v, weight) in enumerate(vd.weights)
        for step in CartesianIndices(@view vd.transitions[:, :, v])
            (leg_in, worm_in) = Tuple(step)
            t = vd.transitions[leg_in, worm_in, v]

            if S.isinvalid(t)
                continue
            end

            probs =
                diff(vcat([0.0], @view vd.transition_cumprobs[t.offset:t.offset+t.length]))
            targets = @view vd.transition_targets[t.offset:t.offset+t.length]
            step_outs = @view vd.transition_step_outs[t.offset:t.offset+t.length]


            for (prob, target, step_out) in zip(probs, targets, step_outs)
                #@test S.vertex_apply_change(vd.leg_states, vd.dims, v, (leg_in, worm_in), step_out) == Int(S.get_vertex_idx(target))
                (leg_out, worm_out) = step_out

                t_inv = vd.transitions[
                    leg_out,
                    S.worm_inverse(
                        worm_out,
                        vd.dims[S.site_of_leg(leg_out, length(vd.dims))],
                    ),
                    S.get_vertex_idx(target),
                ]
                i_inv = findfirst(
                    v .==
                    S.get_vertex_idx.(
                        @view vd.transition_targets[t_inv.offset:t_inv.offset+t_inv.length]
                    ),
                )

                @test i_inv !== nothing
                if i_inv === nothing
                    println("$v, $step -> $(S.get_vertex_idx(target)), $step_out")
                    println("$(vd.leg_states[:,v]) -> $(S.get_leg_state(vd, target))")
                    return nothing
                end

                weight_new = S.get_vertex_weight(vd, target)
                probs_back = diff(
                    vcat(
                        [0.0],
                        @view vd.transition_cumprobs[t_inv.offset:t_inv.offset+t_inv.length]
                    ),
                )
                prob_back = probs_back[i_inv]
                @test prob * weight ≈ prob_back * weight_new
            end
        end
    end

    return nothing
end


@testset "VertexData" begin
    @test S.isinvalid(S.Transition())

    @testset "vertex_change_apply" begin
        dim = 4
        nsites = 2
        leg_states =
            reduce(hcat, collect.(Iterators.product(fill(UInt8.(1:dim), 2 * nsites)...)))

        for step_in in Iterators.product(1:2*nsites, 1:S.worm_count(dim))
            v = rand(1:dim^(2*nsites))
            @test S.vertex_apply_change(
                leg_states,
                Tuple(dim for i = 1:2*nsites),
                v,
                step_in,
                (step_in[1], S.worm_inverse(step_in[2], dim)),
            ) == v
        end
    end

    @testset "S=1/2 Heisenberg" begin
        (splus, sz) = S.spin_operators(2)
        Hbond = kron(sz, sz) + 0.5 * (kron(splus, splus') + kron(splus', splus))

        vd = S.VertexData((2, 2), Hbond; energy_offset_factor = 0.0)

        test_detailed_balance(vd)

        @test vd.energy_offset ≈ -0.25
        @test sum(S.isinvalid.(vd.diagonal_vertices)) == 2
        @test all(vd.weights .≈ 0.5)
        @test vd.transition_cumprobs ≈ ones(4 * 4)

        for vertex = 1:length(vd.weights)
            for leg = 1:4
                @test vd.transition_step_outs[vd.transitions[leg, 1, vertex].offset][1] ==
                      xor(leg - 1, 1) + 1
            end
        end
    end

    @testset "Random Hamiltonian" begin
        dimss = [(4, 4), (2, 4)]

        for dims in dimss
            @testset "dims = $(dims)" begin
                Hbond = rand(prod(dims), prod(dims))
                vd = S.VertexData(dims, Hbond)
                test_detailed_balance(vd)

                for v = 1:size(vd.transitions, 3)
                    leg_state = @view vd.leg_states[:, v]
                    @test all(leg_state .<= tuple(dims..., dims...))
                    if !all(leg_state .<= tuple(dims..., dims...))
                        @show leg_state, dims
                    end
                end

                for state = 1:prod(dims)
                    split = S.split_idx(dims, state)
                    diag = S.get_leg_state(vd, S.get_diagonal_vertex(vd, state))
                    @test diag == UInt8[split..., split...]
                end

                for I in CartesianIndices(vd.transitions)
                    if I[2] > S.worm_count(dims[S.site_of_leg(I[1], length(dims))])
                        @test S.isinvalid(vd.transitions[I])
                    end
                end

                for t in vd.transitions
                    @test S.isinvalid(t) || vd.transition_cumprobs[t.offset+t.length] ≈ 1.0
                end
            end
        end
    end
end
