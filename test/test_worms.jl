@testset "Worms" begin
    for dim in [2, 4, 8]
        @testset "dim = $(dim)" begin
            state_ins = collect(1:dim)

            for state_in in state_ins
                state_outs = Int[]

                for worm = 1:S.worm_count(dim)
                    state_out = S.worm_action(worm, state_in, dim)
                    @test S.worm_action(S.worm_inverse(worm, dim), state_out, dim) ==
                          state_in

                    push!(state_outs, state_out)
                end

                push!(state_outs, state_in)

                sort!(state_outs)
                @test state_outs == state_ins
            end
        end
    end
end
