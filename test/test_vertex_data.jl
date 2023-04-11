
@testset "VertexData" begin
    @test S.isinvalid(S.Transition())

    @testset "S=1/2 Heisenberg" begin
        (splus, sz) = S.spin_operators(2)
        Hbond = kron(sz, sz) + 0.5 * (kron(splus, splus') + kron(splus', splus))

        vd = S.VertexData((2, 2), Hbond; energy_offset_factor = 0.0)

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
end
