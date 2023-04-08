@testset "OperCode" begin
    for bond in rand(1:100*100, 20)
        for vidx in rand(1:77, 20)
            for diagonal in [true, false]
                v = S.VertexCode(diagonal, vidx)

                @test !S.isinvalid(v)

                @test S.isdiagonal(v) == diagonal
                @test S.get_vertex_idx(v) == vidx

                op = S.OperCode(bond, v)

                @test S.get_bond(op) == bond
                @test S.isdiagonal(op) == diagonal
                @test S.get_vertex(op).code == v.code

                @test !S.isidentity(op)
            end
        end
    end

    @test S.isidentity(S.OperCode(S.Identity))
end
