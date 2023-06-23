@testset "VertexList" begin

    site_count = 4
    leg_count = 4
    bonds = [S.SSEBond(1, (1, 2)), S.SSEBond(1, (2, 3)), S.SSEBond(1, (1, 3))]

    v = S.VertexCode(false, 1)
    operators = [
        S.OperCode(S.Identity),
        S.OperCode(S.Identity),
        S.OperCode(1, v),
        S.OperCode(S.Identity),
        S.OperCode(1, v),
        S.OperCode(2, v),
        S.OperCode(3, v),
    ]

    expected = fill((-1, -1), leg_count, length(operators))
    expected[:, 3] = [(3, 7), (3, 6), (1, 5), (2, 5)]
    expected[:, 5] = [(3, 3), (4, 3), (1, 7), (1, 6)]
    expected[:, 6] = [(4, 5), (4, 7), (2, 3), (2, 7)]
    expected[:, 7] = [(3, 5), (4, 6), (1, 3), (2, 6)]

    vl = S.VertexList{leg_count ÷ 2}(site_count)
    S.make_vertex_list!(vl, operators, bonds)

    @test vl.v_last == [(3, 7), (3, 6), (4, 7), (-1, -1)]
    @test vl.v_first == [(1, 3), (2, 3), (2, 6), (-1, -1)]
    @test vl.vertices == expected
end
