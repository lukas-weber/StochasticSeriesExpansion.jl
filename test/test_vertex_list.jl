@testset "VertexList" begin

    site_count = 4
    leg_count = 4
    bonds = [
        S.Bond(1, (1,2)),
        S.Bond(1, (2,3)),
        S.Bond(1, (1,3)),
    ]

    v = S.VertexCode(false, 1)
    operators = [
        S.OperCode(S.Identity),
        S.OperCode(S.Identity),
        S.OperCode(1,v),
        S.OperCode(S.Identity),
        S.OperCode(1,v),
        S.OperCode(2,v),
        S.OperCode(3,v),
    ]

    expected = fill(-1, 2, leg_count, length(operators))
    expected[:, :, 3] = [
        3 3 1 2
        7 6 5 5
    ]
    expected[:, :, 5] = [
        3 4 1 1
        3 3 7 6
    ]
    expected[:, :, 6] = [
        4 4 2 2
        5 7 3 7
    ]
    expected[:, :, 7] = [
        3 4 1 2
        5 6 3 6
    ]

    vl = S.VertexList{leg_count÷2}(site_count)
    S.make_vertex_list!(vl, operators, bonds)

    @test vl.v_last == [
        3 3 4 -1
        7 6 7 -1
    ]
    @test vl.v_first == [
        1 2 2 -1
        3 3 6 -1
    ]
    @test vl.vertices == expected
end
