
@testset "split_idx" begin
    @test S.split_idx((5, 2), 2) == (2, 1)
    @test S.split_idx((5, 2), 6) == (1, 2)

    dims = (6, 3, 1, 5)
    for i = 1:prod(dims)
        @test S.join_idx(dims, S.split_idx(dims, i)) == i
    end

    for dims in eachcol(rand(1:400, 3, 10))
        for i in rand(1:prod(dims), 10)
            @test S.join_idx(Tuple(dims), S.split_idx(Tuple(dims), i)) == i
        end
    end
end
