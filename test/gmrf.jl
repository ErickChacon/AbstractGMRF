@testset "GMRF 1D" begin

    n = 10
    order = [1, 2]
    kappa = [100, 100]
    δ = [0.0001, 1]

    for i = 1:length(order)
        grid = CartesianGrid(n)
        X = GMRF(grid, order[i], kappa[i], δ[i])
        mvn = MvNormalCanon(Matrix(GMRFs.precision(X)))
        # random realizations
        x = rand(X)
        @test length(x) == n
        x = rand(X, 3)
        @test size(x) == (n, 3)
        # logpdf
        lpdf = logpdf(X, x)
        @test length(lpdf) == 3
        @test isapprox(lpdf, logpdf(mvn, x))
    end
end

@testset "GMRF 2D" begin

    n = 7
    order = [1, 2]
    kappa = [100, 100]
    δ = [0.0001, 1]

    for i = 1:length(order)
        grid = CartesianGrid(n,n)
        X = GMRF(grid, order[i], kappa[i], δ[i])
        mvn = MvNormalCanon(Matrix(GMRFs.precision(X)))
        # random realizations
        x = rand(X)
        @test length(x) == n^2
        x = rand(X, 3)
        @test size(x) == (n^2, 3)
        # logpdf
        lpdf = logpdf(X, x)
        @test length(lpdf) == 3
        @test isapprox(lpdf, logpdf(mvn, x))
    end

end

    # gset = GeometrySet(elements(grid))
    # graph = Graphs.SimpleGraphs.grid([10])
    # graphc = Graphs.SimpleGraphs.grid([10], periodic = true)
    # # non circular
    # @test GMRFs.adjacency(grid) == GMRFs.adjacency(gset) == GMRFs.adjacency(graph) ==
    #     spdiagm(1 => repeat([true], 9)) |> xoxt
    # # circular
    # @test GMRFs.adjacency(grid, circular = true) == GMRFs.adjacency(graphc) ==
    #     spdiagm(1 => repeat([true], 9), 9 => [true]) |> xoxt

