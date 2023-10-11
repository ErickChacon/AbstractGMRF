@testset "CGMRF 1D" begin

    n = 100
    order = [1, 2, 3]
    kappa = [5, 11, 105]
    δ = [0.01, 0.01, 0.001]
    grid = CartesianGrid(n)

    for i = 1:length(order)
        X = CGMRF(grid, order[i], kappa[i], δ[i])
        Q = GMRFs.structure(grid, δ = δ[i], order = order[i], circular = true) * GMRFs.scale(X)
        mvn = MvNormalCanon(Matrix(Q))
        # basic methods
        @test length(X) == n
        @test GMRFs.scale(X) == kappa[i]
        @test GMRFs.structure_base(X) == GMRFs.structure_base(grid, order = order[i], δ = δ[i])
        # rand and logpdf: single
        x = rand(X)
        @test length(x) == n
        @test isapprox(logpdf(X, x), logpdf(mvn, x))
        # rand and logpdf: multiple
        x = rand(X, 3)
        @test size(x) == (n, 3)
        @test isapprox(logpdf(X, x), logpdf(mvn, x))
    end
end

@testset "CGMRF 2D" begin

    n1, n2 = 10, 13
    order = [1, 2, 3]
    kappa = [1, 1, 10]
    δ = [0.01, 0.01, 0.001]
    grid = CartesianGrid(n1,n2)

    for i = 1:length(order)
        X = CGMRF(grid, order[i], kappa[i], δ[i])
        Q = GMRFs.structure(grid, δ = δ[i], order = order[i], circular = true) * GMRFs.scale(X)
        mvn = MvNormalCanon(Matrix(Q))
        # basic methods
        @test length(X) == n1 * n2
        @test GMRFs.scale(X) == kappa[i]
        @test GMRFs.structure_base(X) == GMRFs.structure_base(grid, order = order[i], δ = δ[i])
        # rand and logpdf: single
        x = rand(X)
        @test length(x) == n1 * n2
        @test isapprox(logpdf(X, x), logpdf(mvn, x))
        # rand and logpdf: multiple
        x = rand(X, 3)
        @test size(x) == (n1 * n2, 3)
        @test isapprox(logpdf(X, x), logpdf(mvn, x))
    end

end

