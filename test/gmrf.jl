@testset "GMRF 1D" begin

    n = 100
    order = [1, 2, 3]
    kappa = [8, 5, 25]
    δ = [0.01, 0.01, 0.001]
    grid = CartesianGrid(n)

    for i = 1:length(order)
        X = GMRF(grid, order[i], kappa[i], δ[i])
        mvn = MvNormalCanon(Matrix(GMRFs.precision(X)))
        # basic methods
        @test length(X) == n
        @test GMRFs.scale(X) == kappa[i]
        @test GMRFs.structure(X) == GMRFs.structure(grid, δ = δ[i], order = order[i])
        @test GMRFs.precision(X) == structure(grid, δ = δ[i], order = order[i]) * GMRFs.scale(X)
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

@testset "GMRF 2D" begin

    n1, n2 = 10, 13
    order = [1, 2, 3]
    kappa = [1, 1, 7]
    δ = [0.01, 0.01, 0.001]
    grid = CartesianGrid(n1,n2)

    for i = 1:length(order)
        X = GMRF(grid, order[i], kappa[i], δ[i])
        mvn = MvNormalCanon(Matrix(GMRFs.precision(X)))
        # basic methods
        @test length(X) == n1 * n2
        @test GMRFs.scale(X) == kappa[i]
        @test GMRFs.structure(X) == GMRFs.structure(grid, δ = δ[i], order = order[i])
        @test GMRFs.precision(X) == structure(grid, δ = δ[i], order = order[i]) * GMRFs.scale(X)
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

