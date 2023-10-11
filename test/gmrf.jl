@testset "GMRF 1D" begin

    n = 10
    order = [1, 2]
    kappa = [100, 100]
    δ = [0.0001, 1]
    grid = CartesianGrid(n)

    for i = 1:length(order)
        X = GMRF(grid, order[i], kappa[i], δ[i])
        # basic methods
        @test length(X) == n
        @test GMRFs.scale(X) == kappa[i]
        @test GMRFs.structure(X) == GMRFs.structure(grid, δ = δ[i], order = order[i])
        @test GMRFs.precision(X) == structure(grid, δ = δ[i], order = order[i]) * GMRFs.scale(X)
        # rand
        x = rand(X)
        @test length(x) == n
        x = rand(X, 3)
        @test size(x) == (n, 3)
        # logpdf
        lpdf = logpdf(X, x)
        @test length(lpdf) == 3
        mvn = MvNormalCanon(Matrix(GMRFs.precision(X)))
        @test isapprox(lpdf, logpdf(mvn, x))
    end
end

@testset "GMRF 2D" begin

    n1, n2 = 7, 9
    order = [1, 2]
    kappa = [100, 100]
    δ = [0.0001, 1]
    grid = CartesianGrid(n1,n2)

    for i = 1:length(order)
        X = GMRF(grid, order[i], kappa[i], δ[i])
        # basic methods
        @test length(X) == n1 * n2
        @test GMRFs.scale(X) == kappa[i]
        @test GMRFs.structure(X) == GMRFs.structure(grid, δ = δ[i], order = order[i])
        @test GMRFs.precision(X) == structure(grid, δ = δ[i], order = order[i]) * GMRFs.scale(X)
        # random realizations
        x = rand(X)
        @test length(x) == n1 * n2
        x = rand(X, 3)
        @test size(x) == (n1 * n2, 3)
        # logpdf
        lpdf = logpdf(X, x)
        @test length(lpdf) == 3
        mvn = MvNormalCanon(Matrix(GMRFs.precision(X)))
        @test isapprox(lpdf, logpdf(mvn, x))
    end

end

