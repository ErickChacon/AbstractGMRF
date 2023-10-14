@testset "Adjacency matrix" begin

    # 1D
    grid = CartesianGrid(10)
    gset = GeometrySet(elements(grid))
    graph = Graphs.SimpleGraphs.grid([10])
    graphc = Graphs.SimpleGraphs.grid([10], periodic = true)
    # non circular
    @test GMRFs.adjacency(grid) == GMRFs.adjacency(gset) == GMRFs.adjacency(graph) ==
        spdiagm(1 => repeat([true], 9)) |> xoxt
    # circular
    @test GMRFs.adjacency(grid, circular = true) == GMRFs.adjacency(graphc) ==
        spdiagm(1 => repeat([true], 9), 9 => [true]) |> xoxt

    # 2D
    grid = CartesianGrid(3,3)
    gset = GeometrySet(elements(grid))
    graph = Graphs.SimpleGraphs.grid([3,3])
    graphc = Graphs.SimpleGraphs.grid([3,3], periodic = true)
    # non circular
    @test GMRFs.adjacency(grid) == GMRFs.adjacency(graph) ==
        spdiagm(1 => sparse([1, 1, 0, 1, 1, 0, 1, 1] .== 1), 3 => repeat([true], 6)) |> xoxt
    # circular
    @test GMRFs.adjacency(grid, circular = true) == GMRFs.adjacency(graphc) ==
        spdiagm(
            1 => sparse([1, 1, 0, 1, 1, 0, 1, 1] .== 1),
            2 => sparse([1, 0, 0, 1, 0, 0, 1] .== 1),
            3 => repeat([true], 6),
            6 => repeat([true], 3)) |> xoxt

end

@testset "Difference matrix" begin

    # 1D
    grid = CartesianGrid(10)
    graph = Graphs.SimpleGraphs.grid([10])
    graphc = Graphs.SimpleGraphs.grid([10], periodic = true)
    # non circular
    @test GMRFs.difference(grid) == GMRFs.difference(graph) ==
        spdiagm(9, 10, 0 => repeat([-1], 9), 1 => repeat([1], 9))
    # circular
    @test GMRFs.difference(grid, circular = true) ==
        spdiagm(0 => repeat([-1], 10), 1 => repeat([1], 9), -9 => [1])
    @test all(sum(GMRFs.difference(graphc), dims = 2) .== 0)

    # 2D
    grid = CartesianGrid(3,3)
    graph = Graphs.SimpleGraphs.grid([3,3])
    graphc = Graphs.SimpleGraphs.grid([3,3], periodic = true)
    # non circular
    @test all(sum(GMRFs.difference(grid), dims = 2) .== 0)
    @test all(sum(GMRFs.difference(graph), dims = 2) .== 0)
    # circular
    @test all(sum(GMRFs.difference(grid, circular = true), dims = 2) .== 0)
    @test all(sum(GMRFs.difference(graphc), dims = 2) .== 0)

end

@testset "Structure matrix" begin

    # 1D
    grid = CartesianGrid(10)
    graph = Graphs.SimpleGraphs.grid([10])
    graphc = Graphs.SimpleGraphs.grid([10], periodic = true)
    # non circular
    @test structure(grid) == structure(graph) ==
        spdiagm(0 => vcat([1], repeat([2], 8), [1]), 1 => repeat([-1], 9), -1 => repeat([-1], 9))
    @test structure(grid, order = 2) ==
        xpxt(spdiagm(1 => vcat([-3], repeat([-4], 7), [-3]), 2 => repeat([1], 8))) +
            spdiagm(vcat([2], repeat([6], 8), [2]))
    @test structure(grid, order = 3) ==
        xpxt(spdiagm(1 => vcat([-9], repeat([-15], 7), [-9]), 2 => vcat([5], repeat([6], 6), [5]),
            3 => repeat([-1], 7))) +
            spdiagm(vcat([5,19], repeat([20], 6), [19,5]))
    # circular: structure_base
    @test GMRFs.structure_base(grid, order = 1) ==
        sparsevec([1, 2, 10], [2.0, -1, -1])
    @test GMRFs.structure_base(grid, order = 2) ==
        sparsevec([1, 2, 3, 9, 10], [6.0, -4, 1, 1, -4])
    @test GMRFs.structure_base(grid, order = 3) ==
        sparsevec([1, 2, 3, 4, 8, 9, 10], [20.0, -15, 6, -1, -1, 6, -15])
    # circular
    @test structure(grid, circular = true) == structure(graphc) ==
        xpxt(spdiagm(1 => repeat([-1], 9), 9 => [-1])) +
        spdiagm(repeat([2], 10))
    @test structure(grid, order = 2, circular = true) == structure(graphc, order = 2) ==
        xpxt(spdiagm(1 => repeat([-4], 9), 2 => repeat([1], 8), 8 => [1, 1], 9 => [-4])) +
            spdiagm(repeat([6], 10))
    @test structure(grid, order = 3, circular = true) == structure(graphc, order = 3) ==
        xpxt(spdiagm(1 => repeat([-15], 9), 2 => repeat([6], 8), 3 => repeat([-1], 7),
            7 => [-1, -1, -1], 8 => [6, 6], 9 => [-15])) +
            spdiagm(repeat([20], 10))

    # 2D
    grid = CartesianGrid(7,7)
    graph = Graphs.SimpleGraphs.grid([7,7])
    graphc = Graphs.SimpleGraphs.grid([7,7], periodic = true)
    # non circular
    A = GMRFs.adjacency(grid)
    @test structure(grid) == structure(graph) ==
        spdiagm(vec(sum(A, dims = 2))) - A
    # circular: structure_base
    @test GMRFs.structure_base(grid, order = 1) ==
        xpxt(sparse([1, 1, 1], [1, 2, 7], [2.0, -1, -1], 7, 7))
    @test GMRFs.structure_base(grid, order = 2) ==
        xpxt(sparse([1, 1, 1, 1, 1, 2, 2, 7], [1, 2, 3, 6, 7, 2, 7, 7], [10.0, -8, 1, 1, -8, 1, 2, 1]))
    @test GMRFs.structure_base(grid, order = 3) ==
        xpxt(sparse([1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 6, 7],
        [1, 2, 3, 4, 5, 6, 7, 2, 3, 6, 7, 7, 7, 7],
        [56.0, -57, 12, -1, -1, 12, -57, 12, -3, -3, 24, -3, -3, 12]))
    # circular: structure
    S = GMRFs.structure(grid, order = 1, circular = true)
    @test S == GMRFs.structure(graphc, order = 1)
    @test S[1, :] == vec(GMRFs.structure_base(grid, order = 1)')
    S = GMRFs.structure(grid, order = 2, circular = true)
    @test S == GMRFs.structure(graphc, order = 2)
    @test S[1, :] == vec(GMRFs.structure_base(grid, order = 2)')
    S = GMRFs.structure(grid, order = 3, circular = true)
    @test S == GMRFs.structure(graphc, order = 3)
    @test S[1, :] == vec(GMRFs.structure_base(grid, order = 3)')

end
