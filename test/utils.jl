@testset "Adjacency" begin
    # 1D - adjacency
    grid = CartesianGrid(10)
    gset = GeometrySet(elements(grid))
    graph = Graphs.SimpleGraphs.grid([10])
    graphc = Graphs.SimpleGraphs.grid([10], periodic = true)
    # non circular
    @test adjacency(grid) == adjacency(gset) == adjacency(graph) ==
        spdiagm(1 => repeat([true], 9)) |> xxt
    # circular
    @test adjacency(grid, circular = true) == adjacency(graphc) ==
        spdiagm(1 => repeat([true], 9), 9 => [true]) |> xxt

    # 2D - adjacency
    grid = CartesianGrid(3,3)
    gset = GeometrySet(elements(grid))
    graph = Graphs.SimpleGraphs.grid([3,3])
    graphc = Graphs.SimpleGraphs.grid([3,3], periodic = true)
    # non circular
    @test adjacency(grid) == adjacency(graph) ==
        spdiagm(1 => sparse([1, 1, 0, 1, 1, 0, 1, 1] .== 1), 3 => repeat([true], 6)) |> xxt
    # circular
    @test adjacency(grid, circular = true) == adjacency(graphc) ==
        spdiagm(
            1 => sparse([1, 1, 0, 1, 1, 0, 1, 1] .== 1),
            2 => sparse([1, 0, 0, 1, 0, 0, 1] .== 1),
            3 => repeat([true], 6),
            6 => repeat([true], 3)) |> xxt
end

# @testset "Difference and structure matrices" begin
#     # 1d
#     grid = CartesianGrid(Point(-10.0), Point(10.0), dims = (5,))
#     # non circular difference
#     @test difference(grid, order = 1, circular = false) ==
#         spdiagm(4, 5, 0 => repeat([-1], 4), 1 => repeat([1], 4))
#     @test difference(grid, order = 2, circular = false) ==
#         spdiagm(3, 5, 0 => repeat([1], 3), 1 => repeat([-2], 3), 2 => repeat([1], 3))
#     # circular difference
#     @test difference(grid, order = 1, circular = true) ==
#         spdiagm(0 => repeat([-1], 5), 1 => repeat([1], 4), -4 => [1])
#     @test difference(grid, order = 2, circular = true) ==
#         spdiagm(0 => repeat([1], 5), 1 => repeat([-2], 4), 2 => repeat([1], 3),
#         -3 => repeat([1], 2), -4 => [-2])
#
#     # 2d
#     grid = CartesianGrid(Point(-10.0, -10.0), Point(10.0, 10.0), dims = (3,3))
#     # non circular difference
#     # circular difference
# end
#
