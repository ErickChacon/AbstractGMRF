"""
    GGMRF(g::SimpleGraph, order::Integer, δ::Real, κ::Real)

Gaussian Markov random field over the graph `g` with scale parameter `κ` and structure
matrix defined by the penalty of n-th `order` and `δ` parameter.
"""
struct GGMRF <: AbstractGMRF
    g::SimpleGraph
    order::Integer
    δ::Real
    κ::Real
end

GGMRF(d::GeometrySet, order::Integer, δ::Real, κ::Real) = GGMRF(SimpleGraph(adjacency(d)), order, δ, κ)

Base.length(x::GGMRF) = Graphs.nv(x.g)
structure(x::GGMRF) = structure(x.g; δ = x.δ, order = x.order)
scale(x::GGMRF) = x.κ
