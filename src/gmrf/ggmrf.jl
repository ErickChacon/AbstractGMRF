"""
    GGMRF(grid, order, δ, κ)

Regular Gaussian Markov random field with n-th order, δ as and κ as precision.
"""
struct GGMRF <: AbstractGMRF
    g::SimpleGraph
    order::Integer
    δ::Number
    κ::Real
end

GGMRF(m::GeometrySet, order::Integer, δ::Number, κ::Real) = GGMRF(SimpleGraph(adjacency(m)), order, δ, κ)

Base.length(d::GGMRF) = Graphs.nv(d.g)
structure(d::GGMRF) = structure(d.g; δ = d.δ, order = d.order)
scale(d::GGMRF) = d.κ
