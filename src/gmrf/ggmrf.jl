"""
    GGMRF(grid, order, δ, κ)

Regular Gaussian Markov random field with n-th order, δ as and κ as precision.
"""
struct GGMRF <: AbstractGMRF
    g::SimpleGraph
    δ::Number
    κ::Real
end

Base.length(d::GGMRF) = length(Graphs.vertices(g))
structure(d::GGMRF) = structure(d.grid; δ = d.δ, order = d.order)
scale(d::GGMRF) = d.κ
