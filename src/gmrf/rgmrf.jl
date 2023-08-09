"""
    RGMRF(g::CartesianGrid, order::Integer, δ::Real, κ::Real)

Regular Gaussian Markov random field over regular grid `g` with scale parameter `κ` and structure
matrix defined by the penalty of n-th `order` and `δ` parameter.
"""
struct RGMRF <: AbstractGMRF
    g::CartesianGrid
    order::Integer
    δ::Real
    κ::Real
end

Base.length(x::RGMRF) = nelements(x.g)
scale(x::RGMRF) = x.κ
structure(x::RGMRF) = structure(x.g; δ = x.δ, order = x.order, cyclic = false)
