"""
    GMRF(S::AbstractMatrix, κ::Real)

Construct a Gaussian Markov random field  with zero mean and precision matrix `Q = κS`.
"""
struct GMRF <: AbstractGMRF
    S::AbstractMatrix
    κ::Real
end

GMRF(g::SimpleGraph, order::Integer, δ::Real, κ::Real) = GMRF(structure(g; δ = δ, order = order), κ)

GMRF(g::CartesianGrid, order::Integer, δ::Real, κ::Real; circular = false) =
    GMRF(structure(g; δ = δ, order = order, circular = circular), κ)

Base.length(x::GMRF) = size(x.S, 1)
scale(x::GMRF) = x.κ
structure(x::GMRF) = x.S
