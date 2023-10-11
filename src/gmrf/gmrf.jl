"""
    GMRF(S::AbstractMatrix, κ::Real)

Construct a Gaussian Markov random field  with zero mean and precision matrix `Q = κS`.
"""
struct GMRF <: AbstractGMRF
    S::AbstractMatrix
    κ::Real
end

# Constructors

GMRF(d::CartesianGrid, order::Integer, κ::Real, δ::Real; circular = false) =
    GMRF(structure(d; δ = δ, order = order, circular = circular), κ)

GMRF(d::SimpleGraph, order::Integer, κ::Real, δ::Real) =
    GMRF(structure(d; δ = δ, order = order), κ)

GMRF(d::GeometrySet, order::Integer, κ::Real, δ::Real) = GMRF(SimpleGraph(adjacency(d)), order, κ, δ)

# Methods

Base.length(x::GMRF) = size(x.S, 1)
scale(x::GMRF) = x.κ
structure(x::GMRF) = x.S
