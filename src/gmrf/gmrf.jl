"""
    GMRF(S::AbstractMatrix, κ::Real)

Construct a Gaussian Markov random field (GMRF)  with zero mean and precision matrix `Q =
κS`.

    GMRF(d, order::Integer, κ::Real, δ::Real; circular = false)

Construct a Gaussian Markov random field (GMRF) over domain `d`.

The precision matrix `Q` of the GMRF is determined as `κS`, with `κ` serving as the
scaling parameter. The structure matrix `S` is computed using `S = D'D + δI`, where `D`
represents the difference matrix, with its calculation dependent on the order specified
for the penalty or increments used to define the GMRF. The parameter `δ` is introduced to
address positive definiteness issues within the precision matrix. Notably, higher values
of `δ` result in a faster correlation decay, a higher `order` yields smoother processes,
and a greater `κ` leads to reduced process variance. This GMRF construction is applicable
to diverse domains, including `CartesianGrid`, `GeometrySet`, and `SimpleGraph`.
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
