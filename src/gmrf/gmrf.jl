"""
    GMRF(S::AbstractMatrix, κ::Real)

Construct a Gaussian Markov random field  with zero mean and precision matrix `Q = κS`.
"""
struct GMRF <: AbstractGMRF
    S::AbstractMatrix
    κ::Real
end

Base.length(x::GMRF) = size(x.S, 1)
scale(x::GMRF) = x.κ
structure(x::GMRF) = x.S
