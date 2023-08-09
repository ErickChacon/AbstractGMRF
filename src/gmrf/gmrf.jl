"""
    GMRF(S, κ)

Construct a Gaussian Markov random field  with zero mean and precision matrix `Q = κS`.
"""
struct GMRF <: AbstractGMRF
    S::AbstractMatrix
    κ::Real
end

Base.length(x::GMRF) = size(x.S, 1)
structure(x::GMRF) = x.S
scale(x::GMRF) = x.κ
