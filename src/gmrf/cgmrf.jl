"""
    CGMRF(g::CartesianGrid, order::Integer, δ::Real, κ::Real)

Circulant Gaussian Markov random field over the regular grid `g` with scale parameter `κ` and structure
matrix defined by the penalty of n-th `order` and `δ` parameter.
"""
struct CGMRF <: AbstractGMRF
    g::CartesianGrid
    order::Integer
    δ::Number
    κ::Number
end

Base.length(d::CGMRF) = nelements(d.g)
scale(d::CGMRF) = d.κ
structure(d::CGMRF) = structure(d.g; δ = d.δ, order = d.order, circular = true)
structure_base(d::CGMRF) = structure_base(d.g; δ = d.δ, order = d.order)


## Random generator

function Distributions._rand!(rng::AbstractRNG, d::CGMRF, x::AbstractVector{T}) where T<:Real
    # get dimensions
    dims = size(d.g)
    n = prod(dims)

    # compute eigenvalues
    base = structure_base(d)
    λ = FFTW.fft(base)

    # simulate using fft
    z = randn(rng, reverse(dims)) + randn(rng, reverse(dims)) * im
    xaux = real(FFTW.fft(λ .^ (-1/2) .* z))
    if length(dims) == 2
        xaux = transpose(reverse(xaux, dims = 1))
    end
    copyto!(x, xaux)
    ldiv!(sqrt(n * scale(d)), x)

    return x
end

function Distributions._rand!(rng::AbstractRNG, d::CGMRF, x::AbstractArray{<:Real})
    dims = size(d.g)
    n = prod(dims)
    λ = FFTW.fft(structure_base(d))
    @inbounds for xi in Distributions.eachvariate(x, Distributions.variate_form(typeof(d)))
        z = randn(rng, reverse(dims)) + randn(rng, reverse(dims)) * im
        xaux = real(FFTW.fft(λ .^ (-1/2) .* z))
        if length(dims) == 2
            xaux = transpose(reverse(xaux, dims = 1))
        end
        copyto!(xi, xaux)
        ldiv!(sqrt(n * scale(d)), xi)
    end
    return x
end

## Logarithm of the pdf

# function Distributions._logpdf(X::CGMRF, x::AbstractVector{T}) where T<:Real
#     # get dimensions
#     dims = size(X.g)
#     n = prod(dims)
#
#     # compute eigenvalues
#     base = structure_base(X.g; δ = X.δ, order = X.order)
#     λ = fft(base)
#
#     # output = (- n * log(2π) + sum(log.(real(λ)))) / 2
#     output = λ .* ifft()
# end
