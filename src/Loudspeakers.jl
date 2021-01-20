module Loudspeakers

export fillinsert, rfft, irfft
export normdeg, deltadeg
export create_fir_filter, resamplefreq, zerophase, rms_mean

import Base: show, lpad, rpad, getindex, length, lastindex
using DataFrames
using Statistics
import Statistics: mean
using StatsBase
using Reexport
using CSV
using DSP
using FFTW
import FFTW: rfft, irfft
using LinearAlgebra 
using RecipesBase
using IntervalSets
using Printf
using ProgressBars
using StatsBase
using WAV

using Interpolations
import Interpolations.interpolate
export interpolate

using FileIO: load, save, loadstreaming, savestreaming
using Unitful

using SampleArrays

include("utils.jl")

include("smoothing/octave.jl")

include("measurements/peak.jl")

include("measurements/directivity.jl")

include("measurements/spinorama.jl")

include("measurements/impexp.jl")

include("measurements/quasianechoic.jl")

include("measurements/expsweep.jl")

include("measurements/measure.jl")

include("vis/timedomain.jl")

include("vis/freqdomain.jl")

include("vis/directivity.jl")

# TOOLS

"""
    normdeg(α)

Normalize α in degrees to (-180, 180].
"""
function normdeg(α)
    α = mod(α, 360)
    α > 180 ? α - 360 : α
end

"""
    deltadeg(α, β)

Smallest delta degrees between α and β angles, sing positive for α to β CCW.
"""
deltadeg(α, β) = normdeg(β - α)

# CANERASE
# function Base.show(io::IO, ::MIME"text/html", buf::SampleBuf{T, N}) where {T <: Number, N}
#     show(io, MIME"text/html", buf)
# end

"""
    fillinsert(X::AbstractArray{T, 2}, n, idx, v=0, dim=1)

Insert into array `X`'s dimension `d` adding `n` to that dimension size; start after index `idx`, fill with `v`. 
"""
function fillinsert(X::AbstractArray{T}, n, idx, v=zero(T), dim=1) where T
    idx > size(X, dim) && throw(ArgumentError("idx ($idx) > size(X, $dim) ($(size(X, dim)))!"))
    1 <= dim <= ndims(X) || throw(ArgumentError("wrong dim ($dim) should be in {1, ..., $(ndims(X))}"))
    tdims = collect(size(X))
    tdims[dim] += n
    
    t = similar(X, tdims...)
    idxs(sel, dim=dim, ndims=ndims(X)) = [i == dim ? sel : Colon() for i in 1:ndims]
    t[idxs(1:idx)...] .= X[idxs(1:idx)...]
    t[idxs(idx+1:idx+n)...] .= v
    t[idxs(idx+n+1:tdims[dim])...] .= X[idxs(idx+1:size(X, dim))...]
#     @show X[idxs(idx+1:end)...]
    t
end

"""
   lpad(X::AbstractArray{T}, n::Integer, v=0, dim=1)

Prepend array `X` in dimension `dim` with `n` values `v`.
"""
Base.lpad(X::AbstractArray{T}, n::Integer, v=zero(T), dim=1) where T = fillinsert(X, n, 0, v, dim) 

"""
   rpad(X::AbstractArray{T}, n::Integer, v=0, dim=1)

Append `n` values `v` to array `X` in dimension `dim`.
"""
Base.rpad(X::AbstractArray{T}, n::Integer, v=zero(T), dim=1) where T = fillinsert(X, n, size(X, 1), v, dim) 

function create_fir_filter(X::RFFTSpectrumArray, taps)
    @assert nchannels(X) == 1
    x = irfft(data(X), 2nfreqs(X)-1)
    x .= circshift(x, taps÷2)
    xsa = SampleArray(x[1:taps], rate(X)) .* hanning(taps; zerophase=false)
    xsa, rfft(xsa)
end

# TODO: move these to SampleArrays?
Statistics.mean(X::AbstractSpectrumArray{<:SignalElement}) = _mean(X, uweights(Float32, nchannels(X)))
Statistics.mean(X::AbstractSpectrumArray{<:SignalElement}, w::AbstractWeights) = _mean(X, w)

function _mean(X::AbstractSpectrumArray{<:SignalElement}, w::AbstractWeights)
    # complex mean over channels according to
    # Panzer: The use of continuous phase for interpolation, smoothing and forming mean values of complex frequency response curves, 2004
    # println("_mean()")
    # Xu = unwrap(X) # original
    Xu = X # NO UNWRAP
    mags = mean(abs.(Xu), w, dims=2) # original
    
    # phis = mean(angle.(Xu), w, dims=2) # original
    # phis = mean(angle.(Xu), dims=2) # NO WEIGHTING PHASE
    # phis = angle.(Xu)[:, argmax(w)] # PHASE OF A CLOSER
    phis = angle.(mean(data(X), w, dims=2)) # PHASE MEAN based on vector mean
    # phis = angle.(Xu)[:, 1]

    Y = similar(X, nframes(X), 1)
    Y .= (MagPhase(m, p) for (m, p) in zip(mags, phis))
    Y
end

Statistics.mean(X::AbstractSpectrumArray) = _mean(X, uweights(Float32, nchannels(X)))
Statistics.mean(X::AbstractSpectrumArray, w::AbstractWeights) = _mean(X, w)

_mean(X::AbstractSpectrumArray, w::AbstractWeights) = mean(X, w, dims=2)


rms_mean(X::AbstractSpectrumArray{<:SignalElement}) = _rms_mean(X, uweights(Float32, nchannels(X)))
rms_mean(X::AbstractSpectrumArray{<:SignalElement}, w::AbstractWeights) = _rms_mean(X, w)

function _rms_mean(X::AbstractSpectrumArray{<:SignalElement}, w::AbstractWeights)
    error("FIX phase!")
    # RMSE mean magnitude, unwrapped power
    # println("_rms_mean()")
    Xu = unwrap(X)
    # Xu = X
    ms = abs.(Xu)
    mags = sqrt.(mean(ms .* ms, w, dims=2))
    phis = mean(angle.(Xu), w, dims=2)
    Y = similar(X, nframes(X), 1)
    Y .= (MagPhase(m, p) for (m, p) in zip(mags, phis))
    Y
end

end # module
