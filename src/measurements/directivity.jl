export DirectivitySpectra, spectra, angles, maximumamp, minimumamp, normalizeamp, normalizeamp!
export normalizeampangle, normalizeampangle!
export normalizeampfreq, normalizeampfreq!
export relativeto, relativeto!, smooth, closestangle, between, interpolate, mirror, mirror!

struct DirectivitySpectra{T<:AbstractSpectrumArray}
    spectra::T
    angles::Vector{Float64} # in degrees

function DirectivitySpectra(spectra::T, angles::Vector{Float64}) where {E, T <: AbstractSpectrumArray{E}}
        @assert nchannels(spectra) == length(angles)
        @assert length(spectra) > 0
        s = deepcopy(spectra)
        a = normdeg.(Base.vect(angles...))
        @assert length(a) == length(unique(a))
        idxs = sortperm(a)
        new{T}(s[:, idxs], a[idxs])
    end
end

DirectivitySpectra(spectra::AbstractVector{T}, angles::AbstractVector{<:Real}) where {E, T <: AbstractSpectrumArray{E}} = 
    DirectivitySpectra(hcat(spectra...), Float64.(angles))

spectra(ds::DirectivitySpectra) = ds.spectra
angles(ds::DirectivitySpectra) = ds.angles
SampleArrays.data(ds::DirectivitySpectra) = data(ds.spectra)
SampleArrays.rate(ds::DirectivitySpectra) = rate(ds.spectra)
SampleArrays.data_no0(ds::DirectivitySpectra) = data_no0(ds.spectra)
SampleArrays.domain_no0(ds::DirectivitySpectra) = domain_no0(ds.spectra)
SampleArrays.nchannels(ds::DirectivitySpectra) = nchannels(ds.spectra)
SampleArrays.nfreqs(ds::DirectivitySpectra) = nfreqs(ds.spectra)

SampleArrays.domain(ds::DirectivitySpectra) = domain(ds.spectra)
SampleArrays.toindex(ds::DirectivitySpectra, t::Frequency) = toindex(ds.spectra, t)

function Base.similar(ds::DirectivitySpectra, t::Type{T}, dims::Dims) where {E, T <: AbstractSpectrumArray{E}}
    dims[2] != nchannels(ds) && throw(ArgumentError("changing number of angles/channels not supported!"))
    DirectivitySpectra(similar(ds.spectra, E, dims), angles(ds))
end

function Base.setindex!(ds::DirectivitySpectra{T}, v, i::Int) where {E, T <: AbstractSpectrumArray{E}}
    @show "setindex", v, i
    setindex!(data(ds), v, i)::Matrix{E}
end

const RangeIndex = Union{Colon, AbstractRange, Vector{Bool}, Vector{Int}}

# Base.getindex(ds::DirectivitySpectra{T}, i::Int) where {E, T <: AbstractSpectrumArray{E}} = data(ds)[i]::E

Base.getindex(ds::DirectivitySpectra, I::R) where {R <: RangeIndex} =
    DirectivitySpectra(ds.spectra[I, :], angles(ds))

Base.getindex(ds::DirectivitySpectra, I1::R, I2::S) where {R <: RangeIndex, S <: RangeIndex} =
    DirectivitySpectra(ds.spectra[I1, I2], angles(ds)[I2])

Base.getindex(ds::DirectivitySpectra, i::Integer) = ds.spectra[:, i:i]

function Base.getindex(ds::DirectivitySpectra, α::AbstractFloat)
    ch = closestangle(ds, α).idx
    ds[ch]
end

"""
    getindex(ds::DirectivitySpectra, α::AbstractFloat, β::AbstractFloat)

Take indices with angle ∈ [α, β] for α <= β.
"""
function Base.getindex(ds::DirectivitySpectra, α::AbstractFloat, β::AbstractFloat)
    @assert α <= β
    # ds[:, last(between(ds, α)):first(between(ds, β))]
    ds.spectra[:, last(between(ds, α)):first(between(ds, β))]
end

# ERASE?
# Base.getindex(ds::DirectivitySpectra, i::UnitRange{Int64}) = error("CHECK")
# Base.getindex(ds::DirectivitySpectra, i::UnitRange{Int64}) = DirectivitySpectra(ds.spectra[:, i], ds.angles[i])

# Base.getindex(ds::DirectivitySpectra, i::BitArray{1}) = error("CHECK")
# Base.getindex(ds::DirectivitySpectra, i::BitArray{1}) = DirectivitySpectra(ds.spectra[:, i], ds.angles[i])

# Base.getindex(ds::DirectivitySpectra, i::Interval) = DirectivitySpectra(ds.spectra[i, :], ds.angles)

maximumamp(ds::DirectivitySpectra) = maximum(abs.(ds.spectra))
minimumamp(ds::DirectivitySpectra) = minimum(abs.(ds.spectra)) 

# ERASE?
# function SampleArrays.zerophase!(ds::DirectivitySpectra) 
#     zerophase!(ds.spectra)
#     ds
# end

normalizeamp(ds::DirectivitySpectra) = normalizeamp!(deepcopy(ds))

function normalizeamp!(ds::DirectivitySpectra)
    maxamp = maximumamp(ds)
    maxamp == 0 && return ds
    ds.spectra ./= maxamp
    return ds
end

normalizeampangle(ds::DirectivitySpectra) = normalizeampfreq!(deepcopy(ds))

function normalizeampangle!(ds::DirectivitySpectra)
    for α in 1:length(ds)
        maxamp = maximum(abs.(ds.spectra[:, α]))
        maxamp == 0 && continue
        ds.spectra[:, α] ./= maxamp
    end
    return ds
end


normalizeampfreq(ds::DirectivitySpectra) = normalizeampfreq!(deepcopy(ds))

function normalizeampfreq!(ds::DirectivitySpectra)
    for fidx in 1:nfreqs(ds)
        maxamp = maximum(abs.(data(spectra(ds))[fidx, :]))
        maxamp == 0 && continue
        data(spectra(ds))[fidx, :] ./= maxamp
    end
    return ds
end

function between(ds::DirectivitySpectra, α)
    α = normdeg(α)
    nc = nchannels(ds)
    if nc == 1
        return 1, 1
    end
    as = [ds.angles..., ds.angles[1]]
    for i in 1:nc
        j = i % nc + 1
        d1, d2 = deltadeg(ds.angles[i], α), deltadeg(α, ds.angles[j])
        if d1 >= 0 && d2 >= 0
            if d1 == 0
                return i, i
            elseif d2 == 0
                return j, j
            end
            return i, j
        end
    end
end

# function Interpolations.interpolate(ds::DirectivitySpectra, α)
#     error("FIX")
#     i, j = between(ds, α)
#     if i == j
#         return deepcopy(ds.spectra[:, i])
#     end
#     β, γ = ds.angles[i], ds.angles[j]
#     b, g = deltadeg(β, α), deltadeg(α, γ)
#     @assert b >= 0 && g >= 0
#     s = b + g
#     @assert s > 0
#     return mean([2b/s, 2g/s]' .* ds.spectra[:, [i, j]])
# end

relativeto(ds::DirectivitySpectra, to) = relativeto!(deepcopy(ds), to)

function relativeto!(ds::DirectivitySpectra, to::AbstractSpectrumArray)
    @assert nchannels(to) == 1
    to = deepcopy(to.data)
    any(to .== 0.0) && @warn("relativeto!: zero in denominator!")
    ds.spectra.data ./= to
    return ds
end

relativeto!(ds::DirectivitySpectra, α::Real) = relativeto!(ds, ds[Float64(α)])

function closestangle(ds::DirectivitySpectra, α)
    closeness = abs.(ds.angles .- α) 
    _, idx = findmin(closeness) # index of angle spectra closest to the angle
    (idx=idx, α=ds.angles[idx])
end

"""
    mirror(ds::DirectivitySpectra)

Mirror all spectra, i.e., append copies of spectra with negative angles (if not already in the collection).
"""
function mirror(ds::DirectivitySpectra)
    newspectra = []
    newangles = deepcopy(ds.angles)
    for i in 1:nchannels(ds)
        α = normdeg(-newangles[i])
        if α ∉ newangles
            push!(newangles, α)
            push!(newspectra, ds[i])
        end
    end
    newspectra = hcat(ds.spectra, hcat(newspectra...))
    idxs = sortperm(newangles)
    newspectra .= newspectra[:, idxs]
    newangles .= newangles[idxs]
    return DirectivitySpectra(newspectra, newangles)
end

smooth(ds::DirectivitySpectra{T}, octave::Real = 1//3) where {T <: RFFTSpectrumArray} = 
    DirectivitySpectra(smooth(ds.spectra, octave), angles(ds))

mean(ds::DirectivitySpectra) = mean(ds.spectra)
