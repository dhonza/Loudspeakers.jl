export AngularMeasurements, AngularSpectra, parent, angles, angleidx, closest_angle, between_indices

# export maximumamp, minimumamp, normalizeamp, normalizeamp!
# export normalizeampangle, normalizeampangle!
# export normalizeampfreq, normalizeampfreq!
# export relativeto, relativeto!, smooth, interpolate, mirror, mirror!


struct AngularMeasurements{T,P<:AbstractSampleArray{T},A<:AbstractVector{<:Real}} <: AbstractSampleArray{T}
    parent::P
    angles::A
    
    function AngularMeasurements(parent, angles)
        nchannels(parent) == length(angles) || throw(ArgumentError("number of channels ($(nchannels(a))) != number of angles ($(length(angles)))"))
        s = deepcopy(parent) # because we sort it below
        a = normdeg.(Base.vect(angles...)) # normalize degrees to (-180°, 180°]
        length(a) == length(unique(a)) || throw(ArgumentError("non-unique angles!"))
        idxs = sortperm(a)
        new{eltype(parent), typeof(parent), typeof(angles)}(s[:, idxs], collect(a[idxs]))
    end
end

const AngularSpectra{T,P<:AbstractSpectrumArray{T},A<:AbstractVector{<:Real}} = AngularMeasurements{T,P,A}

parent(am::AngularMeasurements) = am.parent
angles(am::AngularMeasurements) = am.angles

SampleArrays.rate(am::AngularMeasurements) = rate(parent(am))
SampleArrays.data(am::AngularMeasurements) = data(parent(am))
SampleArrays.data_no0(as::AngularSpectra) = data_no0(parent(as))

SampleArrays.domain(am::AngularMeasurements) = domain(parent(am))
SampleArrays.domain_no0(as::AngularSpectra) = domain_no0(parent(as))

SampleArrays.names(am::AngularMeasurements) = names(parent(am))
SampleArrays.names!(am::AngularMeasurements) = names!(parent(am))
SampleArrays.nchannels(am::AngularMeasurements) = nchannels(parent(am))
SampleArrays.nframes(am::AngularMeasurements) = nframes(parent(am))

SampleArrays.ntimeframes(as::AngularSpectra) = ntimeframes(parent(as))

function Base.show(io::IO, t::MIME"text/plain", am::AngularMeasurements)
    as = angles(am)
    println(io, "AngularMeasurements: $(minimum(as))°..$(maximum(as))°:")
    show(io, t, parent(am))
end

function angleidx(angles::AbstractVector{<:Real}, α; func=(==))
    # find first angle index (and the actual angle), for which func is true
    for (i, a) in enumerate(angles)
        func(a, α) && return (idx=i, α=a)
    end
    error("angle ($(α)) not found!")
end

angleidx(am::AngularMeasurements, α) = angleidx(angles(am), α) 

function closest_angle(angles::AbstractVector{<:Real}, α)
    closeness = abs.(angles .- α) 
    _, idx = findmin(closeness) # index of angle spectra closest to the angle α
    (idx=idx, α=angles[idx])
end

closest_angle(am::AngularMeasurements, α) = closest_angle(angles(am), α) 

function between_indices(am::AngularMeasurements, α)
    α = normdeg(α)
    nc = nchannels(am)
    if nc == 1
        return 1, 1
    end
    as = [am.angles..., am.angles[1]]
    for i in 1:nc
        j = i % nc + 1
        d1, d2 = deltadeg(am.angles[i], α), deltadeg(α, am.angles[j])
        if d1 >= 0 && d2 >= 0
            if d1 ≈ 0
                return i, i
            elseif d2 ≈ 0
                return j, j
            end
            return i, j
        end
    end
end

function Base.getindex(am::AngularMeasurements, I::R, J::S) where {R <: SampleArrays.FrameIndex, S <: SampleArrays.ChannelIndex}
    idxs = SampleArrays.tochannelidx(parent(am), J)
    p = parent(am)[I, idxs]
    AngularMeasurements(p, angles(am)[idxs])
end
Base.getindex(am::AngularMeasurements, I::R) where {R <: SampleArrays.FrameIndex} = am[I, :]  

Base.getindex(am::AngularMeasurements, α::AbstractFloat) = am[:, angleidx(am, α).idx]
function Base.getindex(am::AngularMeasurements, α::AbstractFloat, β::AbstractFloat)
    α <= β || error("α > β!")
    am[:, last(between_indices(am, α)):first(between_indices(am, β))]
end

function Base.similar(am::AngularMeasurements, t::Type, dims::Dims)
    AngularMeasurements(similar(parent(am), t, dims), angles(am))
end

# maximumamp(ds::DirectivitySpectra) = maximum(abs.(ds.spectra))
# minimumamp(ds::DirectivitySpectra) = minimum(abs.(ds.spectra)) 

# normalizeamp(ds::DirectivitySpectra) = normalizeamp!(deepcopy(ds))

# function normalizeamp!(ds::DirectivitySpectra)
#     maxamp = maximumamp(ds)
#     maxamp == 0 && return ds
#     ds.spectra ./= maxamp
#     return ds
# end

# normalizeampangle(ds::DirectivitySpectra) = normalizeampfreq!(deepcopy(ds))

# function normalizeampangle!(ds::DirectivitySpectra)
#     for α in 1:length(ds)
#         maxamp = maximum(abs.(ds.spectra[:, α]))
#         maxamp == 0 && continue
#         ds.spectra[:, α] ./= maxamp
#     end
#     return ds
# end


# normalizeampfreq(ds::DirectivitySpectra) = normalizeampfreq!(deepcopy(ds))

# function normalizeampfreq!(ds::DirectivitySpectra)
#     for fidx in 1:nfreqs(ds)
#         maxamp = maximum(abs.(data(spectra(ds))[fidx, :]))
#         maxamp == 0 && continue
#         data(spectra(ds))[fidx, :] ./= maxamp
#     end
#     return ds
# end

function Interpolations.interpolate(am::AngularMeasurements, α)
    i, j = between_indices(am, α)
    if i == j
        return am.parent[:, i]
    end
    β, γ = angles(am)[[i, j]]
    b, g = Float32(deltadeg(β, α)), Float32(deltadeg(α, γ)) # TODO better way to keep Float32 of parent(am)!
    @assert b >= 0 && g >= 0
    s = b + g
    @assert s > 0
    return mean([2g/s, 2b/s]' .* parent(am)[:, [i, j]])
end

# relativeto(ds::DirectivitySpectra, to) = relativeto!(deepcopy(ds), to)

# function relativeto!(ds::DirectivitySpectra, to::AbstractSpectrumArray)
#     @assert nchannels(to) == 1
#     to = deepcopy(to.data)
#     any(to .== 0.0) && @warn("relativeto!: zero in denominator!")
#     ds.spectra.data ./= to
#     return ds
# end

# relativeto!(ds::DirectivitySpectra, α::Real) = relativeto!(ds, ds[Float64(α)])

# """
#     mirror(ds::DirectivitySpectra)

# Mirror all spectra, i.e., append copies of spectra with negative angles (if not already in the collection).
# """
# function mirror(ds::DirectivitySpectra)
#     newspectra = []
#     newangles = deepcopy(ds.angles)
#     for i in 1:nchannels(ds)
#         α = normdeg(-newangles[i])
#         if α ∉ newangles
#             push!(newangles, α)
#             push!(newspectra, ds[i])
#         end
#     end
#     newspectra = hcat(ds.spectra, hcat(newspectra...))
#     idxs = sortperm(newangles)
#     newspectra .= newspectra[:, idxs]
#     newangles .= newangles[idxs]
#     return DirectivitySpectra(newspectra, newangles)
# end

# smooth(ds::DirectivitySpectra{T}, octave::Real = 1//3) where {T <: RFFTSpectrumArray} = 
#     DirectivitySpectra(smooth(ds.spectra, octave), angles(ds))

# mean(ds::DirectivitySpectra) = mean(ds.spectra)
