export irmid, irfirstbig

# ---- IR PEAK DETECTION ---------------------------------------------

irmid(x) = round(Int, length(x) / 2, RoundNearestTiesAway) # the first peak is in the middle

function irfirstbig(x::AbstractVector, threshold=0.5) # the first peak is the first one having size greater than threshold (normalized absoulute IR) 
    x2 = abs.(x)
    x2 ./= maximum(x2)
    for i in 1:length(x2)
        if x2[i] > threshold
            return i
        end
    end
    @warn "irfirstbig: not found, returning first index!"
    return 1
end

function irfirstbig(a::SampleArray, threshold=0.5)
    nchannels(a) != 1 && throw(ArgumentError("number of channels != 1 ($(nchannels(a)))"))
    irfirstbig(data(a)[:, 1], threshold)
end

# Base.getindex(x::SampleBuf{T, 1} where T, pidx::Int, pre::Int, post::Int) = x[pidx-pre:pidx+post]
# Base.getindex(x::SampleBuf{T, 1} where T, pidx::Int, pre, post) = x[pidx, _convtime(x, pre), _convtime(x, post)]
# Base.getindex(x::SampleBuf{T, 1} where T, peakf::Function, pre, post) = x[peakf(x), _convtime(x, pre), _convtime(x, post)]