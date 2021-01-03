export smooth
using FFTViews

function smooth(a::RFFTSpectrumArray{T}, octave::Real = 1//3) where {T <: Real}
    # algorithm: Hatziantoniou, Panagiotis D., and John N. Mourjopoulos. "Generalized fractional-octave smoothing of audio and acoustic responses." Journal of the Audio Engineering Society 48.4 (2000): 259-280.
    # but based on description in:
    # Tylka, Joseph G., Braxton B. Boren, and Edgar Y. Choueiri. "A generalized method for fractional-octave smoothing of transfer functions that preserves log-frequency symmetry." Journal of the Audio Engineering Society 65.3 (2017): 239-245.
    # Method 1: Symmetric weights
    # rectangular window
    octave ≥ 0 || error("octave must be ≥ 0")
    n = nframes(a)
    as = similar(a)
    X = data(a)
    Xs = data(as)

    X = [X; conj.(X[end-1:-1:2,:])] # we expect spectra created using rfft(), so extend by negative frequencies

    rQ = 2sinh(octave*log(2)/2) # Q reciprocial: 1/Q CHECK use log2?
    for c in 1:nchannels(a)
        XV = FFTView(view(X, :, c)) # make it circular and starting on zero index
        plb, pub = 0, -1 # previous lower and upper bounds
        msum = 0.0 # moving sum
        for k in 0:n-1
            mk = floor(Int, k*rQ/2)
            lb, ub = k-mk, k+mk
            msum = msum - sum(XV[plb:lb-1]) + sum(XV[pub+1:ub])
            mmean = msum / (2mk + 1)
            Xs[k+1, c] = mmean
            plb, pub = lb, ub
        end
    end
    as
end

function smooth_naive(a::RFFTSpectrumArray{T}, octave::Real = 1//3) where {T <: Real}
    # slower version of smooth - always recomputes window means
    octave ≥ 0 || error("octave must be ≥ 0")
    n = nframes(a)
    as = similar(a)
    X = data(a)
    Xs = data(as)

    X = [X; conj.(X[end-1:-1:2,:])] # we expect spectra created using rfft(), so extend by negative frequencies
    X = FFTView(X) # make it circular and starting on zero index

    rQ = 2sinh(octave*log(2)/2) # Q reciprocial: 1/Q CHECK use log2?
    for c in 1:nchannels(a)
        for k in 0:n-1
            mk = floor(Int, k*rQ/2)
            # Xs[k+1, c] = mean(X[k-mk:k+mk, c])
            Xs[k+1, c] = mean(X[max(k-mk, 0):min(k+mk, n-1), c])
        end
    end
    as
end

function smooth(X::RFFTSpectrumArray{T}, octave::Real = 1//3) where {E, T <: SignalElement{E}}
    Y = similar(X)
    # mags = smooth_naive(abs.(X), octave)
    mags = smooth(abs.(X), octave)
    phis = smooth(unwrap(angle.(X), dims=1), octave)
    Y .= (MagPhase(max(m, 0.0), p) for (m, p) in zip(mags, phis)) # max for numerical instability giving negative magnitude
    Y
end
