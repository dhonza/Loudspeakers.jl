export expsweep

# BASED ON: https://github.com/ssfrr/MeasureIR.jl/blob/master/src/expsweep.jl

"""
    expsweep(length; options...)
Create an impulse response measurement using an exponential sinusoid sweep from
`minfreq` to `maxfreq`, following the procedure from [1]. Durations and
frequencies can be given with or without Unitful units. With the default
`samplerate` of `1` durations can be interpreted as samples and frequencies are
in cycles/sample. The default minimum frequency of 0.0004 cycles/sample
corresponds to about 19Hz at 48kHz.
## Options
- `samplerate=1`
- `minfreq=0.0004`
- `maxfreq=0.5`
- `fadein`
- `fadeout`
- `fade`
- `optimize`
- `func`
## Examples
```julia
using MeasureIR
# generate a 10000-sample sweep
expsweep(10000)
# generate a 10s sweep at 48kHz samplerate
using Unitful: s, kHz
expsweep(10s; samplerate=48kHz)
```
[1]: Farina, Angelo, Simultaneous measurement of impulse response and distortion
with a swept-sine technique, Audio Engineering Society Convention 108 (2000).
"""
function expsweep(length;
        # TODO: consider making 2π the default samplerate, which makes the
        # frequencies interpretable as radians/sample
        samplerate = 1,
        minfreq=0.0004*samplerate,
        # TODO: this probably isn't the right default. for example for
        # meas = expsweep(30s; samplerate=48000Hz), it doesn't roll off on top
        maxfreq=0.5*samplerate,
        # default fades are computed below
        fadein=nothing,
        fadeout=nothing,
        fade=nothing,
        gain=1,
        optimize=true,
        func=sin)

    # TODO: review Novak: https://ant-novak.com/posts/research/2015-10-30_JAES_Swept/
    # it was recommended by Gordon

    w1 = uconvert(NoUnits, minfreq * 2π / samplerate)
    w2 = uconvert(NoUnits, maxfreq * 2π / samplerate)
    L = round(Int, uconvert(NoUnits, length*samplerate))

    if fade !== nothing
        if fadein !== nothing || fadeout !== nothing
            throw(ArgumentError("`fade` argument cannot be used with `fadein` or `fadeout`"))
        end
        fadein = fade
        fadeout = fade
    end

    # default fade times chosen heuristically to be pretty reasonable
    # we could definitely be smarter here, especially for small L
    if fadein === nothing
        f1 = round(Int, min(L/2, 200/w1))
    else
        f1 = round(Int, uconvert(NoUnits, fadein*samplerate))
    end

    if fadeout === nothing
        f2 = round(Int, min(L/4, 800/w2))
    else
        f2 = round(Int, uconvert(NoUnits, fadeout*samplerate))
    end

    if optimize
        w1 = _optimizew1(w1, w2, L)
    end
    sig = _expsweep(func, L, w1, w2)
    winin = 0.5 .- 0.5*cos.(range(0, π; length=f1))
    winout = 0.5 .- 0.5*cos.(range(π, 0; length=f2))
    sig[1:f1] .*= winin
    sig[(end-f2+1):end] .*= winout
    (sig=sig, w1=w1, w2=w2, gain=gain, samplerate=samplerate)
end

function _expsweep(func, L, minfreq, maxfreq)
    # switch to notation from [1]
    w1 = minfreq
    w2 = maxfreq
    T = L
    t = 0:(T-1)
    @. func(w1*T /
           log(w2/w1) * (exp(t/T * log(w2/w1))-1))
end

"""
Compute a new starting freq w < w1 that gives us an integer number
of total cycles.
"""
function _optimizew1(w1, w2, L)
    # scale to [0-1]
    f1 = w1 * L
    f2 = w2 * L
    # compute the total (integer) number of cycles we want
    M = floor((f2 - f1) /
              (2π*log(f2/f1)))

    a(f) = (f2-f)/log(f2/f) - 2π*M
    fnew = find_zero(a, f1)

    fnew/L
end
