export audiofreqlogticks, audiofreqlogtickformatter

function audiofreqlogticks()
    [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 
        11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000]
end

function audiofreqlogtickformatter(x; ticks=[1, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000])
    v = round(Int, x)
    if v ∈ ticks
        if v >= 1000
            v = round(v/1000; digits=1)
            v = v == round(v) ? round(Int, v) : v
            return string(v) * "k"
        end
        return v
    end
    return ""
end

@recipe function f(X::AbstractSpectrumArray, ytype = :magdb; smoothing = nothing)
    xscale --> :log10
    xlims --> (1, Inf)
    xguide --> "Hz"
    size --> (1000, 500)
    xformatter --> audiofreqlogtickformatter
    grid --> true
    xticks --> audiofreqlogticks()
    label --> reshape(replace.(String.(names(X)), "_" => " "), 1, :)

    if !isnothing(smoothing)
        X = smooth(X, smoothing)
    end

    if get(plotattributes, :xscale, nothing) == :log10
        arr = data_no0(X)
        dom = domain_no0(X)
    else
        arr = data(X)
        dom = domain(X)
    end

    if ytype == :raw
        b = arr
    elseif ytype == :magdb
        yguide --> "magnitude (dB)"
        b = amp2db.(abs.(arr))
    elseif ytype == :mag
        yguide --> "magnitude"
        b = abs.(arr)
    elseif ytype == :powdb
        yguide --> "power (dB)"
        b = pow2db.(abs.(arr))
    elseif ytype == :phase
        yguide --> "deg"
        b = rad2deg.(angle.(arr))
    elseif ytype == :uphase # unwrapped
        yguide --> "deg (unwrapped)"
        b = rad2deg.(unwrap(angle.(arr); dims=1))
    elseif ytype == :re
        yguide --> "Re"
        b = real.(arr)
    elseif ytype == :im
        yguide --> "Im"
        b = imag.(arr)
    else
        throw(ArgumentError("unknown ytype = $ytype, use one of :raw, :mag, :magdb, :powdb, :phase, :uphase, :re, :im"))
    end
    dom, b
end
