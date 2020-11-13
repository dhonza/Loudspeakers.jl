@userplot DirectivityMap

@recipe function f(h::DirectivityMap)
    if length(h.args) != 1 || !(typeof(h.args[1]) <: DirectivitySpectra)
        error("Directivity MAp should be given single vector. Got: $(typeof(h.args))")
    end
    ds = h.args[1]

    size --> (1000, 500)
    xlims --> (20, 20000)
    ylims --> (-180,180)
    
    xscale --> :log10
    xticks --> audiofreqlogticks()
    xformatter --> audiofreqlogtickformatter
    yticks --> [-180, -135, -90, -45, 0, 45, 90, 135, 180]

    x = domain(ds)
    y = copy(ds.angles)
    Z = amp2db.(abs.(spectra(ds)))

    @series begin
        seriestype := :contour
        fill := true
        x, y, Z'
    end

    @series begin
        seriestype := :hline
        color := :black
        alpha := 0.3
        label := ""
        [-90, 0, 90]
    end
end

@userplot DirectivityPlot

@recipe function f(h::DirectivityPlot; db=true, minclip=0)
    if length(h.args) != 2 || !(typeof(h.args[1]) <: DirectivitySpectra) || !(typeof(h.args[2]) <: AbstractVector)
        error("Directivity Plot should be given two vectors. Got: $(typeof(h.args))")
    end
    ds, freqs = h.args

    if :ylims ∈ keys(plotattributes)
        ylims := plotattributes[:ylims] .- minclip
    end

    yformatter := y -> begin
        yt = y + minclip
        yt == round(yt) ? round(Int64, yt) : yt
    end

    amplitudes = abs.(spectra(ds))
    anglesrad = deg2rad.(angles(ds))
    for f in freqs
        freqamps = reshape(slice(amplitudes, f), :)
        if db
            freqamps = amp2db.(freqamps)
        end
        freqamps = max.(freqamps, minclip)
        freqamps .-= minclip 
        @assert all(freqamps .>= 0)
        push!(anglesrad, anglesrad[1])
        push!(freqamps, freqamps[1])
        @series begin
            proj := :polar
            label := "$f"
            anglesrad, freqamps
        end
    end
end

@userplot MeanDirectivityPlot

@recipe function f(h::MeanDirectivityPlot; db=true, minclip=0)
    if length(h.args) != 2 || !(typeof(h.args[1]) <: DirectivitySpectra) || !(typeof(h.args[2]) <: AbstractVector)
        error("Directivity Plot should be given two vectors. Got: $(typeof(h.args))")
    end
    ds, angles = h.args

    # mean non-weighted over all angles, and various front angle ranges
    @series begin
        label :="all"
        mean(ds)
    end
    for α in Float64.(angles)
        @series begin
            label :="≤ $(α)"
            mean(ds[-α, α])
        end
    end
end