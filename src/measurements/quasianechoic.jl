export quasi_anechoic, quasi_anechoic_window, boundary_delay

function quasi_anechoic_window(a::SampleArray, fstart=irfirstbig; post, pre=125ms)
    @assert nchannels(a) == 1
    imax = fstart(a)
    ipre = toindexdelta(a, pre) # index deltas
    ipost = toindexdelta(a, post)
    n = ipost + ipre + 1
    # @show imax, n, ipre, ipost
    # TODO make the two-sided Tukey window properly
    win = zeros(Float64, n)    
    win[1:ipre+1] .= tukey(ipre+1, 0.25)
    win[end-ipost:end] .= tukey(ipost+1, 0.25)
    win[ipre÷2:end-ipost÷2] .= 1.0
    fullwin = similar(a)
    fullwin[imax-ipre:imax+ipost, 1:1] .= win
    fullwin[1:imax-ipre-1] .= 0
    fullwin[imax+ipost+1:end] .= 0

    (win=fullwin, imax=imax, ipre=ipre, ipost=ipost)
end

function quasi_anechoic(a::SampleArray, fstart=irfirstbig; post, pre=125ms, cut=true)
    win, imax, ipre, ipost = quasi_anechoic_window(a, fstart; post=post, pre=pre)
    q = a .* win
    if cut
        return q[imax-ipre:imax+ipost, :]
    else
        return q
    end
end

# TODO make c a general constant
function boundary_delay(micheight, micdist; c=342.0)
    refdist = 2sqrt(micheight^2 + (micdist/2)^2) # distance of the reflection
    deltadist = refdist - micdist
    directtime = micdist/c
    reftime = refdist/c
    delay = deltadist/c
    (minfreq=1/delay, delay=delay, directtime=directtime, reftime=reftime, deltadist=deltadist, refdist=refdist)
end
