export record_response_ecasound, measure_response, measurement_info, calibrate_response, compute_ir
export table_measurement
export postprocess_response, postprocess_dir

function record_response_ecasound(fstimulus, fresponse; duration=2s)
    durationsec = ustrip(Float32, u"s", duration)
    cmd = `ecasound -t:$durationsec -f:16,2,44100 -a:1 -i $fstimulus -o alsa,pulse -a:2 -i alsa,pulse -o $fresponse`
    run(pipeline(cmd, devnull))
    readwav(fresponse)
end

function measure_response(stimulus, fstimulus, fresponse; postduration=1s)  
    #TODO fix API so stimulus is not always converted....
    xstimulus = SampleArray(stimulus.sig, stimulus.samplerate)
    duration = getduration(xstimulus)
    xresponse = record_response_ecasound(fstimulus, fresponse; duration=(duration + postduration))
    xresponse
end

function measurement_info(xresponse; c=344u"m/s")
    duration = uconvert(u"s", getduration(xresponse))
    tdiff_samples = DSP.finddelay(xresponse[:, 2:2], xresponse[:, 1:1]);
    tdiff = uconvert(u"ms", tdiff_samples/rateHz(xresponse))
    tdiff_distance = uconvert(u"mm", c * tdiff)
    return (duration=duration, tdiff=tdiff, tdiff_samples=tdiff_samples, tdiff_distance)
end

function calibrate_response(xresponse, calfile="../data/34Q956_cal_0degree.txt")
    # TODO make this more general, e.g., sound card calibration
    Xresponse = rfft(xresponse)
    Xresponse[:, 2:2] .= miccalibrate(Xresponse[:, 2:2], calfile)
    xresponse = irfft(Xresponse)
    xresponse, Xresponse
end

function compute_ir(stimulus, xresponse; noncausal=true)
    # find loopback delay caused by the processing
    # shift signals, keeping only the mic channel
    xstimulus = SampleArray(stimulus.sig, stimulus.samplerate)
    tdiff = DSP.finddelay(xstimulus, xresponse[:, 1:1])
    xresponse = DSP.shiftsignal(xresponse, tdiff)[:, 2:2]
    
    h = analyze(stimulus, xresponse, noncausal=noncausal)
    H = rfft(h)
    h, H
end

function table_measurement(stimulus, fstimulus, rawdir, tableAPI, degs=0:5:355; 
        response_template="response_{deg}.wav",
        rate=44100Hz, duration=5s, postduration=1s, speed=0.1, dampingsleep=2s, c=344m/s)
    tableAPI("speed/$(speed)")
    mkpath(rawdir)
    iter = ProgressBar(degs)
    for deg in iter
        set_description(iter, "$(deg)°")
        tableAPI("rotateto/$(deg)")
        sleep(ustrip(Float32, u"s", dampingsleep)) # wait for damping vibrations
        fresponse = joinpath(rawdir, replace(response_template, "{deg}" => @sprintf("%03d", deg)))
        isfile(fresponse) && error("response file already exists: $(fresponse)!")
        
        xresponse = measure_response(stimulus, fstimulus, fresponse; postduration=postduration)
        minfo = measurement_info(xresponse; c=c)
        set_postfix(iter, Lag="$(minfo.tdiff_samples) spl")
        # @info "recorded $(nframes(xresponse)) samples ($(minfo.duration)); mic lags the loopback by $(minfo.tdiff_samples) samples ≈ $(minfo.tdiff) ≈ $(minfo.tdiff_distance)" 
    end
end

function postprocess_response(infile, outdir, stimulus; 
    response_template=r"response_(\d\d\d).wav",
    shift=0, # shift due to imprecise mic positioning
    calfile=nothing,
    ir_template="ir_{deg}.wav", noncausal=true)

    m = match(response_template, basename(infile))
    if !isnothing(m)
        mkpath(outdir)
        
        deg = parse(Int, m.captures[1])
        if deg < 0
            deg = 360 - deg
        end
        xresponse = readwav(infile)
        @assert nchannels(xresponse) == 2

        # shift the second channel which is the mic recording
        xresponse[:, 2:2] .= DSP.shiftsignal(xresponse[:, 2:2], shift)

        if !isnothing(calfile)
            xresponse, _ = calibrate_response(xresponse, calfile)
        end
        fimpulse = joinpath(outdir, replace(ir_template, "{deg}" => @sprintf("%03d", deg)))
        isfile(fimpulse) && error("impulse file already exists: $(fimpulse)!")

        h, _ = compute_ir(stimulus, xresponse; noncausal=noncausal)
        writewav(h, fimpulse; nbits=16, compression=WAVE_FORMAT_PCM)
        return true
    end
    return false
end

function postprocess_dir(indir, func)
    res = []
    iter = ProgressBar(readdir(indir))
    for f in iter
        set_description(iter, f)
        infile = joinpath(indir, f)
        if isfile(infile)
            push!(res, func(infile))
        end
    end
    @info "processed $(sum(res)) files"
end