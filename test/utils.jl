function compare_magdb(a, b, id, minf, maxf)
    # expects common domains
    @assert nchannels(a) == 1
    @assert nchannels(b) == 1
    @assert domain(a) == domain(b)
    dommask = map(f ->  minf <= f <= maxf, domain(a))
    a, b = SpectrumArray.((a, b))
    
    a_ = data(amp2db.(abs.(a[dommask])))
    b_ = data(amp2db.(abs.(b[dommask])))
    errs = abs.(a_ .- b_)
    errmask = (!).(isinf.(errs) .| isnan.(errs))

    mae = mean(errs[errmask])
    maxerr = maximum(errs[errmask])
    errs_ = SpectrumArray(errs .|> db2amp, rate(a), domain(a)[dommask], [:error])
    (id = id, mae = mae, maxerr = maxerr, errs = errs_)
end

function compare_phasedeg(a, b, id, minf, maxf)
    # TODO: merge with vituix_spinorama
    # expects common domains
    @assert nchannels(a) == 1
    @assert nchannels(b) == 1
    @assert domain(a) == domain(b)
    dommask = map(f ->  minf <= f <= maxf, domain(a))
    a, b = SpectrumArray.((a, b))

    # normalization and conversion complex -> coordinates
    c2v(c) = (nc = c / abs(c); [real(nc), imag(nc)])
    a_ = data(c2v.(a[dommask]))
    b_ = data(c2v.(b[dommask]))
    errs = acosd.(clamp.(dot.(a_, b_), -1.0, 1.0)) # clamp for numerical stability

    errmask = (!).(isinf.(errs) .| isnan.(errs))
    mae = mean(errs[errmask])
    maxerr = maximum(errs[errmask])
    errs_ = SpectrumArray(errs .|> db2amp, rate(a), domain(a)[dommask], [:error])
    (id = id, mae = mae, maxerr = maxerr, errs = errs_)
end

function test_driver_transform(driver_dir, test_dir, freqrange, comparef,
    threshold_max, threshold_mae, matchFIR; 
    template_hor=Regex(raw"VituixCAD_PolarFR hor (.*)\.frd"), template_ver=Regex(raw"VituixCAD_PolarFR ver (.*)\.frd"), 
    transform...)

    driver_spinorama = load_spinorama_from_FRDs(driver_dir, template_hor=template_hor, template_ver=template_ver)
    target_spinorama = load_spinorama_from_FRDs(test_dir, template_hor=template_hor, template_ver=template_ver)

    drivers = [Dict(:spinorama => driver_spinorama, :transform => Dict{Symbol,Any}(transform...))]
    test_transform(drivers, target_spinorama, 
        freqrange, comparef, threshold_max, threshold_mae, matchFIR)
end

function test_transform(drivers, target_spinorama, freqrange, comparef,
    threshold_max, threshold_mae,
    matchFIR) # true -> interpolate both spinorama to the same type and size as H, false -> covert H to match the target_spinorama
    
    minf, maxf = toHz(first(freqrange)), toHz(last(freqrange))

    H1 = get(drivers[1][:transform], :H, nothing)
    if !isnothing(H1)
        if matchFIR
            for driver in drivers
                driver[:spinorama] = interpolate(driver[:spinorama], H1)
            end
            target_spinorama = interpolate(target_spinorama, H1)
        else
            for driver in drivers
                driver[:transform][:H] = interpolate(driver[:transform][:H], parent(target_spinorama.hor))
            end
        end
    end

    errs = []
    for ref in angles(target_spinorama.hor)
        X1 = project_vituix(drivers...; lat=0.0, lon=ref, ldist=2.0)
        X2 = parent(target_spinorama.hor[ref])
        push!(errs, comparef(X1, X2, (:hor, ref), minf, maxf))
    end

    # I found to have considerable discrepancies in my measurements of hor-180 and ver-180, horizontal measurement 
    # is the default one, so remove 180° from vertical
    for ref in filter(a -> a != 180, angles(target_spinorama.ver))
        lat, lon = normlat(ref)
        X1 = project_vituix(drivers...; lat=lat, lon=lon, ldist=2.0)  
        X2 = parent(target_spinorama.ver[ref])
        push!(errs, comparef(X1, X2, (:ver, ref, (lat=lat, lon=lon)), minf, maxf))
    end

    sort!(errs; by=e -> e.maxerr, rev=true)

    worst_maxerr = errs[1].maxerr
    worst_mae = errs[1].mae
    println("\tMAX: $(worst_maxerr), \tMAE: $(worst_mae)\tWORST: ", errs[1].id)
    if worst_maxerr > threshold_max || worst_mae > threshold_mae
        return false
    end
    return true
end

function test_interpolation(input_dir, freqrange, comparef, threshold_max, threshold_mae)
    target_spinorama = load_spinorama_from_FRDs(input_dir)
    
    minf, maxf = toHz(first(freqrange)), toHz(last(freqrange))

    cond(a) = abs(a % 10) == 0 
    select_angles(am, cond) = am[:, map(cond, angles(am))]
    test_spinorama = Spinorama(select_angles(target_spinorama.hor, cond), select_angles(target_spinorama.ver, cond))
    
    errs = []
    
    for ref in angles(select_angles(target_spinorama.hor, c -> !cond(c)))
        X1 = project(test_spinorama, 0.0, ref)
        X2 = parent(target_spinorama.hor[ref])
        push!(errs, comparef(X1, X2, (:hor, ref), minf, maxf))
    end

    for ref in angles(select_angles(target_spinorama.ver, c -> !cond(c)))
        X1 = project(test_spinorama, normlat(ref)...) 
        X2 = parent(target_spinorama.ver[ref])
        push!(errs, comparef(X1, X2, (:ver, ref), minf, maxf))
    end
    
    sort!(errs; by=e -> e.mae, rev=true)

    worst_maxerr = errs[1].maxerr
    worst_mae = errs[1].mae
    println("\tMAX: $(worst_maxerr), \tMAE: $(worst_mae)\tWORST: ", NamedTuple{(:id, :maxerr, :mae)}(errs[1]))
    if worst_maxerr > threshold_max || worst_mae > threshold_mae
        return false
    end
    return true
end