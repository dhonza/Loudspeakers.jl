@testset "Spinorama" begin
    function compare_magdb(a, b, c)
        # expects common domains
        errs = abs.(data(amp2db.(a[:, c]) .- amp2db.(b[:, c])))
        mae = mean(errs)
        maxerr = maximum(errs)
        errs_ = SpectrumArray(errs .|> db2amp, rate(a), domain(a), [:error])
        (id=c, mae=mae, maxerr=maxerr, errs=errs_)
    end

    function compare_curves(a, b, cmpfunc)
        common_channels = collect(intersect(Set(names(a)), Set(names(b)))) |> sort
        res = [cmpfunc(unify_domains(a, b)..., c) for c in common_channels]
        sort!(res; by=e -> e.maxerr, rev=true)
        res
    end

    function test_spinorama(measurement_dir, threshold_max, threshold_mae; 
            template_hor=Regex(raw"VituixCAD_PolarFR hor (.*)\.frd"), template_ver=Regex(raw"VituixCAD_PolarFR ver (.*)\.frd"))
        data_dir = joinpath(test_dir, measurement_dir)
        spinorama = load_spinorama_from_FRDs(data_dir, template_hor=template_hor, template_ver=template_ver)
        spinorama_curves = evaluate(spinorama)
        spinorama_curves_vituix = load_vituix_spinorama_curves(joinpath(data_dir, "spinorama"))
        errs = compare_curves(spinorama_curves, spinorama_curves_vituix, compare_magdb)

        worst_maxerr = errs[1].maxerr
        worst_mae = errs[1].mae
        println("  \"$(measurement_dir)\"")
        println("\tMAX: $(worst_maxerr),\tMAE: $(worst_mae) WORST CASE: ", NamedTuple{(:id, :maxerr, :mae)}(errs[1]))
        if worst_maxerr > threshold_max || worst_mae > threshold_mae
            return false
        end
        return true
    end
    
    println("VITUIX SPINORAMA TESTS")
    @test test_spinorama("T/polar", 0.2, 0.08)
    @test test_spinorama("T/polar_noXO", 0.2, 0.08)
    @test test_spinorama("T/polar_r20", 0.3, 0.08)
    @test test_spinorama("T/polar_noXO_r-80_t-80", 0.15, 0.035)
    @test test_spinorama("T/polar_noXO_r20", 0.25, 0.075)
    @test test_spinorama("T/polar_r10_t15", 0.15, 0.07)
    
    @test test_spinorama("M/polar", 0.15, 0.05)
    @test test_spinorama("M/polar_noXO", 0.15, 0.05)
    @test test_spinorama("M/polar_r10", 0.15, 0.05)
    @test test_spinorama("M/polar_r170", 0.45, 0.075)
    @test test_spinorama("M/polar_y-118", 0.15, 0.035)
    @test test_spinorama("M/polar_y-118_r170", 0.4, 0.075)
    @test test_spinorama("M/polar_y-118_r180", 1.0, 0.15)
    @test test_spinorama("M/polar_noXO_r-170", 0.5, 0.075)
    @test test_spinorama("M/polar_noXO_r170", 0.45, 0.075)
    @test test_spinorama("M/polar_noXO_y-118_r170", 0.4, 0.065)
    @test test_spinorama("M/polar_x95_y-118_z150_r12_t-15", 0.15, 0.05)
    
    @test test_spinorama("W/polar", 0.8, 0.055)
    @test test_spinorama("W/polar_noXO", 0.8, 0.055)
    @test test_spinorama("W/polar_y-348", 0.6, 0.055)
    @test test_spinorama("W/polar_y-348_r14", 0.5, 0.045)
    @test test_spinorama("W/polar_noXO_y-348_r14", 0.5, 0.045)
    @test test_spinorama("W/polar_x-150_y-348_z-35_r-90_t3", 0.15, 0.08)
    
    @test test_spinorama("total/polar", 0.15, 0.045)
    @test test_spinorama("total_rotated/polar", 0.15, 0.04)
    @test test_spinorama("total_tilted/polar", 0.15, 0.045)
    @test test_spinorama("total_rotated_tilted/polar", 0.15, 0.04)
    
    @test test_spinorama("total_noXO/polar", 0.15, 0.05)
    @test test_spinorama("total_noXO_rotated/polar", 0.4, 0.06)
    @test test_spinorama("total_noXO_tilted/polar", 0.15, 0.05)
    @test test_spinorama("total_noXO_rotated_tilted/polar", 0.4, 0.05)

end