@testset "Multiple Drivers" begin
    println("VITUIX MULTIPLE DRIVER TESTS")

    Hs = load_vituix_FIRs((joinpath(test_dir, "FIR", wav) for wav in ["T.wav", "M.wav", "W.wav"])...; Δt=-371.519ms, names=[:FIR_T, :FIR_M, :FIR_W])

    template_hor = Regex(raw"VituixCAD_PolarFR hor (.*)\.frd")
    template_ver = Regex(raw"VituixCAD_PolarFR ver (.*)\.frd")
    
    matchFIR = true

    test_set = [
            (src=[
                Dict(:name => "T", :spinorama => "T/polar_noXO", :transform => Dict(:H => Hs[:, :FIR_T], :y => 0.0)),
                Dict(:name => "M", :spinorama => "M/polar_noXO", :transform => Dict(:H => Hs[:, :FIR_M], :y => -0.118)),
                Dict(:name => "W", :spinorama => "W/polar_noXO", :transform => Dict(:H => Hs[:, :FIR_W], :y => -0.348)),
                ], target="total/polar", freqrange=(0Hz, 20kHz), T=Dict(), max_mag=14, mae_mag=0.2, max_phi=180.0, mae_phi=2.5),
            (src=[
                Dict(:name => "T", :spinorama => "T/polar_noXO", :transform => Dict(:H => Hs[:, :FIR_T], :r => 10.0, :t => 15.0)),
                Dict(:name => "M", :spinorama => "M/polar_noXO", :transform => Dict(:H => Hs[:, :FIR_M], :x => 0.095, :y => -0.118, :z => 0.150, :r => 12.0, :t => -15.0)),
                Dict(:name => "W", :spinorama => "W/polar_noXO", :transform => Dict(:H => Hs[:, :FIR_W], :x => -0.15, :y => -0.348, :z => -0.035, :r => -90.0, :t => 3.0)),
                ], target="total_rotated_tilted2/polar", freqrange=(0Hz, 20kHz), T=Dict(), max_mag=21.0, mae_mag=0.7, max_phi=180.0, mae_phi=10.5),
        ]
    
    for matchFIR in [true, false]
        println("match FIR: $matchFIR")
        for e in test_set
            sources = String[]
            drivers = []
            for s in e.src
                push!(sources, s[:spinorama])
                push!(drivers, Dict(
                    :spinorama => load_spinorama_from_FRDs(joinpath(test_dir, s[:spinorama]), template_hor=template_hor, template_ver=template_ver),
                    :transform => s[:transform]
                    ))
            end
            println("  \"$(sources)\" -> \"$(e.target)\"")

            target_spinorama = load_spinorama_from_FRDs(joinpath(test_dir, e.target), template_hor=template_hor, template_ver=template_ver)

            print("\tA:")
            @test test_transform(drivers, target_spinorama, e.freqrange, compare_magdb, e.max_mag, e.mae_mag, matchFIR)
            print("\tΦ:")
            @test test_transform(drivers, target_spinorama, e.freqrange, compare_phasedeg, e.max_phi, e.mae_phi, matchFIR)
        end
    end
end