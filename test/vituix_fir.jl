@testset "FIR" begin
    println("VITUIX FIR TESTS")

    Hs = load_vituix_FIRs((joinpath(test_dir, "FIR", wav) for wav in ["T.wav", "M.wav", "W.wav"])...; Δt=-371.519ms, names=[:FIR_T, :FIR_M, :FIR_W])

    test_set = [
        (src="T/polar_noXO", target="T/polar", freqrange=(1kHz, 20kHz), T=Dict(:H => Hs[:, :FIR_T]), max_mag=0.6, mae_mag=0.15, max_phi=180.0, mae_phi=2.5),
        (src="T/polar_noXO_r20", target="T/polar_r20", freqrange=(1kHz, 20kHz), T=Dict(:H => Hs[:, :FIR_T]), max_mag=0.6, mae_mag=0.15, max_phi=180.0, mae_phi=3.7),
        (src="T/polar_noXO", target="T/polar_r10_t15", freqrange=(1kHz, 20kHz), T=Dict([:H => Hs[:, :FIR_T], :r => 10.0, :t => 15.0]), max_mag=4.7, mae_mag=0.5, max_phi=180.0, mae_phi=11.0),

        (src="M/polar_noXO", target="M/polar", freqrange=(0kHz, 20kHz), T=Dict(:H => Hs[:, :FIR_M]), max_mag=1.9, mae_mag=0.06, max_phi=2, mae_phi=1.0),
        (src="M/polar_noXO_r170", target="M/polar_r170", freqrange=(0kHz, 20kHz), T=Dict(:H => Hs[:, :FIR_M]), max_mag=1.9, mae_mag=0.06, max_phi=180, mae_phi=2.2),
        (src="M/polar_noXO_y-118_r170", target="M/polar_y-118_r170", freqrange=(0kHz, 20kHz), T=Dict(:H => Hs[:, :FIR_M]), max_mag=1.9, mae_mag=0.06, max_phi=2, mae_phi=1.0),
        (src="M/polar_noXO", target="M/polar_x95_y-118_z150_r12_t-15", freqrange=(0kHz, 20kHz), T=Dict([:H => Hs[:, :FIR_M], :x => 0.095, :y => -0.118, :z => 0.150, :r => 12.0, :t => -15.0]), max_mag=6.0, mae_mag=0.6, max_phi=180.0, mae_phi=13.0),

        (src="W/polar_noXO", target="W/polar", freqrange=(0kHz, 2kHz), T=Dict(:H => Hs[:, :FIR_W]), max_mag=0.28, mae_mag=0.04, max_phi=180.0, mae_phi=1.2),
        (src="W/polar_noXO_y-348_r14", target="W/polar_y-348_r14", freqrange=(0kHz, 2kHz), T=Dict(:H => Hs[:, :FIR_W]), max_mag=0.25, mae_mag=0.03, max_phi=180.0, mae_phi=1.4),
        (src="W/polar_noXO", target="W/polar_x-150_y-348_z-35_r-90_t3", freqrange=(0kHz, 2kHz), T=Dict([:H => Hs[:, :FIR_W], :x => -0.15, :y => -0.348, :z => -0.035, :r => -90.0, :t => 3.0]), max_mag=4.5, mae_mag=0.3, max_phi=180.0, mae_phi=1.5),
    ]
    
    for matchFIR in [true, false]
        println("match FIR: $matchFIR")
        for e in test_set
            println("  \"$(e.src)\" -> \"$(e.target)\"")
            src, target = joinpath(test_dir, e.src), joinpath(test_dir, e.target)
            print("\tA:")
            @test test_driver_transform(src, target, e.freqrange, compare_magdb, e.max_mag, e.mae_mag, matchFIR; e.T...)
            print("\tΦ:")
            @test test_driver_transform(src, target, e.freqrange, compare_phasedeg, e.max_phi, e.mae_phi, matchFIR; e.T...)
        end
    end
end