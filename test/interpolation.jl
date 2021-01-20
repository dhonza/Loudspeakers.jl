@testset "Interpolation" begin
    println("INTERPOLATION TEST")

    test_set = [
        (input="T/polar", freqrange=(1kHz, 12kHz), max_mag=30.0, mae_mag=5.0, max_phi=180.0, mae_phi=60.0),
        (input="T/polar_noXO", freqrange=(1kHz, 12kHz), max_mag=30.0, mae_mag=5.0, max_phi=180.0, mae_phi=60.0),
        
        (input="M/polar", freqrange=(200Hz, 12kHz), max_mag=35.0, mae_mag=3.0, max_phi=180.0, mae_phi=35.0),
        (input="M/polar_noXO", freqrange=(200Hz, 12kHz), max_mag=35.0, mae_mag=3.0, max_phi=180.0, mae_phi=35.0),
        
        (input="W/polar", freqrange=(0Hz, 1kHz), max_mag=2.0, mae_mag=0.5, max_phi=20.0, mae_phi=3.0),
        (input="W/polar_noXO", freqrange=(0Hz, 1kHz), max_mag=2.0, mae_mag=0.5, max_phi=20.0, mae_phi=3.0),
    ]
    # NOTE: errors are quite large, mainly for higher frequencies
    for e in test_set
        println("  \"$(e.input)\"")
        input_dir = joinpath(test_dir, e.input)
        print("\tA:")
        @test test_interpolation(input_dir, e.freqrange, compare_magdb, e.max_mag, e.mae_mag)
        print("\tΦ:")
        @test test_interpolation(input_dir, e.freqrange, compare_phasedeg, e.max_phi, e.mae_phi)
    end
end