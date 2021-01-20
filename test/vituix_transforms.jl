@testset "Transforms" begin
    println("VITUIX TRANSFORM TESTS")

    test_set = [
        (src="T/polar", target="T/polar_r20", freqrange=(1kHz, 20kHz), T=Dict(:r => 20.0), max_mag=3.2, mae_mag=0.4, max_phi=180.0, mae_phi=6.5),
        (src="T/polar", target="T/polar_r10_t15", freqrange=(1kHz, 20kHz), T=Dict(:r => 10.0, :t => 15.0), max_mag=4.7, mae_mag=0.3, max_phi=180.0, mae_phi=8.5),
        (src="T/polar_noXO", target="T/polar_noXO_r20", freqrange=(1kHz, 20kHz), T=Dict(:r => 20.0), max_mag=3.2, mae_mag=0.4, max_phi=180.0, mae_phi=6.3),
        (src="T/polar_noXO", target="T/polar_noXO_r-80_t-80", freqrange=(1kHz, 20kHz), T=Dict(:r => -80.0, :t => -80.0), max_mag=8.5, mae_mag=0.4, max_phi=180.0, mae_phi=6.0),
        
        (src="M/polar", target="M/polar_r10", freqrange=(0kHz, 20kHz), T=Dict(:r => 10.0), max_mag=4.3, mae_mag=0.25, max_phi=180.0, mae_phi=3.2),
        (src="M/polar", target="M/polar_r170", freqrange=(0kHz, 20kHz), T=Dict(:r => 170.0), max_mag=5.0, mae_mag=0.25, max_phi=180.0, mae_phi=4.0),
        (src="M/polar", target="M/polar_y-118", freqrange=(0kHz, 20kHz), T=Dict(:y => -0.118), max_mag=13.0, mae_mag=0.25, max_phi=180.0, mae_phi=0.7),
        (src="M/polar", target="M/polar_y-118_r170", freqrange=(0kHz, 20kHz), T=Dict(:y => -0.118, :r => 170.0), max_mag=13.0, mae_mag=0.25, max_phi=180.0, mae_phi=2.2),
        (src="M/polar", target="M/polar_y-118_r180", freqrange=(0kHz, 20kHz), T=Dict(:y => -0.118, :r => 180.0), max_mag=13.0, mae_mag=0.25, max_phi=180.0, mae_phi=2.2),
        (src="M/polar", target="M/polar_x95_y-118_z150_r12_t-15", freqrange=(0kHz, 20kHz), T=Dict(:x => 0.095, :y => -0.118, :z => 0.150, :r => 12.0, :t => -15.0), max_mag=6.0, mae_mag=0.3, max_phi=180.0, mae_phi=4.0),
        (src="M/polar_noXO", target="M/polar_noXO_r170", freqrange=(0kHz, 20kHz), T=Dict(:r => 170.0), max_mag=5.0, mae_mag=0.25, max_phi=180.0, mae_phi=4.0),
        (src="M/polar_noXO", target="M/polar_noXO_r-170", freqrange=(0kHz, 20kHz), T=Dict(:r => -170.0), max_mag=5.0, mae_mag=0.25, max_phi=180.0, mae_phi=3.5),
        (src="M/polar_noXO", target="M/polar_noXO_y-118_r170", freqrange=(0kHz, 20kHz), T=Dict(:y => -0.118, :r => 170.0), max_mag=13.0, mae_mag=0.25, max_phi=180.0, mae_phi=2.2),
        
        (src="W/polar", target="W/polar_y-348", freqrange=(0kHz, 2kHz), T=Dict(:y => -0.348), max_mag=3.5, mae_mag=0.15, max_phi=180.0, mae_phi=1.2),
        (src="W/polar", target="W/polar_y-348_r14", freqrange=(0kHz, 2kHz), T=Dict(:y => -0.348, :r => 14.0), max_mag=5.6, mae_mag=0.15, max_phi=180.0, mae_phi=1.2),
        (src="W/polar", target="W/polar_x-150_y-348_z-35_r-90_t3", freqrange=(0kHz, 2kHz), T=Dict(:x => -0.15, :y => -0.348, :z => -0.035, :r => -90.0, :t => 3.0), max_mag=4.5, mae_mag=0.15, max_phi=180.0, mae_phi=1.5),
        (src="W/polar_noXO", target="W/polar_noXO_y-348_r14", freqrange=(0kHz, 2kHz), T=Dict(:y => -0.348, :r => 14.0), max_mag=5.6, mae_mag=0.15, max_phi=180.0, mae_phi=1.2),
    ]

    for e in test_set
        println("  \"$(e.src)\" -> \"$(e.target)\"")
        src, target = joinpath(test_dir, e.src), joinpath(test_dir, e.target)
        print("\tA:")
        @test test_driver_transform(src, target, e.freqrange, compare_magdb, e.max_mag, e.mae_mag, false; e.T...)
        print("\tΦ:")
        @test test_driver_transform(src, target, e.freqrange, compare_phasedeg, e.max_phi, e.mae_phi, false; e.T...)
    end
end