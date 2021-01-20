
using Test
using DSP
using LinearAlgebra
using Statistics
using Unitful
using Unitful: ns, ms, µs, s, Hz, kHz, MHz, GHz, THz, m
using SampleArrays
using Loudspeakers


test_dir = "../data/test/vituix"

@testset "Loudspeakers.jl" begin
    include("utils.jl")

    # check interpolation on orbitals, e.g., 10° and 20° to estimate 15°
    include("interpolation.jl")

    # test spinorama curves (listening window, power response, ...)
    include("vituix_spinorama.jl")

    # test single driver transforms (translation, rotation, tilt, listening distance)
    # TODO: different listening distances
    include("vituix_transforms.jl")

    # test FIR filters applied to single drivers
    include("vituix_fir.jl")

    # test multiple drivers
    # TODO: check different listening distances - it seems there some differences (try even < 1m and < 2m)
    include("vituix_multi_driver.jl")

    # TODO: projection point transforms? How to get data from Vituix? Use reference angle, or Room?
end
