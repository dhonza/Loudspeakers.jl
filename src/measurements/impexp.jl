export miccalibrate, miccalibrate!, load_directivity_spectra

function miccalibrate!(a::AbstractSpectrumArray, calfile; skipto=0)
    df = CSV.read(calfile; header=false, skipto=skipto)
    freqs, magdb = df.Column1, df.Column2
    mag = LinearInterpolation(freqs, db2amp.(magdb); extrapolation_bc = Flat())(domain(a))
    a ./ Complex.(mag)
end

miccalibrate(a::AbstractSpectrumArray, calfile; skipto=0) = miccalibrate!(deepcopy(a), calfile; skipto=skipto)

"""
    load_directivity_spectra(rootdir, dirtemplate, loader, select=:all)

Load all directivity spectra in `rootdir` with each spectrum in a separate subdirectory matching regex `dirtemplate`. 

Use `loader(dir)` function to get the `SampleArray`. Imported spectra can be limited to only `:positive` or `:negative` angles using `select`. 
Degrees are infered as to be the first group of `dirtemplate`. Example: `r"LEFT_SPEAKER_(\\d\\d\\d)_RATE48000"`.
"""
function load_directivity_spectra(rootdir, dirtemplate, loader, select=:all)
    @assert select ∈ [:all, :positive, :negative]
    spectra = []
    angles = Float64[]
    for f in readdir(rootdir)
        fpath = joinpath(rootdir, f)
        if isdir(fpath)
            m = match(dirtemplate, f)
            if !isnothing(m)
                angle = normdeg(parse(Float64, m.captures[1]))
                if select == :positive && angle < 0
                    continue
                end
                if select == :negative && angle > 0
                    continue
                end
                push!(spectra, loader(fpath))
                push!(angles, angle)
            end
        end
    end
    DirectivitySpectra(hcat(spectra...), angles)
end