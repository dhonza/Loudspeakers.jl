export miccalibrate, miccalibrate!, load_angular_measurements

function miccalibrate!(a::AbstractSpectrumArray, calfile; skipto=0)
    # df = CSV.read(calfile, DataFrame; header=false, skipto=skipto)
    df = CSV.read(calfile, DataFrame; header=false)
    freqs, magdb = df.Column1, df.Column2
    mag = LinearInterpolation(freqs, db2amp.(magdb); extrapolation_bc = Flat())(domain(a))
    a ./ Complex.(mag)
end

miccalibrate(a::AbstractSpectrumArray, calfile; skipto=0) = miccalibrate!(deepcopy(a), calfile; skipto=skipto)

"""
    load_directivity_spectra(rootdir, dirtemplate, loader, select=:all)

Load all directivity spectra in `irdir` with each spectrum in a separate IR file/directory matching regex `ir_template`. 

Use `loader(f)` function to get the a subtype of the `AbstractSpectrumArray`. Imported spectra can be limited to only `:positive` or `:negative` angles using `select`. 
Degrees are infered as to be the first group of `ir_template`. Example: `r"ir_(\\d\\d\\d).wav"`.
"""
function load_angular_measurements(irdir, ir_template, loader, select=:all)
    @assert select ∈ [:all, :positive, :negative]
    spectra = []
    angles = Float64[]
    for f in readdir(irdir)
        fpath = joinpath(irdir, f)
        m = match(ir_template, f)
        if !isnothing(m)
            angle = normdeg(parse(Float64, m.captures[1]))
            if select == :positive && angle < 0
                continue
            end
            if select == :negative && angle > 0
                continue
            end
            spectrum = loader(fpath)
            names!(spectrum, [Symbol(round(angle) ≈ angle ? round(Int, angle) : angle)])
            push!(spectra, spectrum)
            push!(angles, angle)
        end
    end
    length(spectra) > 0 || error("now file imported!")
    catspectra = hcat(spectra...)
    AngularMeasurements(catspectra, angles)
end