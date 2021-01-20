export miccalibrate, miccalibrate!, load_angular_measurements, load_vituix_spinorama_curves, load_spinorama_from_FRDs, load_vituix_FIRs

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

Load all directivity spectra in `irdir` with each spectrum in a separate IR file/directory matching regex `template`. 

Use `loader(f)` function to get the a subtype of the `AbstractSpectrumArray`. Imported spectra can be limited to only `:positive` or `:negative` angles using `select`. 
Degrees are infered as to be the first group of `ir_template`. Example: `r"ir_(\\d\\d\\d).wav"`.
"""
function load_angular_measurements(irdir, template, loader, select=:all)
    @assert select ∈ [:all, :positive, :negative]
    spectra = []
    angles = Float64[]
    for f in readdir(irdir)
        fpath = joinpath(irdir, f)
        m = match(template, f)
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
    length(spectra) > 0 || error("no file imported!")
    catspectra = hcat(spectra...)
    AngularMeasurements(catspectra, angles)
end

function load_vituix_spinorama_curves(rootdir; samplerate=44100*Unitful.Hz)
    function load(file_name, channel_name)
        readfrd(joinpath(rootdir, file_name), samplerate, magf=db2amp, name=channel_name) .|> abs
    end
    pairs = [
        "VituixCAD_CTA-2034-A Directivity Index.frd" => :SPDI,
        "VituixCAD_CTA-2034-A ER Ceiling.frd" => :early_reflections_ceiling,
        "VituixCAD_CTA-2034-A ER Floor.frd" => :early_reflections_floor,
        "VituixCAD_CTA-2034-A ER Front.frd" => :early_reflections_front,
        "VituixCAD_CTA-2034-A ER Horizontal.frd" => :early_reflections_horizontal,
        "VituixCAD_CTA-2034-A ER Rear.frd" => :early_reflections_rear,
        "VituixCAD_CTA-2034-A ER Side.frd" => :early_reflections_side,
        "VituixCAD_CTA-2034-A ER Total.frd" => :early_reflections,
        "VituixCAD_CTA-2034-A ER Vertical.frd" => :early_reflections_vertical,
        "VituixCAD_CTA-2034-A ERDI Horizontal.frd" => :ERDI_horizontal,
        "VituixCAD_CTA-2034-A ERDI Total.frd" => :ERDI,
        "VituixCAD_CTA-2034-A ERDI Vertical.frd" => :ERDI_vertical,
        "VituixCAD_CTA-2034-A In-room response.frd" => :in_room,
        "VituixCAD_CTA-2034-A Listening window.frd" => :listening_window,
        "VituixCAD_CTA-2034-A Power response.frd" => :power_response,
        "VituixCAD_CTA-2034-A Reference angle.frd" => :reference,
    ]
    hcat(map(p -> load(p...), pairs)...)
end

function load_spinorama_from_FRDs(rootdir; 
    template_hor=Regex(raw"VituixCAD_PolarFR hor (.*)\.frd"), template_ver=Regex(raw"VituixCAD_PolarFR ver (.*)\.frd"), samplerate=44100*Unitful.Hz)
    loader = f->readfrd(f, samplerate, magf=db2amp)
    hor = load_angular_measurements(rootdir, template_hor, loader, :all)
    ver = load_angular_measurements(rootdir, template_ver, loader, :all)
    return Spinorama(hor, ver)
end

function load_vituix_FIRs(wavs...; Δt=0.0, names=nothing)
    if !isnothing(names) && length(wavs) != length(names)
        throw(ArgumentError("the number of WAV files does not match the number of names!"))
    end
    names_ = isnothing(names) ? [Symbol("FIR_$i") for i in 1:length(wavs)] : names
    Hs = hcat((readwav_rfft(wav)[2] for wav in wavs)...)
    names!(Hs, names_)
    Hmax = maximum(abs.(Hs))
    Δt ≈ 0.0 || delay!(Hs, Δt)
    return Hs ./ Hmax
end