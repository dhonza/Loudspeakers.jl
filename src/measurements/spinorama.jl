export Spinorama, evaluate, power_response, listening_window, early_reflections, project, project_vituix, orbit_evaluator, normlat

struct Spinorama{T<:AngularSpectra}
    hor::T
    ver::T
    step::Int # angle step
    lw_hor # listening window
    lw_ver
    pr_floor # power response
    pr_ceiling
    pr_front_wall
    pr_side_wall
    pr_rear_wall
    pr_rear_wall_standard

    function Spinorama(hor::T, ver::T; 
        lw_hor=30,
        lw_ver=10,
        pr_floor=-40..(-20),
        pr_ceiling=40..60,
        pr_front_wall=-30..30,
        pr_side_wall=[-80..(-40), 40..80],
        pr_rear_wall=[-175..(-90), 90..180],
        pr_rear_wall_standard=[180, -90, 90],
        ) where T
        lw_hor >= 0 || error("hor < 0")
        lw_ver >= 0 || error("ver < 0")

        step_hor = angular_step(hor.angles)
        step_ver = angular_step(ver.angles)
        step_hor ≈ step_ver || error("horizontal step ($(step_hor)) != vetical step ($(step_ver))")

        new{T}(hor, ver, step_hor, 
            lw_hor, lw_ver, 
            pr_floor, pr_ceiling, pr_front_wall, pr_side_wall, pr_rear_wall, pr_rear_wall_standard)
    end
end

function Interpolations.interpolate(spinorama::Spinorama, template::T) where {U, T <: AbstractSampleArray{U}}
    Spinorama(
        interpolate(spinorama.hor, template),
        interpolate(spinorama.ver, template);
        lw_hor=spinorama.lw_hor,
        lw_ver=spinorama.lw_ver,
        pr_floor=spinorama.pr_floor,
        pr_ceiling=spinorama.pr_ceiling,
        pr_front_wall=spinorama.pr_front_wall,
        pr_side_wall=spinorama.pr_side_wall,
        pr_rear_wall=spinorama.pr_rear_wall,
        pr_rear_wall_standard=spinorama.pr_rear_wall_standard,
    ) 
end

function angular_step(angles)
    # expects sorted, normalized to (-180°, 180°] (see normdeg), non-duplicate angles
    # full 360° expected for now
    idx, α = angleidx(angles, 0.0)
    α ≈ 0 || error("missing zero angle!")
    collect(angles) == sort(angles) || error("angles not sorted!")
    step = 360.0 / length(angles) # expected angular step
    for i in 2:length(angles)
        diff = angles[i] - angles[i-1]
        diff ≈ step || error("irregular angle steps! $diff != $step")
    end
    diff = angles[end] + angles[1]
    diff ≈ step || error("irregular angle steps! $diff != $step")
    return step
end

function evaluate(spinorama::Spinorama)
    # note that spinorama results involve magnitude only

    ref = abs.(parent(spinorama.hor[0.0]))
    names!(ref, [:reference])

    pr = power_response(spinorama)    
    lw = listening_window(spinorama)    
    er = early_reflections(spinorama)

    spdi = lw[:, :listening_window] ./ pr[:, :power_response]
    names!(spdi, [:SPDI])

    erdi = lw[:, :listening_window] ./ er[:, :early_reflections]
    names!(erdi, [:ERDI])

    erdi_horizontal = lw[:, :listening_window] ./ er[:, :early_reflections_horizontal]
    names!(erdi_horizontal, [:ERDI_horizontal])

    erdi_vertical = lw[:, :listening_window] ./ er[:, :early_reflections_vertical]
    names!(erdi_vertical, [:ERDI_vertical])

    erdi_standard = lw[:, :listening_window] ./ er[:, :early_reflections_standard]
    names!(erdi_standard, [:ERDI_standard])

    hcat(ref, pr, lw, er, spdi, erdi, erdi_horizontal, erdi_vertical, erdi_standard)
end

function power_response(spinorama::Spinorama)
    # weights according to:
    # Tylka, Joseph G., and Edgar Y. Choueiri. "On the calculation of full and partial directivity indices." 3D Audio and Applied Acoustics Laboratory Technical Report (Princeton University Press, Princeton, NJ) (2014).
    function orbit_weights(orbit)
        angles = orbit.angles
        Δθ = deg2rad(angular_step(angles))

        weights = Float64[]
        for θ in angles
            if θ ≈ 0.0 || θ ≈ 180.0
                w = (1 - cos(Δθ/2))/4 
            else
                θ = deg2rad(θ)
                w = abs(cos(θ + Δθ/2) - cos(θ - Δθ/2))/8
            end
            push!(weights, w)
        end
        weights
    end
    
    pr = sqrt.(+([sum(abs.(parent(orbit)).^2 .* reshape(orbit_weights(orbit), 1, :), dims=2) 
        for orbit in [spinorama.hor, spinorama.ver]]...))
    names!(pr, [:power_response])
    pr
end

function listening_window(spinorama::Spinorama)
    # according to ANSI-CTA-2034-A
    hor, ver, step = spinorama.lw_hor, spinorama.lw_ver, spinorama.step
    
    hor_angle_idxs = [angleidx(spinorama.hor, α).idx for α in -hor:step:hor]
    ver_angle_idxs = [angleidx(spinorama.ver, α).idx for α in -ver:step:ver if α != 0]
    h = parent(spinorama.hor[:, hor_angle_idxs])
    v = parent(spinorama.ver[:, ver_angle_idxs])        
    lw = sqrt.(mean(Float64.(abs.(hcat(h, v)).^2), dims=2)) # RMS
    names!(lw, [:listening_window])
    lw
end

function early_reflections(spinorama::Spinorama)
    # according to ANSI-CTA-2034-A
    # Devantier, Allan. "Characterizing the Amplitude Response of Loudspeaker Systems." Audio Engineering Society Convention 113. Audio Engineering Society, 2002.
    step = spinorama.step

    get_angle_idxs(orbit, α) = [angleidx(orbit, α).idx]

    get_angle_idxs(orbit, αs::AbstractVector) = vcat([get_angle_idxs(orbit, α) for α in αs]...)

    get_angle_idxs(orbit, αs::ClosedInterval) = get_angle_idxs(orbit, [α for α in infimum(αs):step:supremum(αs)])
    
    get_angles(orbit, αs) = parent(orbit)[:, get_angle_idxs(orbit, αs)]
    
    floor_ = get_angles(spinorama.ver, spinorama.pr_floor)
    ceiling_ = get_angles(spinorama.ver, spinorama.pr_ceiling)
    front_wall_ = get_angles(spinorama.hor, spinorama.pr_front_wall)
    side_wall_ = get_angles(spinorama.hor, spinorama.pr_side_wall)
    rear_wall_ = get_angles(spinorama.hor, spinorama.pr_rear_wall)
    rear_wall_standard_ = get_angles(spinorama.hor, spinorama.pr_rear_wall_standard)

    mse_mean(selection...) = mean(Float64.(abs.(hcat(selection...)).^2), dims=2)

    function sqrt_mean(name, mse_means...) 
        er =  sqrt.(mean([mse_means...]))
        names!(er, [name])
        er
    end
    
    floor_mean = mse_mean(floor_)
    ceiling_mean = mse_mean(ceiling_)
    front_wall_mean = mse_mean(front_wall_)
    side_wall_mean = mse_mean(side_wall_)
    rear_wall_mean = mse_mean(rear_wall_)
    rear_wall_mean_standard = mse_mean(rear_wall_standard_)
        
    er = hcat(
        sqrt_mean(:early_reflections, floor_mean, ceiling_mean, front_wall_mean, side_wall_mean, rear_wall_mean),
        sqrt_mean(:early_reflections_horizontal, front_wall_mean, side_wall_mean, rear_wall_mean),
        sqrt_mean(:early_reflections_vertical, floor_mean, ceiling_mean),
        sqrt_mean(:early_reflections_ceiling, ceiling_mean),
        sqrt_mean(:early_reflections_floor, floor_mean),
        sqrt_mean(:early_reflections_front, front_wall_mean),
        sqrt_mean(:early_reflections_rear, rear_wall_mean),
        sqrt_mean(:early_reflections_side, side_wall_mean),
        sqrt_mean(:early_reflections_standard, floor_mean, ceiling_mean, front_wall_mean, side_wall_mean, rear_wall_mean_standard),
        sqrt_mean(:early_reflections_horizontal_standard, front_wall_mean, side_wall_mean, rear_wall_mean_standard),
        sqrt_mean(:early_reflections_rear_standard, rear_wall_mean_standard),
    )
    er
end

# see http://www.movable-type.co.uk/scripts/latlong-vectors.html
# According to Vituix x: left-right (speaker width), y: down-top (heigth), z: front-back (negative z ahead speaker)
# x -> z, y-> x, z -> y
# 
#  y    z
#  |   /
#  |  /
#  | /
#  +----- x
# I use standard coordinates:
# x: back-front (positive x ahead speaker), y: left-right, z: down-top
# 
#   z
#   |
#   |
#   |
#   +----- y
#  /
# x
function checklatlon(lat, lon)
    @assert lat >= -90 && lat <= 90
    @assert lon > -180 && lon <= 180
end

function tovect(lat::T, lon::U) where {T <: Real,U <: Real}
    checklatlon(lat, lon)
    [cosd(lat) * cosd(lon), cosd(lat) * sind(lon), sind(lat)]
end

function tolatlon(x::T, y::U, z::V) where {T <: Real, U <: Real, V <: Real}
    lat = atand(z, sqrt(x^2 + y^2))
    lon = atand(y, x >= 0 ? abs(x) : x) # check for -0.0
    checklatlon(lat, lon)
    lat, lon
end

function tolatlon(vec)
    length(vec) == 3 || throw(ArgumentError("vector length(vec) = $(length) != 3!"))
    tolatlon(vec[1], vec[2], vec[3])
end

# right-hand rotation matrices
Rx(α) = [
    1 0 0;
    0 cosd(α) -sind(α);
    0 sind(α) cosd(α);
]

Ry(α) = [
    cosd(α) 0 sind(α);
    0 1 0;
    -sind(α) 0 cosd(α);
]

Rz(α) = [
    cosd(α) -sind(α) 0;
    sind(α) cosd(α) 0 ;
    0 0 1;
]

function sort_closest_points(orbit_eval, lat, lon, k=2)
    # checklatlon(lat, lon)
    dst(lat2, lon2) = norm(tovect(lat, lon)-tovect(lat2, lon2))
    dists = [dst(p.latlon...) for p in orbit_eval]
    idxs = sortperm(dists)[1:k]
    orbit_eval[idxs], dists[idxs]
end

function orbit_evaluator(spinorama::Spinorama, ::Val{:horizontal})
    [(α=α, latlon=(0, α), f=() -> parent(spinorama.hor)[:, i]) for (i, α) in enumerate(angles(spinorama.hor))], 0.0
end

function normlat(α)
    # normalized latitude, longitude for the vertical orbit
    @assert -180 < α && α <= 180
    α > 90 && return 180 - α, 180
    α < -90 && return -180 - α, 180
    return α, 0
end
function orbit_evaluator(spinorama::Spinorama, ::Val{:vertical})
    [(α=α, latlon=normlat(α), f=() -> parent(spinorama.ver)[:, i]) for (i, α) in enumerate(angles(spinorama.ver))], 0.0
end

function orbit_evaluator(spinorama::Spinorama, lat, lon)
    function latlon_to_orbit_angle(lat, lon)
        # checklatlon(lat, lon)
        x, y, z = tovect(lat, lon)
        p = [0, y, z] # projection to n = [1, 0, 0] , i.e, Y-Z plane
        p = p[2] / norm(p) # = sum(p .* [0, 1, 0]) / norm(p)
        sign(z) * acosd(p)
    end
    β = latlon_to_orbit_angle(lat, lon)
    b = abs(β)
    
    w =  b >= 90 ? [b-90, 180-b] : [90-b, b]        
    
    # change normal ordering of angles
    npts = nchannels(spinorama.hor)
    
    hidx(i) = -90 < b < 90 ? i : (i < npts ? npts - i : npts)
    vidx(i) = β < 0 ? (i < npts ? npts - i : npts) : i

    # hidx(i) = i
    # vidx(i) = (i < npts ? npts - i : npts)
    
    orbit = [(α=α, latlon=tolatlon(Rx(β) * tovect(0, α)),
            f=() -> mean(
                        hcat(
                            parent(spinorama.hor)[:, hidx(i)],
                            parent(spinorama.ver)[:, vidx(i)]
                        ), weights(w)
                    )
        )
        for (i, α) in enumerate(angles(spinorama.hor))]
    orbit, β
end

function interpolation_orbit_evaluator(spinorama::Spinorama, lat, lon)
    # checklatlon(lat, lon)
    lat ≈ 0 && return orbit_evaluator(spinorama, Val(:horizontal))
    lon ≈ 0 && return orbit_evaluator(spinorama, Val(:vertical))
    return orbit_evaluator(spinorama, lat, lon)
end

function project(spinorama::Spinorama, lat, lon)
    # TODO: don't have to construct whole orbit to find its two closest points
    orbit_eval, β = interpolation_orbit_evaluator(spinorama, lat, lon)
    # @show β
    selpts, dists = sort_closest_points(orbit_eval, lat, lon, 2)
    pts = hcat([o.f() for o in selpts]...)
    # pts = selpts[1].f()
    # pts = selpts[1].f()
    # return pts

    mean(pts, weights([dists[2], dists[1]]))
    # rms_mean(pts, weights([dists[2], dists[1]]))
end

function project(spinorama::Spinorama; lat=0.0, lon=0.0, ldist=2.0, H=nothing, x=0.0, y=0.0, z=0.0, r=0.0, t=0.0, c=344.0)
    # meters, degrees, speed of sound in m/s
    x, y, z, ldist = tom.((x, y, z, ldist))
    c = tomps(c)
    driver_center = [x, y, z]
    # @show driver_center

    target = tovect(lat, lon) .* ldist
    tvec = target - driver_center
    # @show tvec, norm(tvec)

    # driver_trans = Ry(t) * Rz(-r)
    # driver_invtrans = Rz(-r) * Ry(t)
    driver_invtrans = Ry(t) * Rz(-r) # OK

    tvec_invtrans = driver_invtrans * tvec
    # @show tvec_invtrans, norm(tvec_invtrans)

    lat_, lon_ = tolatlon(tvec_invtrans)
    # @show lat_, lon_

    tdist = norm(tvec)
    dratio = ldist/tdist
    # @show dratio

    response = project(spinorama, lat_, lon_) .* dratio
    delay_dist = (tdist - ldist)
    # @show delay_dist/c
    delay!(response, delay_dist / c)
    names!(response, [:response])
    
    if isnothing(H)
        return response
    else
        return H .* response
    end
end

# Vituix axes: x -> z, y-> x, z -> y
function project_vituix(spinorama::Spinorama; lat=0.0, lon=0.0, ldist=2.0, H=nothing, x=0.0, y=0.0, z=0.0, r=0.0, t=0.0, c=344.0)
    project(spinorama; lat=lat, lon=lon, ldist=ldist, H=H, x=-z, y=x, z=y, r=r, t=t, c=c)
end

# for now
const DriverDefinition = Dict{Symbol, Any}

function project(drivers::DriverDefinition...; lat=0.0, lon=0.0, ldist=2.0)
    Y = [project_vituix(d[:spinorama]; lat=lat, lon=lon, ldist=ldist, d[:transform]...) for d in drivers]
    sum(hcat(Y...); dims=2)
end

function project_vituix(drivers::DriverDefinition...; lat=0.0, lon=0.0, ldist=2.0)
    Y = [project_vituix(d[:spinorama]; lat=lat, lon=lon, ldist=ldist, d[:transform]...) for d in drivers]
    sum(hcat(Y...); dims=2)
end
