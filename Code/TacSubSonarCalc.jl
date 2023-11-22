#!/bin/julia
# TacSubSonarCalc.jl
# 2023-09-17
# curtis
# Sonar equation calculations for a tactical submarine game

# ArgParse: A package for handling command line arguments
using ArgParse

# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION -------------------------

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "activetablecsv"
            arg_type = String
            help = "File where active sonar table is saved (a CSV)"
            required = true
        "passivetablecsv"
            arg_type = String
            help = "File where passive sonar table is saved (a CSV)"
            required = true
        "--nlmean", "-n"
            arg_type = Float64
            help = "Mean noise level"
            default = 72.0
        "--nlsd", "-s"
            arg_type = Float64
            help = "Noise level standard deviation"
            default = 10.0
        "--maxrange", "-R"
            arg_type = Float64
            help = "Maximum range tracked for calculations"
            default = 30.0
        "--rangeinc", "-r"
            arg_type = Float64
            help = "Increments of range at which calculations are made"
            default = 4.0
        "--ddepthshallow", "-d"
            arg_type = Float64
            help = "Depth of the detector when at shallow depth"
            default = 200.0
        "--ddepthdeep", "-D"
            arg_type = Float64
            help = "Depth of the detector when at deep depth"
            default = 1200.0
        "--minangle", "-a"
            arg_type = Float64
            help = "Minimum angle considered for sound propagation"
            default = -10.0
        "--maxangle", "-A"
            arg_type = Float64
            help = "Maximum angle considered for sound propagation"
            default = 10.0
        "--stepangle", "-t"
            arg_type = Float64
            help = "Increments made to angle"
            default = 0.5
        "--emtdepthshallow", "-e"
            arg_type = Float64
            help = "Emitter depth when shallow"
            default = 260.0
        "--emtdepthdeep", "-E"
            arg_type = Float64
            help = "Emitter depth when deep"
            default = 1210.0
        "--maxdepth", "-X"
            arg_type = Float64
            help = "Maximum depth considered for ocean sound propagation"
            default = 18000.0
        "--freq", "-f"
            arg_type = Float64
            help = "Frequency of sound considered"
            default = 150.0
        "--svpstep", "-v"
            arg_type = Float64
            help = "SVP incrememnts for raytracing"
            default = 1.0
        "--dtpassive", "-P"
            arg_type = Float64
            help = "Detection threshold for passive detection"
            default = 15.0
        "--dipassive", "-i"
            arg_type = Float64
            help = "Directivity index for passive sonar; if not set, elements and spacing may be used for setting up a line array"
            default = nothing
        "--elements", "-l"
            arg_type = Float64
            help = "Number of elements of line array sonar used for passive sonar"
            default = nothing
        "--spacing", "-c"
            arg_type = Float64
            help = "Spacing of elements of line array sonar used for passive sonar"
            default = nothing
        "--slcreep", "-C"
            arg_type = Float64
            help = "Source level at creep speed"
            default = 110.0
        "--slslow", "-S"
            arg_type = Float64
            help = "Source level at slow speed"
            default = 120.0
        "--slfast", "-F"
            arg_type = Float64
            help = "Source level at fast speed"
            default = 130.0
        "--slflank", "-L"
            arg_type = Float64
            help = "Source level at flank speed"
            default = 140.0
        "--tssub", "-u"
            arg_type = Float64
            help = "Target strength of submarine target"
            default = 15.0
        "--tssurf", "-U"
            arg_type = Float64
            help = "Target strength of surface target"
            default = 25.0
        "--slactive", "-V"
            arg_type = Float64
            help = "Sound level of active sonar"
            default = 210.0
        "--dtactive", "-T"
            arg_type = Float64
            help = "Detection threshold for active sonar"
            default = 50.0
        "--pistondiameter", "-o"
            arg_type = Float64
            help = "Piston diameter for piston sonar, overridden by --diactive"
            default = nothing
        "--diactive", "-I"
            arg_type = Float64
            help = "Directivity index of active sonar; if not set, piston sonar used"
            default = nothing
        "--svpcsv", "-Y"
            arg_type = Float64
            help = "CSV file for SVP; if not set, default SVP used"
            default = nothing
    end

    return Dict([(Symbol(key), val) for (key, val) in parse_args(s)])
end

if !isinteractive()
    parsed_args = parse_commandline()
end

# PACKAGES ---------------------------------------------------------------------

using SimpleSonar
using Distributions
using DataFrames
using Interpolations
using Plots
using StatsPlots
using DataFramesMeta
using Chain
using Statistics
using CSV
using Tables

# EXECUTABLE SCRIPT MAIN FUNCTIONALITY -----------------------------------------

"""
See [`parse_commandline`](@parse_commandline) for argument description
"""
function main(;
              activetablecsv::String,
              passivetablecsv::String,
              nlmean::Real = 72.0,
              nlsd::Real = 10.0,
              maxrange::Real = 30.0,
              rangeinc::Real = 4.0,
              ddepthshallow::Real = 200.0,
              ddepthdeep::Real = 1200.0,
              minangle::Real = -10.0,
              maxangle::Real = 10.0,
              stepangle::Real = 0.5,
              emtdepthshallow::Real = 260.0,
              emtdepthdeep::Real = 1210.0,
              maxdepth::Real = 18000.0,
              freq::Real = 150.0,
              svpstep::Real = 1.0,
              dtpassive::Real = 15.0,
              dipassive::Union{Nothing, Real} = nothing,
              elements::Union{Nothing, Integer} = nothing,
              spacing::Union{Nothing, Real} = nothing,
              slcreep::Real = 110.0,
              slslow::Real = 120.0,
              slfast::Real = 130.0,
              slflank::Real = 140.0,
              tssub::Real = 15.0,
              tssurf::Real = 25.0,
              slactive::Real = 210.0,
              dtactive::Real = 50.0,
              pistondiameter::Union{Nothing, Real} = nothing,
              diactive::Union{Nothing, Real} = nothing,
              svpcsv::Union{Nothing, String} = nothing)
    if !isnothing(svpcsv)
        svpmat = Tables.matrix(CSV.File(svpcsv))
    else
        svpmat = [# depth    velocity 
                        0      1540.4 
                       10      1540.5 
                       20      1540.7 
                       30      1534.4 
                       50      1523.3 
                       75      1519.6 
                      100      1518.5 
                      125      1517.9 
                      150      1517.3 
                      200      1516.6 
                      250      1516.5 
                      300      1516.2 
                      400      1516.4 
                      500      1517.2 
                      600      1518.2 
                      700      1519.5 
                      800      1521.0 
                      900      1522.6 
                     1000      1524.1 
                     1100      1525.7 
                     1200      1527.3 
                     1300      1529.0 
                     1400      1530.7 
                     1500      1532.4 
                     1750      1536.7 
                     2000      1541.0 
                 ].* 3.28084
    end

    # Common parameters
    common_params = Dict(
        "nlmean"           => nlmean,
        "nlsd"             => nlsd,
        "drange"           => [0:rangeinc:maxrange...].* 6000.0,
        "ddepth_shallow"   => ddepthshallow,
        "ddepth_deep"      => ddepthdeep,
        "minangle"         => minangle,
        "maxangle"         => maxangle,
        "stepangle"        => stepangle,
        "emtdepth_shallow" => emtdepthshallow,
        "emtdepth_deep"    => emtdepthdeep,
        "maxdepth"         => maxdepth,
        "freq"             => freq,
        "svpmat"           => svpmat,
        "svpstep"          => svpstep,
        "maxrange"         => maxrange * 6000
    )

    svp_obj = svp(common_params["svpmat"][:,1], common_params["svpmat"][:,2])
    plot(svp_obj.velocity, svp_obj.depth, yflip = true, legend = false)
    svp_itp = svp_to_interpolation(svp_obj, Gridded(Linear()), Linear())
    velocity = svp_itp(common_params["ddepth_shallow"])

    fine_svp_obj = svp_refine(svp_obj, max_depth = common_params["maxdepth"])
    plot(fine_svp_obj.velocity, fine_svp_obj.depth, yflip = true,
         legend = false)

    wavelength = freq_to_wavelength(common_params["freq"], velocity)

    if isnothing(dipassive)
        if !isnothing(elements) && !isnothing(spacing)
            dipassive = line_di(Unsigned(elements), spacing, wavelength)
        elseif isnothing(elements) && isnothing(spacing)
            dipassive = line_di(Unsigned(100), 0.5 / 12, wavelength) + 40
        else
            throw(ArgumentError("Cannot have dipassive be nothing and only one of elements and spacing be nothing"))
        end
    end

    if isnothing(diactive)
        if !isnothing(pistondiameter)
            diactive = piston_di(pistondiameter, wavelength)
        else
            diactive = piston_di(18.0, wavelength)
        end
    end

    # Passive parameters
    passive_params = Dict(
        "di"       => dipassive,
        "dt"       => dtpassive,
        "sl_creep" => slcreep,
        "sl_slow"  => slslow,
        "sl_fast"  => slfast,
        "sl_flank" => slflank
    )
    
    # Active parameters
    active_params = Dict(
        "di"      => diactive,
        "ts_sub"  => tssub,
        "ts_surf" => tssurf,
        "sl"      => slactive,
        "dt"      => dtactive
    )

    d6 = DiscreteUniform(1, 6)
    mean_2d6 = 2 * mean(d6)
    sd_2d6 = sqrt(2) * std(d6)
    
    angles = [(common_params["minangle"]:common_params["stepangle"]:common_params["maxangle"])...].*π./180
    ray_combined_result_shallow = raytrace_angle_df(fine_svp_obj,
                                                    common_params["emtdepth_shallow"],
                                                    angles,
                                                    common_params["maxrange"])
    ray_combined_shallow_df = ray_combined_result_shallow[:ray]
    @df ray_combined_shallow_df plot(:range / 6000, :depth, group = :angle,
                                     yflip = true,
                                     legend=false)

    ray_combined_result_deep = raytrace_angle_df(fine_svp_obj,
                                                 common_params["emtdepth_deep"],
                                                 angles,
                                                 common_params["maxrange"])
    ray_combined_deep_df = ray_combined_result_deep[:ray]
    @df ray_combined_deep_df plot(:range / 6000, :depth, group=:angle,
                                  yflip=true,
                                  legend=false, xformatter=:plain)

    tl = vcat(
        DataFrame(detector="shallow", emitter = "shallow",
                  tl = ray_df_to_tl.([ray_combined_shallow_df],
                                          common_params["drange"],
                                          [common_params["ddepth_shallow"]]),
                  range = common_params["drange"]),
        DataFrame(detector="shallow", emitter = "deep",
                  tl = ray_df_to_tl.([ray_combined_deep_df],
                                      common_params["drange"],
                                      [common_params["ddepth_shallow"]]),
                  range = common_params["drange"]),
        DataFrame(detector="deep", emitter = "shallow",
                  tl = ray_df_to_tl.([ray_combined_shallow_df],
                                      common_params["drange"],
                                      [common_params["ddepth_deep"]]),
                  range = common_params["drange"]),
        DataFrame(detector="deep", emitter = "deep",
                  tl = ray_df_to_tl.([ray_combined_deep_df],
                                      common_params["drange"],
                                      [common_params["ddepth_deep"]]),
                  range = common_params["drange"]),
    )
    
    passive_speed_decoder = Dict("creep" => passive_params["sl_creep"],
                             "slow" => passive_params["sl_slow"],
                             "fast" => passive_params["sl_fast"],
                             "flank" => passive_params["sl_flank"])
    active_source_decoder = Dict("sub"  => active_params["ts_sub"],
                                 "surf" => active_params["ts_surf"])
    detection_table_passive = crossjoin(DataFrame(speed = ["creep", "slow",
                                                          "fast", "flank"]),
                                        tl)
    detection_table_active = crossjoin(DataFrame(source = ["sub", "surf"]),
                                        tl)
    detection_table_passive = @chain detection_table_passive begin
        @transform(:se=sonar_passive.(getindex.(Ref(passive_speed_decoder),
                                                :speed),
                                      :tl,
                                      Ref(Normal(common_params["nlmean"],
                                                 common_params["nlsd"])),
                                      Ref(passive_params["di"]),
                                      Ref(passive_params["dt"])))
    end
    detection_table_passive.se_threshold = [sonar_threshold(x) for x in detection_table_passive.se]
    detection_table_passive.detection_prob = [detection_prob(x) for x in detection_table_passive.se]
    detection_table_active = @chain detection_table_active begin
        @transform(:se=sonar_noise.(Ref(active_params["sl"]),
                                    :tl,
                                    getindex.(Ref(active_source_decoder),
                                        :source),
                                    Ref(Normal(common_params["nlmean"],
                                               common_params["nlsd"])),
                                    Ref(active_params["di"]),
                                    Ref(active_params["dt"])))
    end
    detection_table_active.se_threshold = [sonar_threshold(x) for x in detection_table_active.se]
    detection_table_active.detection_prob = [detection_prob(x) for x in detection_table_active.se]
    
    detection_table_passive.raw_modifier = sd_2d6/common_params["nlsd"] .* (
        detection_table_passive.se_threshold .- common_params["nlmean"])
    detection_table_active.raw_modifier = sd_2d6/common_params["nlsd"] .* (
        detection_table_active.se_threshold .- common_params["nlmean"])

    detection_threshold = DataFrame(
        class = ["passive", "active"],
        threshold = (sd_2d6/common_params["nlsd"]) .*
            [passive_params["dt"], active_params["dt"]] .+ mean_2d6)

    detection_threshold_overall = mean(detection_threshold.threshold)
    detection_threshold_adjust = detection_threshold.threshold .-
        mean(detection_threshold.threshold)
    detection_table_passive.modifier = detection_table_passive.raw_modifier .-
        detection_threshold_adjust[1]
    detection_table_active.modifier = detection_table_active.raw_modifier .-
        detection_threshold_adjust[2]

    CSV.write(activetablecsv, detection_table_passive)
    CSV.write(passivetablecsv, detection_table_active)
end

if !isinteractive()
    main(; parsed_args...)

    exit()
end
