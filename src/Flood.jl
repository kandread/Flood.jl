module Flood

# constants
const g = 9.81
const alpha = 0.7
const depth_thresh = 0.001

# boundary conditions
abstract type BC end

struct HFIX <: BC
    values :: Array{Real, 1}
end

const QFIX = HFIX

const FREE = HFIX

struct QVAR <: BC
    values :: Array{Real, 1}
    time :: Array{Int, 1}
    prev_time :: Int
end

const HVAR = QVAR

# domain data structure
struct Domain
    xul :: Float32
    yul :: Float32
    xres :: Float32
    yres :: Float32
    nrows :: Int
    ncols :: Int
end

# parameter structure
struct Params
    dem_file :: String
    sim_time :: Float32
    init_tstep :: Float32
    fpfric :: Float32
    saveint :: Float32
    bci_file :: String
    bdy_file :: Nullable{String}
    output_dir :: Nullable{String}

end

function read_params(filename::String)
    data = Dict{String, String}()
    open(filename, "r") do f
        for line in eachline(f)
            if !startswith(line, "#")
                tokens = split(line)
                if length(tokens) > 1
                    data[tokens[1]] = tokens[2]
                end
            end
        end
    end
    dem_file = haskey(data, "DEMfile") ? data["DEMfile"] : nothing
    sim_time = haskey(data, "sim_time") ? parse(Float32, data["sim_time"]) : nothing
    init_tstep = haskey(data, "initial_tstep") ? parse(Float32, data["initial_tstep"]) : nothing
    fpfric = haskey(data, "fpfric") ? parse(Float32, data["fpfric"]) : nothing
    bci_file = haskey(data, "bcifile") ? data["bcifile"] : nothing
    bdy_file = haskey(data, "bdyfile") ? data["bdyfile"] : nothing
    saveint = haskey(data, "saveint") ? parse(Float32, data["saveint"]) : nothing
    outdir = haskey(data, "dirroot") ? data["dirroot"] : nothing
    return Params(dem_file, sim_time, init_tstep, fpfric, saveint, bci_file, bdy_file)
end

end
