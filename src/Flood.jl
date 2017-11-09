module Flood

import GDAL

# constants
const g = 9.81
const alpha = 0.7
const depth_thresh = 0.001

# grid types
@enum Grid H QX QY

# boundary conditions
abstract type BoundaryCondition end

struct HFIX <: BoundaryCondition
    values :: Array{Real, 1}
end

const QFIX = HFIX

const FREE = HFIX

struct QVAR <: BoundaryCondition
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

function read_domain(filename::String)
    GDAL.registerall()
    ds = GDAL.open(filename, GDAL.GA_ReadOnly)
    nrows = GDAL.getrasterysize(ds)
    ncols = GDAL.getrasterxsize(ds)
    gt = zeros(6)
    GDAL.getgeotransform(ds, gt)
    dom = Domain(gt[1], gt[4], gt[2], gt[6], nrows, ncols)
    return dom
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
    return Params(dem_file, sim_time, init_tstep, fpfric, saveint, bci_file, bdy_file, outdir)
end

# raster input and output
function read_raster(filename::String)
    GDAL.registerall()
    ds = GDAL.open(filename, GDAL.GA_ReadOnly)
    nrows = GDAL.getrasterysize(ds)
    ncols = GDAL.getrasterxsize(ds)
    band = GDAL.getrasterband(ds, 1)
    dtype = GDAL.getrasterdatatype(band)
    jtype = eval(parse(GDAL.getdatatypename(GDAL.getrasterdatatype(band))))
    data = Array{jtype}(nrows*ncols)
    GDAL.rasterio(band, GDAL.GF_Read, 0, 0, ncols, nrows, data, ncols, nrows, dtype, 0, 0)
    GDAL.close(ds)
    return data
end

function write_raster(filename::String, data::Array{Float32}, dom::Domain, edges::Grid)
    GDAL.registerall()
    gt = Array([dom.xul, dom.xres, 0., dom.yul, 0., dom.yres])
    if edges == QX
        nrows = dom.nrows
        ncols = dom.ncols + 1
    elseif edges == QY
        nrows = dom.nrows + 1
        ncols = dom.ncols
    else
        nrows = dom.nrows
        ncols = dom.ncols
    end
    driver = GDAL.getdriverbyname("GTiff")
    ds = GDAL.create(driver, filename, ncols, nrows, 1, GDAL.GDT_Float32, C_NULL)
    GDAL.setgeotransform(ds, gt)
    band = GDAL.getrasterband(ds)
    GDAL.rasterio(band, GDAL.GF_Write, 0, 0, ncols, nrows, data, ncols, nrows, GDAL.GDT_Float32, 0, 0)
    GDAL.close(ds)
end

# boundary conditions
abstract type Boundary end

type Point <: Boundary end

type West <: Boundary end

type East <: Boundary end

type North <: Boundary end

type South <: Boundary end

function get_boundary_type(str::String)
    bctype = Dict("P" => Point, "W" => West, "E" => East, "N" => North, "S" => South)
    return bctype[str.split()[1]]
end

function set_boundary!(bci::Array{BoundaryCondition}, bctype::Point, defn::String)
end

function read_bci(filename::String, dom::Domain)
    bci = Array{BoundaryCondition}(dom.nrows * dom.ncols)
    open(filename, "r") do f
        for line in eachline(f)
            bctype = get_boundary_type(line)
            set_boundary!(bci, bctype, line)
        end
    end
    return bci
end

function read_bdy(filename::String)
    bdy = Dict{String, Array{Float32, 2}}()
    name = ""
    t = 0
    open(filename, "r") do f
        readline(f)
        for line in eachline(f)
            if length(line) > 0
                tokens = split(line)
                if length(tokens) < 2
                    name = tokens[1]
                else
                    if isa(parse(tokens[2]), Number)
                        bdy[name][t, :] = [parse(Float32, s) for s in tokens]
                        t += 1
                    else
                        duration = parse(Int, tokens[1])
                        bdy[name] = Array{Float32, 2}(duration, 2)
                        t = 1
                    end
                end
            end
        end
    end
    return bdy
end

# main function
function run(paramfile::String)
    params = read_params(paramfile)
    z = read_raster(params.dem_file)
    domain = read_domain(params.dem_file)
    n = params.fpfric
    dx = domain.xres
    dy = -domain.yres
end

end
