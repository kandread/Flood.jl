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

type NONE <: BoundaryCondition end

struct HFIX <: BoundaryCondition
    values :: Array{Float32, 1}
end

function hfix(defn::String)
    values = [parse(Float32, split(defn)[5])]
    return HFIX(values)
end

const QFIX = HFIX

const qfix = hfix

struct FREE <: BoundaryCondition
end

struct QVAR <: BoundaryCondition
    name :: String
    prev_time :: Int
end

function qvar(defn::String)
    return QVAR(split(defn)[5], 0)
end

const HVAR = QVAR

const hvar = qvar

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
    return bctype[split(str)[1]]
end

function set_boundary!(bci::Array{BoundaryCondition}, bctype::Point, defn::String, dom::Domain)
    tokens = split(defn)
    x = parse(Float32, tokens[2])
    y = parse(Float32, tokens[3])
    BC = eval(parse(lowercase(tokens[4])))
    j = trunc(Int, (x - dom.xul) / dom.xres) + 1
    i = trunc(Int, (y - dom.yul) / dom.yres) + 1
    pos = dom.nrows * (j - 1) + i
    bci[pos] = BC(defn)
end

function set_boundary!(bci::Array{BoundaryCondition}, bctype::West, defn::String, dom::Domain)
    tokens = split(defn)
    y1 = parse(Float32, tokens[2])
    y2 = parse(Float32, tokens[3])
    i1 = trunc(Int, (y1 - dom.yul) / dom.yres) + 1
    i2 = trunc(Int, (y2 - dom.yul) / dom.yres) + 1
    for i in i1:i2
        pos = i
        bci[pos] = FREE()
    end
end

function set_boundary!(bci::Array{BoundaryCondition}, bctype::East, defn::String, dom::Domain)
    tokens = split(defn)
    y1 = parse(Float32, tokens[2])
    y2 = parse(Float32, tokens[3])
    i1 = trunc(Int, (y1 - dom.yul) / dom.yres) + 1
    i2 = trunc(Int, (y2 - dom.yul) / dom.yres) + 1
    for i in i1:i2
        pos = dom.nrows * (dom.ncols - 1) + i
        bci[pos] = FREE()
    end
end

function set_boundary!(bci::Array{BoundaryCondition}, bctype::North, defn::String, dom::Domain)
    tokens = split(defn)
    x1 = parse(Float32, tokens[2])
    x2 = parse(Float32, tokens[3])
    j1 = trunc(Int, (x1 - dom.xul) / dom.xres) + 1
    j2 = trunc(Int, (x2 - dom.xul) / dom.xres) + 1
    for j in j1:j2
        pos = dom.nrows * (j - 1) + 1
        bci[pos] = FREE()
    end
end

function set_boundary!(bci::Array{BoundaryCondition}, bctype::South, defn::String, dom::Domain)
    tokens = split(defn)
    x1 = parse(Float32, tokens[2])
    x2 = parse(Float32, tokens[3])
    j1 = trunc(Int, (x1 - dom.xul) / dom.xres) + 1
    j2 = trunc(Int, (x2 - dom.xul) / dom.xres) + 1
    for j in j1:j2
        pos = dom.nrows * (j - 1) + dom.nrows - 1
        bci[pos] = FREE()
    end
end

function read_bci(filename::String, dom::Domain)
    bci = Array{BoundaryCondition}(dom.nrows * dom.ncols)
    bci[:] = NONE
    open(filename, "r") do f
        for line in eachline(f)
            bctype = get_boundary_type(line)
            set_boundary!(bci, bctype(), line)
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
