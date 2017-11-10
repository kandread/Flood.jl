module Flood

import GDAL

# constants
const g = 9.81
const alpha = 0.7
const depth_thresh = 0.001
const max_hflow = 10.0

# grid types
@enum Grid H QX QY

# boundary conditions
abstract type BoundaryCondition end

struct HFIX <: BoundaryCondition
    values :: Array{Float32, 1}
end

function hfix(defn::String)
    values = [parse(Float32, split(defn)[5])]
    return HFIX(values)
end

struct QFIX <: BoundaryCondition
    values :: Array{Float32, 1}
end

function qfix(defn::String)
    values = [parse(Float32, split(defn)[5])]
    return QFIX(values)
end

struct FREE <: BoundaryCondition
    value :: Float32
end

function free(val=0.0)
    return FREE(val)
end

mutable struct QVAR <: BoundaryCondition
    name :: String
    prev_time :: Int
    time :: Array{Float32}
    values :: Array{Float32}
end

function qvar(defn::String)
    return QVAR(split(defn)[5], 0, [], [])
end

mutable struct HVAR <: BoundaryCondition
    name :: String
    prev_time :: Int
    time :: Array{Float32}
    values :: Array{Float32}
end

function hvar(defn::String)
    return HVAR(split(defn)[5], 0, [], [])
end

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
    bdy_file :: String
    output_dir :: String

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
    bci_file = haskey(data, "bcifile") ? data["bcifile"] : ""
    bdy_file = haskey(data, "bdyfile") ? data["bdyfile"] : ""
    saveint = haskey(data, "saveint") ? parse(Float32, data["saveint"]) : nothing
    outdir = haskey(data, "dirroot") ? data["dirroot"] : "."
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

function write_raster(filename::String, data::Array{Float32}, dom::Domain, edges::Grid, format="GTiff")
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
    driver = GDAL.getdriverbyname(format)
    ds = GDAL.create(driver, filename, ncols, nrows, 1, GDAL.GDT_Float32, C_NULL)
    GDAL.setgeotransform(ds, gt)
    band = GDAL.getrasterband(ds, 1)
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

function set_boundary!(bci::Dict{Int, BoundaryCondition}, bctype::Point, defn::String, dom::Domain)
    tokens = split(defn)
    x = parse(Float32, tokens[2])
    y = parse(Float32, tokens[3])
    BC = eval(parse(lowercase(tokens[4])))
    j = trunc(Int, (x - dom.xul) / dom.xres) + 1
    i = trunc(Int, (y - dom.yul) / dom.yres) + 1
    pos = dom.nrows * (j - 1) + i
    bci[pos] = BC(defn)
end

function set_boundary!(bci::Dict{Int, BoundaryCondition}, bctype::West, defn::String, dom::Domain)
    tokens = split(defn)
    y1 = parse(Float32, tokens[2])
    y2 = parse(Float32, tokens[3])
    i1 = trunc(Int, (y1 - dom.yul) / dom.yres) + 1
    i2 = trunc(Int, (y2 - dom.yul) / dom.yres) + 1
    pos = min(i1, i2)
    dist = 0.0
    while pos <= dom.nrows && dist < abs(dom.yres)
        bci[pos] = free()
        pos +=1
        dist += abs(dom.yres)
    end
end

function set_boundary!(bci::Dict{Int, BoundaryCondition}, bctype::East, defn::String, dom::Domain)
    tokens = split(defn)
    y1 = parse(Float32, tokens[2])
    y2 = parse(Float32, tokens[3])
    i1 = trunc(Int, (y1 - dom.yul) / dom.yres) + 1
    i2 = trunc(Int, (y2 - dom.yul) / dom.yres) + 1
    pos = min(i1, i2)
    dist = 0.0
    while pos <= dom.nrows && dist < abs(dom.yres)
        bci[pos+(dom.ncols-1)*dom.nrows] = free()
        pos +=1
        dist += abs(dom.yres)
    end
end

function set_boundary!(bci::Dict{Int, BoundaryCondition}, bctype::North, defn::String, dom::Domain)
    tokens = split(defn)
    x1 = parse(Float32, tokens[2])
    x2 = parse(Float32, tokens[3])
    j1 = trunc(Int, (x1 - dom.xul) / dom.xres) + 1
    j2 = trunc(Int, (x2 - dom.xul) / dom.xres) + 1
    pos = min(j1, j2)
    dist = 0.0
    while pos <= dom.ncols && dist < abs(dom.xres)
        bci[1+(pos-1)*dom.nrows] = free()
        pos += 1
        dist += abs(dom.xres)
    end
end

function set_boundary!(bci::Dict{Int, BoundaryCondition}, bctype::South, defn::String, dom::Domain)
    tokens = split(defn)
    x1 = parse(Float32, tokens[2])
    x2 = parse(Float32, tokens[3])
    j1 = trunc(Int, (x1 - dom.xul) / dom.xres) + 1
    j2 = trunc(Int, (x2 - dom.xul) / dom.xres) + 1
    pos = min(j1, j2)
    dist = 0.0
    while pos <= dom.ncols && dist < abs(dom.xres)
        bci[dom.nrows*pos] = free()
        pos += 1
        dist += abs(dom.xres)
    end
end

function read_bci(filename::String, dom::Domain)
    bci = Dict{Int, BoundaryCondition}()
    open(filename, "r") do f
        for line in eachline(f)
            bctype = get_boundary_type(line)
            set_boundary!(bci, bctype(), line, dom)
        end
    end
    return bci
end

function read_bdy!(bci::Dict{Int64, BoundaryCondition}, filename::String)
    name = ""
    pos = 0
    t = 0
    open(filename, "r") do f
        readline(f)
        for line in eachline(f)
            if length(line) > 0
                tokens = split(line)
                if length(tokens) < 2
                    name = tokens[1]
                    pos = [k for (k, v) in bci if isa(v, Union{QVAR, HVAR}) && v.name == name][1]
                else
                    if isa(parse(tokens[2]), Number)
                        bci[pos].time[t] = parse(Float32, tokens[2])
                        bci[pos].values[t] = parse(Float32, tokens[1])
                        t += 1
                    else
                        duration = parse(Int, tokens[1])
                        bci[pos].time = zeros(duration)
                        bci[pos].values = zeros(duration)
                        t = 1
                    end
                end
            end
        end
    end
end

# floodplain flow
function dry_check!(h::Array{Float32}, Qx::Array{Float32}, Qy::Array{Float32}, dom::Domain, dt::Float32)
    dA = dom.xres * (-dom.yres)
    for i in 1:dom.nrows
        for j in 1:dom.ncols
            dh = dt * (Qx[i+(j-1)*dom.nrows] - Qx[i+j*dom.nrows] + Qy[i+(j-1)*(dom.nrows+1)] - Qy[i+1+(j-1)*(dom.nrows+1)]) / dA
            if (dh + h[i+(j-1)*dom.nrows]) < 0.0
                q_adjust = -(h[i+(j-1)*dom.nrows] / dh)
                Qx[i+(j-1)*dom.nrows] *= q_adjust
                Qx[i+j*dom.nrows] *= q_adjust
                Qy[i+(j-1)*(dom.nrows+1)] *= q_adjust
                Qy[i+1+(j-1)*(dom.nrows+1)] *= q_adjust
            end
        end
    end
end

function calc_sx!(Sx::Array{Float32}, h::Array{Float32}, z::Array{Float32}, dom::Domain, bci::Dict{Int64, BoundaryCondition})
    dx = dom.xres
    # inner domain
    for i in 1:dom.nrows
        for j in 2:dom.ncols
            Sx[i+(j-1)*dom.nrows] = -(h[i+(j-2)*dom.nrows] + z[i+(j-2)*dom.nrows] - h[i+(j-1)*dom.nrows] - z[i+(j-1)*dom.nrows]) / dx
        end
    end
    # west and east edges
    for i in 1:dom.nrows
        if isa(get(bci, i, 0), FREE)
            Sx[i] = (bci[i].value == 0.0) ? Sx[i+dom.nrows] : bci[i].value
        end
        if isa(get(bci, i+dom.ncols*dom.nrows, 0), FREE)
            Sx[i+dom.ncols*dom.nrows] = (bci[i+dom.ncols*dom.nrows].value == 0.0) ? Sx[i+(dom.ncols-1)*dom.nrows] : bci[i+dom.ncols*dom.nrows].value
        end
    end
end

function calc_sy!(Sy::Array{Float32}, h::Array{Float32}, z::Array{Float32}, dom::Domain, bci::Dict{Int64, BoundaryCondition})
    dy = -dom.yres
    # inner domain
    for i in 2:dom.nrows
        for j in 1:dom.ncols
            Sy[i+(j-1)*(dom.nrows+1)] = -(h[i-1+(j-1)*dom.nrows] + z[i-1+(j-1)*dom.nrows] - h[i+(j-1)*dom.nrows] - z[i+(j-1)*dom.nrows]) / dy
        end
    end
    # north and south edges
    for j in 1:dom.ncols
        if isa(get(bci, 1+(j-1)*dom.nrows, 0), FREE)
            Sy[1+(j-1)*(dom.nrows+1)] = (bci[1+(j-1)*dom.nrows].value == 0.0) ? Sy[2+(j-1)*(dom.nrows+1)] : bci[1+(j-1)*dom.nrows].value
        end
        if isa(get(bci, 1+j*(dom.nrows+1), 0), FREE)
            Sy[1+j*(dom.nrows+1)] = (bci[1+j*dom.nrows].value == 0.0) ? Sy[j*(dom.nrows+1)] : bci[1+j*dom.nrows].value
        end
    end
end


function calc_qx!(Qx::Array{Float32}, Sx::Array{Float32}, h::Array{Float32}, z::Array{Float32}, dom::Domain, dt::Float32, n::Float32)
    dx = dom.xres
    for i in 1:dom.nrows
        for j in 1:dom.ncols+1
            q = Qx[i+(j-1)*dom.nrows] / dx
            if j == 1
                hflow = h[i+(j-1)*dom.nrows]
            elseif j == dom.ncols+1
                hflow = h[i+(j-2)*dom.nrows]
            else
                hflow = max(h[i+(j-1)*dom.nrows] + z[i+(j-1)*dom.nrows], h[i+(j-2)*dom.nrows] + z[i+(j-2)*dom.nrows]) - max(z[i+(j-1)*dom.nrows], z[i+(j-2)*dom.nrows])
            end
            hflow = min(hflow, max_hflow)
            if hflow > 0
                Qx[i+(j-1)*dom.nrows] = ((q - g * hflow * dt * Sx[i+(j-1)*dom.nrows]) / (1.0 + g * dt * n * n * abs(q) / hflow^(7/3))) * dx
            else
                Qx[i+(j-1)*dom.nrows] = 0.0
            end
        end
    end
end

function calc_qy!(Qy::Array{Float32}, Sy::Array{Float32}, h::Array{Float32}, z::Array{Float32}, dom::Domain, dt::Float32, n::Float32)
    dy = -dom.yres
    for i in 1:dom.nrows+1
        for j in 1:dom.ncols
            q = Qy[i+(j-1)*(dom.nrows+1)] / dy
            if i == 1
                hflow = h[1+(j-1)*dom.nrows]
            elseif i == dom.nrows+1
                hflow = h[j*dom.nrows]
            else
                hflow = max(h[i+(j-1)*dom.nrows] + z[i+(j-1)*dom.nrows], h[i-1+(j-1)*dom.nrows] + z[i-1+(j-1)*dom.nrows]) - max(z[i+(j-1)*dom.nrows], z[i-1+(j-1)*dom.nrows])
            end
            hflow = min(hflow, max_hflow)
            if hflow > 0
                Qy[i+(j-1)*(dom.nrows+1)] = ((q - g * hflow * dt * Sy[i+(j-1)*(dom.nrows+1)]) / (1.0 + g * dt * n * n * abs(q) / hflow^(7/3))) * dy
            else
                Qy[i+(j-1)*(dom.nrows+1)] = 0.0
            end
        end
    end
end

function calc_h!(h::Array{Float32}, Qx::Array{Float32}, Qy::Array{Float32}, dom::Domain, bci::Dict{Int64, BoundaryCondition}, t::Real, dt::Real)
    dx = dom.xres
    dy = -dom.yres
    for i in 1:dom.nrows
        for j in 1:dom.ncols
            if isa(get(bci, i+(j-1)*dom.nrows, 0), HFIX) || isa(get(bci, i+(j-1)*dom.nrows, 0), HVAR)
                # FIXME: Interpolate time series
                # h[i+(j-1)*dom.nrows] = interpolate_value(bci[i+(j-1)*dom.nrows], t)
                h[i+(j-1)*dom.nrows] = bci[i+(j-1)*dom.nrows].values[1]
            else
                h[i+(j-1)*dom.nrows] += (dt * (Qx[i+(j-1)*dom.nrows] - Qx[i+j*dom.nrows] + Qy[i+(j-1)*(dom.nrows+1)] - Qy[i+1+(j-1)*(dom.nrows+1)]) / (dx * dy))
            end
            if isa(get(bci, i+(j-1)*dom.nrows, 0), QFIX) || isa(get(bci, i+(j-1)*dom.nrows, 0), QVAR)
                # FIXME: Interpolate time series
                # h[i+(j-1)*dom.nrows] += interpolate_value(bci[i+(j-1)*dom.nrows], t) * dt / dx
                h[i+(j-1)*dom.nrows] = h[i+(j-1)*dom.nrows] + (bci[i+(j-1)*dom.nrows].values[1]) * dt / dx
            end
            if h[i+(j-1)*dom.nrows] < depth_thresh
                h[i+(j-1)*dom.nrows] = 0.0
            end
        end
    end
end


# main function
function run(paramfile::String)
    params = read_params(paramfile)
    z = read_raster(params.dem_file)
    z = convert(Array{Float32}, z)
    dom = read_domain(params.dem_file)
    bci = read_bci(params.bci_file, dom)
    if !isempty(params.bdy_file)
        read_bdy!(bci, params.bdy_file)
    end
    n = params.fpfric
    h = zeros(Float32, dom.nrows*dom.ncols)
    Qx = zeros(Float32, dom.nrows*(dom.ncols+1))
    Qy = zeros(Float32, (dom.nrows+1)*dom.ncols)
    Sx = zeros(Float32, dom.nrows*(dom.ncols+1))
    Sy = zeros(Float32, (dom.nrows+1)*dom.ncols)
    cur_save = 0.0
    t = 0.0
    while t < params.sim_time
        dt = params.init_tstep
        calc_sx!(Sx, h, z, dom, bci)
        calc_sy!(Sy, h, z, dom, bci)
        calc_qx!(Qx, Sx, h, z, dom, dt, n)
        calc_qy!(Qy, Sy, h, z, dom, dt, n)
        dry_check!(h, Qx, Qy, dom, dt)
        calc_h!(h, Qx, Qy, dom, bci, t, dt)
        t += dt
        if t >= cur_save
            write_raster(@sprintf("h-%04d.tif", cur_save), h, dom, H)
            cur_save += params.saveint
        end
    end
end

end
