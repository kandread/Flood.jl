module Flood

# constants
const g = 9.81
const alpha = 0.7
const depth_thresh = 0.001
const max_hflow = 10.0

# boundary conditions
abstract type BoundaryCondition end

struct HFIX <: BoundaryCondition
    values :: Array{Float32, 1}
    HFIX(defn::String) = new([parse(Float32, split(defn)[5])])
end

struct QFIX <: BoundaryCondition
    values :: Array{Float32, 1}
    QFIX(defn::String) = new([parse(Float32, split(defn)[5])])
end

struct FREE <: BoundaryCondition
    value :: Float32
    FREE(val=0.0) = new(val)
end

mutable struct QVAR <: BoundaryCondition
    name :: String
    prev_time :: Int
    time :: Array{Float32}
    values :: Array{Float32}
    QVAR(defn::String) = new(split(defn)[5], 1, [], [])
end

mutable struct HVAR <: BoundaryCondition
    name :: String
    prev_time :: Int
    time :: Array{Float32}
    values :: Array{Float32}
    HVAR(defn::String) = new(split(defn)[5], 1, [], [])
end

function interpolate_value(bc::Union{QFIX, HFIX}, t::Float64)
    return bc.values[1]
end

function interpolate_value(bc::Union{QVAR, HVAR}, t::Float64)
    if t <= bc.time[1]
        val = bc.values[1]
    elseif t >= bc.time[end]
        val = bc.values[end]
    else
        ti = bc.prev_time
        while bc.time[ti] <= t
            ti += 1
        end
        dt = bc.time[ti] - bc.time[ti-1]
        a = (t - bc.time[ti-1]) / dt
        val = (1 - a) * bc.values[ti-1] + a * bc.values[ti]
    end
    return val
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
    f = open(filename)
    ncols = parse(Int, split(readline(f))[2])
    nrows = parse(Int, split(readline(f))[2])
    xll = parse(Float32, split(readline(f))[2])
    yll = parse(Float32, split(readline(f))[2])
    cellsize = parse(Float32, split(readline(f))[2])
    nodata = parse(Float32, split(readline(f))[2])
    dom = Domain(xll, yll+nrows*cellsize, cellsize, -cellsize, nrows, ncols)
    close(f)
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
    f = open(filename)
    ncols = parse(Int, split(readline(f))[2])
    nrows = parse(Int, split(readline(f))[2])
    xll = parse(Float32, split(readline(f))[2])
    yll = parse(Float32, split(readline(f))[2])
    cellsize = parse(Float32, split(readline(f))[2])
    nodata = parse(Float32, split(readline(f))[2])
    data = zeros(nrows, ncols)
    for i in 1:nrows
        line = readline(f)
        data[i, :] = [parse(Float32, s) for s in split(line)]
    end
    close(f)
    return data
end

function write_raster(filename::String, data::Array{Float32, 2}, dom::Domain)
    nrows, ncols = size(data)
    f = open(filename, "w")
    @printf(f, "ncols\t%d\n", dom.ncols)
    @printf(f, "nrows\t%d\n", dom.nrows)
    @printf(f, "xllcorner\t%.4f\n", dom.xul)
    yll = dom.yul + dom.nrows * dom.yres
    @printf(f, "yllcorner\t%.4f\n", yll)
    @printf(f, "cellsize\t%.4f\n", dom.xres)
    write(f, "NODATA_value\t-9999\n")
    for i in 1:nrows
        for j in 1:ncols
            @printf(f, "%.3f ", data[i, j])
        end
        write(f, "\n")
    end
    close(f)
end

# boundary conditions
abstract type Boundary end

type Point <: Boundary end

type West <: Boundary end

type East <: Boundary end

type North <: Boundary end

type South <: Boundary end

const Index = Tuple{Int, Int}

function get_boundary_type(str::String)
    bctype = Dict("P" => Point, "W" => West, "E" => East, "N" => North, "S" => South)
    return bctype[split(str)[1]]
end

function set_boundary!(bci::Dict{Index, BoundaryCondition}, bctype::Point, defn::String, dom::Domain)
    tokens = split(defn)
    x = parse(Float32, tokens[2])
    y = parse(Float32, tokens[3])
    BC = eval(parse(tokens[4]))
    j = trunc(Int, (x - dom.xul) / dom.xres) + 1
    i = trunc(Int, (y - dom.yul) / dom.yres) + 1
    bci[(i, j)] = BC(defn)
end

function set_boundary!(bci::Dict{Index, BoundaryCondition}, bctype::West, defn::String, dom::Domain)
    tokens = split(defn)
    y1 = parse(Float32, tokens[2])
    y2 = parse(Float32, tokens[3])
    i1 = trunc(Int, (y1 - dom.yul) / dom.yres) + 1
    i2 = trunc(Int, (y2 - dom.yul) / dom.yres) + 1
    pos = min(i1, i2)
    dist = 0.0
    while pos <= dom.nrows && dist < abs(dom.yres)
        bci[(pos, 1)] = FREE()
        pos +=1
        dist += abs(dom.yres)
    end
end

function set_boundary!(bci::Dict{Index, BoundaryCondition}, bctype::East, defn::String, dom::Domain)
    tokens = split(defn)
    y1 = parse(Float32, tokens[2])
    y2 = parse(Float32, tokens[3])
    i1 = trunc(Int, (y1 - dom.yul) / dom.yres) + 1
    i2 = trunc(Int, (y2 - dom.yul) / dom.yres) + 1
    pos = min(i1, i2)
    dist = 0.0
    while pos <= dom.nrows && dist < abs(dom.yres)
        bci[(pos, dom.ncols)] = FREE()
        pos +=1
        dist += abs(dom.yres)
    end
end

function set_boundary!(bci::Dict{Index, BoundaryCondition}, bctype::North, defn::String, dom::Domain)
    tokens = split(defn)
    x1 = parse(Float32, tokens[2])
    x2 = parse(Float32, tokens[3])
    j1 = trunc(Int, (x1 - dom.xul) / dom.xres) + 1
    j2 = trunc(Int, (x2 - dom.xul) / dom.xres) + 1
    pos = min(j1, j2)
    dist = 0.0
    while pos <= dom.ncols && dist < abs(dom.xres)
        bci[(1, pos)] = FREE()
        pos += 1
        dist += abs(dom.xres)
    end
end

function set_boundary!(bci::Dict{Index, BoundaryCondition}, bctype::South, defn::String, dom::Domain)
    tokens = split(defn)
    x1 = parse(Float32, tokens[2])
    x2 = parse(Float32, tokens[3])
    j1 = trunc(Int, (x1 - dom.xul) / dom.xres) + 1
    j2 = trunc(Int, (x2 - dom.xul) / dom.xres) + 1
    pos = min(j1, j2)
    dist = 0.0
    while pos <= dom.ncols && dist < abs(dom.xres)
        bci[(dom.nrows, pos)] = FREE()
        pos += 1
        dist += abs(dom.xres)
    end
end

function read_bci(filename::String, dom::Domain)
    bci = Dict{Index, BoundaryCondition}()
    open(filename, "r") do f
        for line in eachline(f)
            bctype = get_boundary_type(line)
            set_boundary!(bci, bctype(), line, dom)
        end
    end
    return bci
end

function read_bdy!(bci::Dict{Index, BoundaryCondition}, filename::String)
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
function dry_check!(h::Array{Float32, 2}, Qx::Array{Float32, 2}, Qy::Array{Float32, 2}, dom::Domain, dt::Float32)
    dA = dom.xres * (-dom.yres)
    for i in 1:dom.nrows
        for j in 1:dom.ncols
            dh = dt * (Qx[i, j] - Qx[i, j+1] + Qy[i, j] - Qy[i+1, j]) / dA
            if (dh + h[i, j]) < 0.0
                q_adjust = -(h[i, j] / dh)
                Qx[i, j] *= q_adjust
                Qx[i, j+1] *= q_adjust
                Qy[i, j] *= q_adjust
                Qy[i+1, j] *= q_adjust
            end
        end
    end
end

function calc_sx!(Sx::Array{Float32, 2}, h::Array{Float32, 2}, z::Array{Float32, 2}, dom::Domain, bci::Dict{Index, BoundaryCondition})
    dx = dom.xres
    # inner domain
    for i in 1:dom.nrows
        for j in 2:dom.ncols
            Sx[i, j] = -(h[i, j-1] + z[i, j-1] - h[i, j] -z[i, j]) / dx
        end
    end
    # west and east edges
    for i in 1:dom.nrows
        if isa(get(bci, (i, 1), 0), FREE)
            Sx[i, 1] = bci[(i, 1)].value == 0.0 ? Sx[i, 2] : bci[(i, 1)].value
        else
            Sx[i, 1] = 0.0
        end
        if isa(get(bci, (i, dom.ncols), 0), FREE)
            Sx[i, dom.ncols+1] = bci[(i, dom.ncols)].value == 0.0 ? Sx[i, dom.ncols] : bci[(i, dom.ncols)].value
        else
            Sx[i, dom.ncols+1] = 0.0
        end
    end
end

function calc_sy!(Sy::Array{Float32, 2}, h::Array{Float32, 2}, z::Array{Float32, 2}, dom::Domain, bci::Dict{Index, BoundaryCondition})
    dy = -dom.yres
    # inner domain
    for i in 2:dom.nrows
        for j in 1:dom.ncols
            Sy[i, j] = -(h[i-1, j] + z[i-1, j] - h[i, j] - z[i, j]) / dy
        end
    end
    # north and south edges
    for j in 1:dom.ncols
        if isa(get(bci, (1, j), 0), FREE)
            Sy[1, j] = bci[(1, j)].value == 0.0 ? Sy[2, j] : bci[(1, j)].value
        else
            Sy[1, j] = 0.0
        end
        if isa(get(bci, (dom.nrows, j), 0), FREE)
            Sy[dom.nrows+1, j] = bci[(dom.nrows, j)].value == 0.0 ? Sy[dom.nrows, j] : bci[(dom.nrows, j)].value
        else
            Sy[dom.nrows+1, j] = 0.0
        end
    end
end


function calc_qx!(Qx::Array{Float32, 2}, Sx::Array{Float32, 2}, h::Array{Float32, 2}, z::Array{Float32, 2}, dom::Domain, bci::Dict{Index, BoundaryCondition}, dt::Float32, n::Float32)
    dx = dom.xres
    for i in 1:dom.nrows
        if isa(get(bci, (i, 1), 0), FREE)
            hflow = min(h[i, 1], max_hflow)
            q = Qx[i, 1] / dx
            Qx[i, 1] = (hflow > 0) ? -((abs(q) + abs(g * hflow * dt * Sx[i, 1])) / (1 + g * dt * n * n * abs(q) / hflow^(7/3))) * dx : 0.0
        else
            Qx[i, 1] = 0.0
        end
        if isa(get(bci, (i, dom.ncols), 0), FREE)
            hflow = min(h[i, dom.ncols], max_hflow)
            q = Qx[i, dom.ncols+1] / dx
            Qx[i, dom.ncols+1] = (hflow > 0) ? ((abs(q) + abs(g * hflow * dt * Sx[i, dom.ncols+1])) / (1 + g * dt * n * n * abs(q) / hflow^(7/3))) * dx : 0.0
        else
            Qx[i, dom.ncols+1] = 0.0
        end
    end
    for i in 1:dom.nrows
        for j in 2:dom.ncols
            q = Qx[i, j] / dx
            hflow = max(h[i, j] + z[i, j], h[i, j-1] + z[i, j-1]) - max(z[i, j], z[i, j-1])
            hflow = min(hflow, max_hflow)
            if hflow > 0
                Qx[i, j] = ((q - g * hflow * dt * Sx[i, j]) / (1.0 + g * dt * n * n * abs(q) / hflow^(7/3))) * dx
            else
                Qx[i, j] = 0.0
            end
        end
    end
end

function calc_qy!(Qy::Array{Float32, 2}, Sy::Array{Float32, 2}, h::Array{Float32, 2}, z::Array{Float32, 2}, dom::Domain, bci::Dict{Index, BoundaryCondition}, dt::Float32, n::Float32)
    dy = -dom.yres
    for j in 1:dom.ncols
        if isa(get(bci, (1, j), 0), FREE)
            hflow = min(h[1, j], max_hflow)
            q = Qy[1, j] / dy
            Qy[1, j] = (hflow > 0) ? -((abs(q) + abs(g * hflow * dt * Sy[1, j])) / (1 + g * dt * n * n * abs(q) / hflow^(7/3))) * dy : 0.0
        else
            Qy[1, j] = 0.0
        end
        if isa(get(bci, (dom.nrows, j), 0), FREE)
            hflow = min(h[dom.nrows, j], max_hflow)
            q = Qy[dom.nrows+1, j] / dy
            Qy[dom.nrows+1, j] = (hflow > 0) ? ((abs(q) + abs(g * hflow * dt * Sy[dom.nrows+1, j])) / (1 + g * dt * n * n * abs(q) / hflow^(7/3))) * dy : 0.0
        else
            Qy[dom.nrows+1, j] = 0.0
        end
    end
    for i in 2:dom.nrows
        for j in 1:dom.ncols
            q = Qy[i, j] / dy
            hflow = max(h[i, j] + z[i, j], h[i-1, j] + z[i-1, j]) - max(z[i, j], z[i-1, j])
            hflow = min(hflow, max_hflow)
            if hflow > 0
                Qy[i, j] = ((q - g * hflow * dt * Sy[i, j]) / (1.0 + g * dt * n * n * abs(q) / hflow^(7/3))) * dy
            else
                Qy[i, j] = 0.0
            end
        end
    end
end

function calc_h!(h::Array{Float32, 2}, Qx::Array{Float32, 2}, Qy::Array{Float32, 2}, dom::Domain, bci::Dict{Index, BoundaryCondition}, t::Real, dt::Real)
    dx = dom.xres
    dy = -dom.yres
    for i in 1:dom.nrows
        for j in 1:dom.ncols
            if isa(get(bci, (i, j), 0), HFIX) || isa(get(bci, (i, j), 0), HVAR)
                h[i, j] = interpolate_value(bci[(i, j)], t)
            else
                h[i, j] += (dt * (Qx[i, j] - Qx[i, j+1] + Qy[i, j] - Qy[i+1, j]) / (dx * dy))
            end
            if isa(get(bci, (i, j), 0), QFIX) || isa(get(bci, (i, j), 0), QVAR)
                h[i, j] += interpolate_value(bci[(i, j)], t) * dt  / dx
            end
            if h[i, j] < depth_thresh
                h[i, j] = 0.0
            end
        end
    end
end

"""Calculate time step for numerically stable solution based on Courant–Friedrichs–Lewy condition."""
function calc_timestep(h::Array{Float32, 2}, dom::Domain, params::Params)
    dx = dom.xres
    hmax = maximum(h)
    if hmax > depth_thresh
        dt = alpha * dx / sqrt(hmax * g)
    else
        dt = params.init_tstep
    end
    return dt
end


# main function
function run(paramfile::String)
    params = read_params(paramfile)
    dom = read_domain(params.dem_file)
    z = read_raster(params.dem_file)
    bci = read_bci(params.bci_file, dom)
    if !isempty(params.bdy_file)
        read_bdy!(bci, params.bdy_file)
    end
    n = params.fpfric
    h = zeros(Float32, dom.nrows, dom.ncols)
    Qx = zeros(Float32, dom.nrows, (dom.ncols+1))
    Qy = zeros(Float32, (dom.nrows+1), dom.ncols)
    Sx = zeros(Float32, dom.nrows, (dom.ncols+1))
    Sy = zeros(Float32, (dom.nrows+1), dom.ncols)
    cur_save = 0.0
    t = 0.0
    while t < params.sim_time
        dt = params.init_tstep
        calc_sx!(Sx, h, z, dom, bci)
        calc_sy!(Sy, h, z, dom, bci)
        calc_qx!(Qx, Sx, h, z, dom, bci, dt, n)
        calc_qy!(Qy, Sy, h, z, dom, bci, dt, n)
        dry_check!(h, Qx, Qy, dom, dt)
        calc_h!(h, Qx, Qy, dom, bci, t, dt)
        t += dt
        if t >= cur_save
            write_raster(@sprintf("h-%04d.asc", cur_save / params.saveint), h, dom)
            write_raster(@sprintf("qx-%04d.asc", cur_save / params.saveint), Qx, dom)
            write_raster(@sprintf("qy-%04d.asc", cur_save / params.saveint), Qy, dom)
            cur_save += params.saveint
        end
    end
end

end
