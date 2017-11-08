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
    xul :: Float
    yul :: Float
    xres :: Float
    yres :: Float
    nrows :: Int
    ncols :: Int
end

# parameter structure
struct Params
    dem_file :: String
    sim_time :: Float
    init_tstep :: Float
    fpfric :: Float
    saveint :: Float
    bci_file :: String
    bdy_file :: String
end


end
