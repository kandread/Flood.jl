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


end
