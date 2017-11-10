using Base.Test
using Flood

z = Flood.read_raster("dtm.asc")
@test z ≈ map(Int32, [20  20  10  20  20  20  20  10  20  20  10  10  10  10  10  20  20  10  20  20  20  20  10  20  20][1, :])
