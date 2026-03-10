ENV["COLORSCHEME"] = "CloudySky" # or others, see `plottheme.jl`
ENV["BGCOLOR"] = :white

using MakieForProjects
using CairoMakie

figwidth = 1200
figheight = 300 # the height of a single panel
update_theme!(;
    size = (figwidth, 2figheight),
)

function color_from_u0(u0::Dict; Cg = [0.05, 0.5, 0.95], SSTg = [285.0, 295, 305])
    # D = length(Cg)*length(SSTg)
    C = u0[:C]; SST = u0[:SST]
    idxC = findmin(abs2.(C .- Cg))[2]
    idxS = findmin(abs2.(SST .- SSTg))[2]
    # convert to linear
    l = LinearIndices((length(Cg), length(SSTg)))[idxC, idxS]
    c = to_color(COLORS[l])
    # add some random variability
    r() = 1 + 2e-2randn()
    c = RGBf(c.r*r(), c.g*r(), c.b*r())
end

color_from_u0(::AbstractArray) = to_color(COLORS[rand(Int)])

"convenience function for x-axis of densities of observables."
function density_limits(symbol)
    get(Dict(
        :SST => (285, 305),
        :C => (0, 1.2),
        :q_b => (5, 20),
        :LHF => (25, 225),
        :T₊ => (275, 305),
        :Ld => (330, 450),
        :RCT => (-0.05, 0.6),
        :SHF => (-10, 30),
        :Lnet => (0, 80),
        :z_b => (0, 2000),
        :EIS => (-1, 12),
    ), symbol, (nothing, nothing))
end

# so that we can label axesgrid easier via `transpose`:
import LinearAlgebra
LinearAlgebra.transpose(ax::Axis) = ax
