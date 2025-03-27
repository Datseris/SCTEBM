using DrWatson
@quickactivate
using ClimateBase
using Statistics
using CairoMakie; CairoMakie.activate!()
include(srcdir("theme.jl"))
include(srcdir("observations_analysis.jl"))

year = 2014
region = "Namibia"
sampling = "monthly"
filectrc = datadir("CTRC", "CTRC_ZhengGRL_2014_$(region)_$(sampling).nc")
input = @dict region sampling

# Also load C from ERA5 to compare with
name(D) = "$(D)D_$(region)_$(sampling)_$(year).nc"
file2 = datadir("ERA5", name(2))
file3 = datadir("ERA5", name(3))

# Notice that if we do not resolve diurnal cycle,
# then it is crucial to use the daily mean field value for SW
SW = ncread(filectrc, "SW_heat_dmean")
CTRC = -ncread(filectrc, "CTRC") # we want it positive definite
LW = CTRC .- SW
I = ncread(file2, "mtdwswrf")
CFZheng = ncread(filectrc, "CF")

CTRC = CFZheng .* CTRC # here we make the cloud fraction all sky
SW = CFZheng .* SW
LW = CFZheng .* LW
SST = ncread(file2, "sst")

# %%

# we assume that the absorbsion of shortwave depends
# on cloud fraction and insolation as a simple coefficient:
# SW = γ*C*I.
# Let's see:

i = copy(vec(I))
c = copy(vec(CFZheng))
s = copy(vec(SW))
l = copy(vec(LW))
r = copy(vec(CTRC))
t = copy(vec(SST))
nans = sort!(union([findall(isnan, x) for x in (i, c, s, r, skipmissing(t))]...))

for x in (i, c, l, s, r, t); deleteat!(x, nans); end

fig, axs = axesgrid(3, 3;
sharex = true,
    size = (800, 800), ylabels = ["CTRCsw", "CTRClw", "CTRC"], xlabels = ["CF", "I", "SST"], sharey = true
)

for (k, x) in enumerate((c, i, t))
    for (j, y) in enumerate((s, l, r))
        ax = axs[j, k]
        scatter!(ax, x, y; color = (COLORS[1], 0.01), strokecolor = (COLORS[1], 0.1), strokewidth = 1)
        linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y
        o, γ = linreg(x, y)
        xmin, xmax = extrema(x)
        lines!(ax, [xmin, xmax], o .+ γ .* [xmin, xmax]; linestyle = :dash, color = "red")
        text!(ax, 0, 0.9; text = "off = $(round(o; sigdigits = 2)), slope = $(round(γ; sigdigits = 2))",
        offset = (2, 2), space = :relative)
    end
end

display(fig)

wsave(plotsdir("ctbbl", "observations", savename("ZhengGRL_dependence", input, "png")), fig)

# %% Just CTRT vs cloud fraction

fig = Figure()
ax = Axis(fig[1,1])
x = c
y = r
scatter!(ax, x, y; color = (COLORS[1], 0.01), strokecolor = (COLORS[1], 0.1), strokewidth = 1)
linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y
o, γ = linreg(x, y)
xmin, xmax = extrema(x)
lines!(ax, [xmin, xmax], o .+ γ .* [xmin, xmax]; linestyle = :dash, color = "red")
text!(ax, 0, 0.9; text = "offset = $(round(o; sigdigits = 2)), slope = $(round(γ; sigdigits = 2))",
offset = (2, 2), space = :relative)
ax.xlabel = "C from Zheng2021"
ax.ylabel = "all-sky CTRC (normalized by C)"

fig

wsave(plotsdir("ctbbl", "observations", savename("CTRC_vs_C_ZhengGRL", input, "png")), fig)
