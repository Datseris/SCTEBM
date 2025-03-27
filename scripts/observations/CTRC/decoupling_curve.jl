using DrWatson
@quickactivate
using ClimateBase
using Statistics
using Dates
using CairoMakie; CairoMakie.activate!()
using ConceptualClimateModels
include(srcdir("theme.jl"))
include(srcdir("observations_analysis.jl"))

# I probably need daily data for this to be more accurate but oh well,
# it is decent enough of a fit already!

year = 2014
region = "Namibia"
sampling = "monthly"
filectrc = datadir("CTRC", "CTRC_ZhengGRL_2014_$(region)_$(sampling).nc")
CTRC = ncread(filectrc, "CTRC")
CFZheng = ncread(filectrc, "CF")

# we need all sky CTRC so we redefine
CTRC = CTRC .* CFZheng

# Also load C from ERA5
name(D) = "$(D)D_$(region)_$(sampling)_$(year).nc"
file2 = datadir("ERA5", name(2))
file3 = datadir("ERA5", name(3))
LCC = ncread(file2, "lcc")
TCC = ncread(file2, "tcc")

# %% Which cloud fraction to use...?
# Should we use the cloud fraction from ZhengGRL or the ERA5?

fig, axs = axesgrid(2, 3;
    size = (800, 400), xlabels = "CTRC",
    ylabels = ["LCC ERA5", "CF Zheng"], sharey = true, sharex = true,
    titles = ("All", "January", "July")
)

for (i, C) in enumerate((LCC, CFZheng))
    for (j, acc) in enumerate((1:12, 1, 6))
        scatter!(axs[i, j], vec(CTRC[Tim(acc)]), vec(C[Tim(acc)]);
        color = (COLORS[1], 0.25), strokewidth = 1, strokecolor = COLORS[1])
    end
end

figuretitle!(fig, "Which cloud fraction to use?")
# wsave(plotsdir("ctbbl", "observations", "CTRC_ZhengGRL_whichC_$region.png"), fig)

fig

# The CF_Zheng has a much clearer relation to CTRC, and is linear,
# as one would expect. Since cooling would be additively integrated
# over an area, and cloud fraction satisfies the same property.

# %% Decoupling curve
# we need the following ingredients:

LHF = ncread(file2, "mslhf")
CBH = ncread(file2, "cbh")
BLH = ncread(file2, "blh")
CLT = @. clamp((BLH - CBH)/(BLH), 0, 1)
D = ncread(file3, "d") ./ 1e-6 # large scale subsidence; use only positive!
D₊ = D[Pre(At(1000))]
RH_b = dropagg(mean, RH[Pre(Between(1000, 950))], Pre) # near surface average should be most accurate

Tsfc = ncread(file2, "t2m")

𝒟 = @. CLT*LHF/CTRC

idxs = findall(x -> ismissing(x) || x==0, CLT)

fig, ax = density(vec(collect(skipmissing(CLT))); label = "CLT: CBH")
density!(vec(collect(skipmissing(CLT2))); label = "CLT2: LCL")
axislegend(ax)
ax.title = "Which definition for CLT to use...?"
# wsave(plotsdir("ctbbl", "observations", "CTRC_ZhengGRL_whichCLT_$region.png"), fig)

display(fig)


# %%

fig, axs = axesgrid(1, 3;
    size = (1000, 400), xlabels = "𝒟",
    ylabels = ["C"], sharey = true, sharex = true,
    titles = ("All months", "Summer", "Winter")
)

accesses = [
    # X -> timemean(X), # this is problematic as it will make countless NaN values...
    X -> X,
    X -> X[Tim(DimensionalData.Where(t -> month(t) ∈ [1, 2, 12]))],
    X -> X[Tim(DimensionalData.Where(t -> month(t) ∈ 6:8))],
]

𝒟version = "Bretherton1997"
# 𝒟version = "Chung2012"

CLTversion = "CBH"

left = [1.0, 1, 1]
right = [0.2, 0.15, 0.25]
scale = [1.0, 1, 1.]
start = [0.5, 0.5, 0.5]
expscale = 0.6 .* [1, 1, 1] .* 2
expstart = [1, 1, 1] .- 0.5
powerstart = [1, 1, 1] .- 0.5
powerscale = [1.2, 1.7, 1.2]

function fitexp(D, Cmax, Cmin, start, scale)
    D < start && return Cmax
    return (Cmax - Cmin)*exp(-scale*(D-start)) + Cmin
end

function fitpower(D, Cmax, Cmin, start, scale)
    D < start && return Cmax
    return 1/(D-start+1)^scale
end

for (j, acc) in enumerate(accesses)
    # access data in appropriate timeframe
    clt = acc(CLT)
    ctrc = acc(CTRC)
    lhf = acc(LHF)
    cf = acc(CFZheng)
    cf2 = acc(LCC)
    d = acc(D₊)

    # the CLT, as estimated by ERA5, has countless 0 entries.
    # We remove them for now by skipping the 0 entries of 𝒟
    # Clean missing/incorrect data
    if 𝒟version == "Bretherton1997"
        𝒹 = @. clt*lhf/ctrc
        idxs = findall(x -> !(ismissing(x) || x==0 || isnan(x)), 𝒹)
        sidxs = findall(>(0), d)
        idxs = intersect(idxs, sidxs)
    elseif 𝒟version == "Chung2012"
        # or, we use an alternative form inspired by chung and teixeira
        𝒹 = @. lhf/ctrc
        idxs = findall(>(0), d)
    end

    # scatter!(axs[j], vec(𝒹[idxs]), vec(cf[idxs]);
    # color = (COLORS[1], 0.05), strokewidth = 0, strokecolor = (COLORS[1], 0.01))
    hexbin!(axs[j], vec(𝒹[idxs]), vec(cf[idxs]); colormap = [:white, :black], bins = (80, 40))

    # color = (COLORS[1], 0.1), strokewidth = 1, strokecolor = COLORS[1])


    # 𝒹 = @. clt2*lhf/ctrc
    # idxs = findall(x -> !(ismissing(x) || x==0), 𝒹)

    # scatter!(axs[j], vec(𝒹[idxs]), vec(cf[idxs]);
    # color = (COLORS[2], 0.25), strokewidth = 1, strokecolor = COLORS[2], marker = MARKERS[2])

    # scatter!(axs[j], 2vec(𝒹[idxs]), vec(cf[idxs]);
    # color = (COLORS[2], 0.25), strokewidth = 1, strokecolor = COLORS[2])

    # plot curve

    dmax = 6
    driver = range(0, dmax; length = 101)
    # case 1:
    reference = start + scale/2
    sigmoid = @. ConceptualClimateModels.sigmoid_expression(driver, left[j], right[j], scale[j], reference[j])
    lines!(axs[j], driver, sigmoid; color = Cycled(1), linestyle = :dash, linewidth = 3)

    # # case 2:
    # fit = @. fitexp(driver, left[j], right[j], expstart[j], expscale[j])
    # lines!(axs[j], driver, fit; color = Cycled(4), linestyle = :dot, linewidth = 3)

    # # case 3:
    # fit = @. fitpower(driver, left[j], right[j], powerstart[j], powerscale[j])
    # lines!(axs[j], driver, fit; color = Cycled(6), linestyle = :dash, linewidth = 3)

    xlims!(axs[j], 0, dmax)
    # maximum(driver))
    # axs[j].xticks = 0:0.5:1.5
end

cb = Colorbar(fig[1, 4]; colormap = [:white, :black], label = "data density")
cb.ticks = ([0, 1], ["low", "high"])
cb.labelpadding = -25
input = @dict region 𝒟version CLTversion

wsave(papersdir("figures", "decoupling_fit.png"), fig)

figuretitle!(fig, savename("all-sky CTRC of Zheng2021", input; connector = ", ", equals = ": "))

display(fig)

# wsave(plotsdir("ctbbl", "observations", savename("CTRC_ZhengGRL", input, "png")), fig)
