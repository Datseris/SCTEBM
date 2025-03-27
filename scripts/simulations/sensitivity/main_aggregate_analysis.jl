using DrWatson
@quickactivate
using DataFrames
using CairoMakie
using ClipData # allows copy pasting dataframes into Excel
using OrderedCollections
include(srcdir("theme.jl"))
include(srcdir("simulations_process.jl"))

# First, collect all outputs into a single dataframe
foldergroup = ["sims", "random_params"]
folder = datadir(foldergroup...)
dffull = collect_results(folder)

# Do any filtering we want if necessary
df = dffull
df = filter(row -> row.cooling == :q_x, df)
df = filter(row -> row.entrain == :Gesso2014, df)

# and rename columns that used an old notation
replace!(df[!, "ftrgrad"], :strong => :weak, missing => :weak)

# and finally sort dataframe by the most important aspects
sort!(df, ["ΔF"]; by = x -> string(x)[2], rev = true)

# a sorted parameter container for parameters to use
parameters = [:U, :D, :RH₊, :δ_Δ₊T, :δ_FTR, :CO2]
pnames = ["U", "D", "RH₊", "δ_Δ₊T", "δ_FTR", "CO₂"] # We correct the missing entries as overtime simulations progressed

# and a sorted aspect container for the aspects to use
# that must be in the same sorting as the options structure
aspects = collect(keys(ASPECTS_OPTIONS))

# Colormap of MI values
function rmi_matrix(df, var::Symbol, parameters)
    M = zeros(size(df, 1), length(parameters))
    for i in axes(df, 1)
        rmi = df[i, "rMI"][var]
        for (j, par) in enumerate(parameters)
            M[i, j] = rmi[par]
        end
    end
    return M
end

rmi_sst = rmi_matrix(df, :SST, parameters)
rmi_c = rmi_matrix(df, :C, parameters)

# create a plot from the dataframe
# max 3 unique options per aspect
function variants_description!(ax::Axis, leg, df, aspects)
    N = length(aspects)
    ax.yticks = (1:N, replace(aspects, "ftrgrad" => "tftr"))
    for (i, row) in enumerate(eachrow(df))
        option_idxs = aspects_to_option_idx(row)
        color = [COLORS[o] for o in option_idxs]
        marker = [MARKERS[o] for o in option_idxs]
        scatter!(ax, fill(i, N), 1:N; color, marker, markersize = 20)
    end
    Legend(leg, [MarkerElement(color = COLORS[i], marker = MARKERS[i], markersize = 25) for i in 1:3], ["1", "2", "3"], "option")
    return ax
end

# Large plot of aggregate results
xsize = min(size(df, 1)*40 + 60, 1500)

fig, axs = axesgrid(7, 1;
    xlabels = "model variant", ylabels = ["rMI: SST", "rMI: C", "aspect", "SST", "C", "q_b", "dynamics %"],
    sharey = true, size = (xsize, 1500), sharex = true,
)

# Mutual information plots
using Statistics: mean, median, std
L = size(df, 1)

# Make nonuniform colorscale
upper = 3
lower = 0
midpoint = 1
ratio = (upper - midpoint) / (midpoint - lower)
colorscale = Makie.ReversibleScale(
    x -> x > midpoint ? (x - midpoint) / ratio + midpoint : x,
    x -> x > midpoint ? (x - midpoint) * ratio + midpoint : x
)

for (i, m) in enumerate((rmi_sst, rmi_c))
    ax = axs[i]
    # add one more row that is the mean of the MI values
    mcol = median(m; dims = 1)
    m = vcat(m, mcol)
    hm = heatmap!(ax, m; colorscale, colorrange = (lower, upper),
        colormap = Reverse(:balance),
    )
    ax.yticks = (1:length(pnames), pnames)
    ax.xticks = (1:L+1, string.([1:L..., "m"]))
    # rectangles highlighting the columns
    rect = Rect(L+1 - 0.5, 0.5, 1, length(pnames))
    poly!(ax, rect; color = ("green", 0.0), strokecolor = "white", strokewidth = 5)

    if i == 2
        Colorbar(fig[1:2, 2], hm, ticks = [0, 0.5, 1, 2, 3], width = 25, label = "rel. mutual information")
    end
end

# model variants
variants_description!(axs[3], fig[3, 2], df, aspects)

# physics plots
states = ("limitcycle%", "cumulus%", "stratocumulus%", "diverged%", "multistable%",)
S = length(states)

for (i,row) in enumerate(eachrow(df))
    for (j, quartiles) in enumerate((row.quartiles_SST, row.quartiles_C, row.quartiles_q))
        scatter!(axs[3+j], fill(i, 3), quartiles; color = COLORS[4:6], marker = MARKERS[4:6], markersize = 22)
    end
    # 3rd axis gets dynamical information. I won't plot the diverged in the paper,
    # but we got a surprisingly high number of up to 40% depending on variant.
    dynamical = [getindex(row, idx) for idx in states]
    scatter!(axs[7], fill(i, S), dynamical; strokewidth = 1.5, strokecolor = COLORS[(6-S+1):6], color = COLORS[(6-S+1):6], alpha = 0.25, marker = MARKERS[(6-S+1):6], markersize = 20)
end

Legend(fig[4:6, 2], [MarkerElement(color = COLORS[i+3], marker = MARKERS[i+3], markersize = 25) for i in 1:3], ["25%", "50%", "75%"], "quartile")
Legend(fig[7, 2], [MarkerElement(strokewidth = 2, strokecolor = COLORS[(6-S+1):6], color = (COLORS[i+(6-S)], 0.25), marker = MARKERS[i+(6-S)], markersize = 25) for i in 1:5], ["LC", "Cu", "Sc", "F", "MS"], "state")

xlims!(axs[end], 0.5, size(df, 1)+1.5)

label_axes!(axs; halign = 1.015)

display(fig)

# wsave(plotsdir("sims", "aggregate", foldergroup[end]), fig)
wsave(papersdir("figures", "aggregate"), fig)

# %% Any sorts of subsequent analysis here
function remove_ending_percentage!(df)
    n = names(df)
    problematic = findall(endswith("%"), n)
    renames = [p => p[1:end-1] for p in n[problematic]]
    rename!(df, renames...)
end
remove_ending_percentage!(df)

using DataFramesMeta
sim_options = Symbol.(aspects)

# Find 10 most stable simulations
x = @chain df begin
    DataFramesMeta.@orderby(:diverged)
    DataFramesMeta.@rsubset(:diverged < 5)
    DataFramesMeta.@select($(sim_options...))
end
display(x)

# Okay it becomes clear that the most stable model configuration is with the stevens
# entrainmenta and the difference inversion condition.