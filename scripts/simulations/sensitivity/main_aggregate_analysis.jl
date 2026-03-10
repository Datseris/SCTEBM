using DrWatson
@quickactivate
using DataFrames
using OrderedCollections
include(srcdir("theme.jl"))
include(srcdir("aspects_definition.jl"))

# First, collect all outputs into a single dataframe
foldergroup = ["sims", "random_params"]
folder = datadir(foldergroup...)
df = collect_results(folder)

# fix missing entries for options that were not specifiedin some simulations
replace!(df[!, "cloud_rad"], missing => :const)

# Do any filtering we want if necessary.
# Here we filter the inversion fixed condition
filter!(row -> !(
    row.invfix == :temperature
), df)

# # Here we filter out all rows that have ftrgrad == :weak and invfix == :temperature
# # but only when ΔF_s is NOT :three_layer
# filter!(row -> !(
#     row.ftrgrad == :weak && row.invfix == :temperature && row."ΔF_s" ≠ :three_layer
# ), df)

# # we then filter again so that we remove the same combos for ΔF_s == :three_layer
# # and cloud_rad == :const
# filter!(row -> !(
#     row.ftrgrad == :weak && row.invfix == :temperature && row."ΔF_s" == :three_layer && row.cloud_rad == :const
# ), df)

# and finally sort dataframe by the most important aspects;
# the later in the order the most prioritized the sorting is
# (the dictionary maps aspect to name letter to use for sorting)
sorting_aspects = OrderedDict(
    "co2" => 1, "ftrgrad" => 1, "entrain" => 3, "Ld" => 2, "cloud_rad" => 1, "ΔF_s" => 4
)
sorting_indexes = ones(length(sorting_aspects))
for (aspect, sidx) in sorting_aspects
    sort!(df, [aspect]; by = x -> string(x)[sidx])
end

# %%

# a sorted parameter container for parameters to use
parameters = [:U, :D, :RH₊, :δ_Δ₊T, :δ_FTR, :CO2]
pnames = ["U", "D", "RH₊", "δ_Δ₊T", "δ_FTR", "CO₂"]

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
        scatter!(ax, fill(i, N), 1:N; color, marker, markersize = 15)
    end
    Legend(leg, [MarkerElement(color = COLORS[i], marker = MARKERS[i], markersize = 25) for i in 1:3], ["1", "2", "3"], "option")
    return ax
end

# Large plot of aggregate results
xsize = min(size(df, 1)*30 + 40, 1200)
fig, axs = axesgrid(7, 1;
    xlabels = "model variant", ylabels = ["rMI: SST", "rMI: C", "aspect", "SST [K]", "C", "q_b [g/kg]", "dynamics %"],
    sharex = true, size = (figwidth, figheight*4),
    figure_padding = (5, -20, 5, 10),
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

mean_gap = 3

for (i, m) in enumerate((rmi_sst, rmi_c))
    ax = axs[i]
    # add one more column that is the mean of the MI values
    # (and some empty columns in between, it looks nicer than a rectangle)
    mcol = median(m; dims = 1)
    m = vcat(m, fill(NaN, mean_gap-1, length(mcol)), mcol)
    hm = heatmap!(ax, m; colorscale, colorrange = (lower, upper),
        colormap = Reverse(:balance),
    )
    ax.yticks = (1:length(pnames), pnames)
    if i == 2
        Colorbar(fig[1:2, 2], hm, ticks = [0, 0.5, 1, 2, 3], width = 25, label = "rel. mutual information")
    end
end

# model variants
variants_description!(axs[3], fig[3, 2], df, aspects)

# physics plots
states = ("cumulus%", "multistable%", "stratocumulus%", "limitcycle%", "diverged%", )
S = length(states)

for (i,row) in enumerate(eachrow(df))
    for (j, quantiles) in enumerate((row.quantiles_SST, row.quantiles_C, row.quantiles_q_b))
        scatter!(axs[3+j], fill(i, 3), quantiles; strokecolor = COLORS[4:6], strokewidth = 1.5, color = COLORS[4:6], alpha = 0.25, marker = MARKERS[4:6], markersize = 15)
    end
    # 3rd axis gets dynamical information. I won't plot the diverged in the paper,
    # but we got a surprisingly high number of up to 40% depending on variant.
    dynamical = [getindex(row, idx) for idx in states]
    scatter!(axs[7], fill(i, S), dynamical; strokewidth = 1.5, strokecolor = COLORS[(5-S+1):5], color = COLORS[(5-S+1):5], alpha = 0.25, marker = MARKERS[((6-S+1):6) .- 0], markersize = 15)
end

# also plot the mean:
for (j, quantiles) in enumerate((:quantiles_SST, :quantiles_C, :quantiles_q_b))
    q = reduce(hcat, df[!, quantiles])
    med = vec(mapslices(mean, q; dims = 2))
    scatter!(axs[3+j], fill(L+mean_gap, 3), med; strokecolor = COLORS[4:6], strokewidth = 2, color = COLORS[4:6], alpha = 1.0, marker = MARKERS[4:6], markersize = 20)
end
# and the medians of the dynamical states
dynamical = [mean(df[!, idx]) for idx in states]
scatter!(axs[7], fill(L+mean_gap, S), dynamical; strokewidth = 1.5, strokecolor = COLORS[(5-S+1):5], color = COLORS[(5-S+1):5], alpha = 0.75, marker = MARKERS[((6-S+1):6) .- 0], markersize = 20)

# Add legends
Legend(fig[4:6, 2], [MarkerElement(color = COLORS[i+3], marker = MARKERS[i+3], markersize = 25) for i in 1:3], ["10%", "50%", "90%"], "quantile")
Legend(fig[7, 2], [MarkerElement(color = (COLORS[i+(6-S) .- 1], 1.0), marker = MARKERS[i+(6-S) .- 0], markersize = 25) for i in 1:5], ["Cu", "MS", "Sc", "LC", "F"], "state")

# and final layout stuff:
ticks = 1:20:size(df, 1)
ticks = vcat(ticks, size(df, 1)+mean_gap)
for ax in axs
    ax.xticks = (ticks, string.([ticks[1:end-1]... , "m"]))
end
for ax in axs[1:3]; ax.ygridvisible = false; end

# add vertical lines to separate the most important aspect, the ΔF_s
x = df[!, "ΔF_s"]
idx = unique(z -> x[z], 1:length(x))[2:end] .- 0.5
for ax in axs
    vlines!(ax, idx; color = "black", linewidth = 3, linestyle = :dash)
end

label_axes!(axs; halign = 1.015)
ylims!(axs[4], 282, 308)
ylims!(axs[5], 0.15, 1.05)
xlims!(axs[end], 0.5, size(df, 1)+mean_gap + 1)
axs[3].ylabelpadding = -15
align_labels!(axs[4:7], :y)
display(fig)

# wsave(plotsdir("sims", "aggregate", foldergroup[end]), fig)
wsave(papersdir("figures", "aggregate"), fig; px_per_unit = 4)
