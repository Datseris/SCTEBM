# This script loads the continuations (climate change scenario simulations)
# and plots them. The plotting is rather advanced and should be highlighted in Attractors.jl

using DrWatson
@quickactivate
using DynamicalSystems
using Statistics
include(srcdir("theme.jl"))
include("continuation_plotting_helpers.jl")

# Inputs to visualize:
cooling = :q_x
invfix = :difference
Ld = :three_layer # top bottom may be nonsensical; I should study it in the GUI
ΔF = :Gesso2014
co2 = 2
entrain = :Stevens2006
ftrgrad = :none
input = @dict(cooling, invfix, Ld, ΔF, co2, entrain, ftrgrad)

# prepare all columns and rows of the multiplot
ΔFs = [:ctrc, :three_layer, :Gesso2014]
all_input = [(x = copy(input); x[:ΔF] = df; x) for df in ΔFs]
all_used = sort.([[:CO2], [:D], [:D, :CO2], [:D, :CO2, :U]])
titles = map(used -> join(string.(used), ", "), all_used)

# This is how you plot a particular combination:
# fig = Figure(size = (figwidth/4, figheight))
# load_n_plot_continuation!(fig, all_input[1], all_used[2]; imax = 1)
# display(fig)
# wsave(papersdir("figures", "continuation_example_1"), fig)

# and this is how you plot all of them!
fig, axs = axesgrid(length(all_input), length(all_used);
    size = (figwidth, figheight*2.5), titles, xlabels = "time (a. u)",
    ylabels = ["ΔFₛ option = $(i)\nC" for i in 1:3],
    sharex = true, sharey = true,
)

for (i, input) in enumerate(all_input)
    for (j, used) in enumerate(all_used)
        axc = axs[i, j]
        load_n_plot_continuation_only_c!(axc, input, used)
    end
end

# add legend
state_names = sort(collect(keys(state_colors)))
elements = [
   [LineElement(color = state_colors[k][][1]),
    MarkerElement(color = state_colors[k][][1], marker = state_markers[k], markersize = 25)]
    for k in state_names
]

Legend(fig[2, 5], elements, state_names)

display(fig)
# wsave(papersdir("figures", "continuations"), fig)
