# plot limit cycle timeseries to elucidate how it exists
using DrWatson
@quickactivate
using DynamicalSystems
using ConceptualClimateModels
using Statistics: mean
using ProgressMeter
include(srcdir("sctebm_setups.jl"))
include(srcdir("theme.jl"))

# setup simulation
cooling = :q_x
invfix = :difference
Ld = :three_layer
entrain = :Stevens2006
ΔF = :ctrc
invdec = true # having this as `true` promotes the LC stability up to very high τ_C
input = @dict(cooling, invfix, Ld, ΔF, entrain, invdec)

starting_parameters = Dict(
    :D => 2.5e-6, :U => 7, :q_x_rate => 2, :RH₊ => 0.1,
    :τ_C => 2.0, :τ_SST => 50.0
)

ds, eqs = sctebm_setup(; input..., starting_parameters)

obs = [:C, :SST, :q_b, :𝒟, :CTRC, :RCT, :LHF]

X, t = trajectory(ds, 100.0; Δt = 0.01, save_idxs = obs, Ttr = 200.0)

# estimate period
x = X[:, 1]
t0 = t[1]
T_i = round(Int, estimate_period(X[:, 1], :zerocrossing))
T = t[T_i] - t0

fig = Figure(size = (figwidth, 2figheight))
gl1 = GridLayout(fig[1,1])

axs = axesgrid!(gl1, length(obs), 1;
sharex = true, xlabels = "time (day)",
ylabels = string.(obs))

tscales = [current_parameter(ds, :τ_C),
1/current_parameter(ds, :D)/CTMLM.sec_in_day,
1200/(current_parameter(ds, :U)*current_parameter(ds, :c_d) + observe_state(ds, :w_e))/CTMLM.sec_in_day,
]
tlabels = ["τ_C", "τ_z", "τ_q"]

for (i, ax) in enumerate(axs)
    lines!(ax, t .- t0, X[:, i]; color = Cycled(i))
    vlines!(ax, [T]; linestyle = :dash, linewidth = 2, color = "black")

    # add timescale lines
    # for (j, τ) in enumerate(tscales)
    #     vlines!(ax, [τ, T+τ]; linestyle = :dash, linewidth = 2,
    #     label = tlabels[j], color = Cycled(j+1))
    # end
    ax.xticks = WilkinsonTicks(10; k_min = 8, k_max = 15)
    ax.yticks = WilkinsonTicks(2; k_min = 2, k_max = 3)
end
xlims!(axs[end], 0, t[end] - t0)
rowgap!(gl1, 0)

# A nice subplot of the limit cycle for different timescales?
# one panel period and the other ...? To show that it does not depend on timescales?

axs = axesgrid!(fig[1,2], 2, 1;
sharex = true, sharey = true, xlabels = L"\tau_C", ylabels = L"\tau_\mathrm{SST}")
colsize!(fig.layout, 2, Relative(0.33))

display(fig)

# %%

empty!.(axs)

τCs = 1:0.25:6.5
τSSTs = 10:5:100
periods = randn(length.((τCs, τSSTs)))
TsDiff = randn(length.((τCs, τSSTs)))

@showprogress for (i, τC) in enumerate(τCs)
    for (j, τSST) in enumerate(τSSTs)
        set_parameter!(ds, :τ_C, τC)
        set_parameter!(ds, :τ_SST, τSST)
        Δt = 0.05
        X, t = trajectory(ds, 2000.0; Δt, save_idxs = [:SST], Ttr = 1000.0)
        x = X[:, 1]
        TsDiff[i, j] = maximum(x) - minimum(x)
        if TsDiff[i, j] < 0.1
            period = NaN
        else
            period = estimate_period(x .- mean(x), :periodogram)*Δt
        end
        periods[i,j] = period
    end
end

hmp = heatmap!(axs[1], τCs, τSSTs, periods; colormap = :dense)
Colorbar(fig[1,2][1, 2], hmp, label = "period (day)")
hmt = heatmap!(axs[2], τCs, τSSTs, TsDiff; colormap = :dense)
Colorbar(fig[1,2][2, 2], hmt, label = "SST max - min")

display(fig)

wsave(papersdir("figures", "limit_cycle"), fig)
