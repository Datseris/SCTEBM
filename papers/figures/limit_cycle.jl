using DrWatson
@quickactivate
using DynamicalSystems
using OrdinaryDiffEq
using ConceptualClimateModels
using CairoMakie
using OrderedCollections: OrderedDict
import ClimateBase

include(srcdir("theme.jl"))
include(srcdir("CloudToppedBulkBoundaryLayer", "CloudToppedBulkBoundaryLayer.jl"))
include(srcdir("ctbbl_setups.jl"))

# This setup, with the parameters set below, results in a limit cycle
w_m = true
ftrgrad = :weak
cfversion = :Singer2023
co2_sst = false
co2_ΔT = true
radiation = :detailed # detailed or merged; for merged we need smaller drag coefficient
q_x = radiation == :detailed
input = @strdict(w_m, q_x, cfversion, co2_sst, radiation, co2_ΔT, ftrgrad)

starting_parameters = Dict(:V => 6.68, :T_FTR_0 => 289.0)

ds, eqs = ctbbl_setup(input; starting_parameters)

# %% plot limit cycle
# Choose the observables
obs = [:C, :SST, :qᵦ, :zᵦ,  :Λ, :CTRC, :CLT, :LHF]

X, t = trajectory(ds, 100.0; Δt = 0.01, save_idxs = obs, Ttr = 200.0)
x = X[:, 1]

tf_i = findfirst(j -> (x[j] > x[j+1]) && (x[j] > x[j-1]), 2:length(x)-1)
tf = t[tf_i]
T_i = round(Int, estimate_period(X[:, 1], :zerocrossing))
T = t[T_i] - t[1]

fig, axs = axesgrid(length(obs)÷2, 2;
size = (1200, 400), sharex = true, xlabels = "time (day)")

tscales = [current_parameter(ds, :τ_C),
1/current_parameter(ds, :D)/CTMLM.sec_in_day,
1200/(current_parameter(ds, :V)*current_parameter(ds, :d) + observe_state(ds, :w_e))/CTMLM.sec_in_day,
]
tlabels = ["τ_C", "τ_z", "τ_q"]

for (i, ax) in enumerate(axs)
    lines!(ax, t .- tf, X[:, i])
    # add timescale lines
    for (j, τ) in enumerate(tscales)
        vlines!(ax, [tf+τ, tf+T+τ] .- tf; linestyle = :dash, linewidth = 2,
        label = tlabels[j], color = Cycled(j+1))
    end
    axs[i].xticks = WilkinsonTicks(10; k_min = 8, k_max = 15)
    axs[i].ylabel = string(obs[i])
end
axislegend(axs[4]; nbanks = 2)
xlims!(axs[4], 0, 2T)
xlims!(axs[8], 0, 2T)
display(fig)

wsave(plotsdir("ctbbl", "limitcycle.png"), fig)


# %% Plot CLT as a function of SST and q_b
fig = Figure()
ax = Axis(fig[1,1])
ax.xlabel = "SST"
ax.ylabel = "q_b"
sstgrid = 295:0.1:305.0
qgrid = 10:0.01:22.0
lclm = let
    RH = @. clamp(qgrid' / CTMLM.q_saturation(sstgrid), 0.0, 1.0)
    Tadj = @. sstgrid - 55.0
    c = CTMLM.cₚ/CTMLM.g
    @. c*(Tadj - (1/Tadj - log(RH)/2840.0)^(-1))
end

hm = heatmap!(ax, sstgrid, qgrid, lclm; colorrange = (500, 2000))
cb = Colorbar(fig[1,2], hm; colorrange = (500, 2000))
cb.label = "LCL"

# overplot limitcycle
x, y = X[1:T_i, 3], X[1:T_i, 5]
lines!(ax, x, y; color = "red", linestyle = :dash)

# scatter triangles that show direction
# using Statistics
# ir = round.(Int, range(1, length(x); length = 20))
# xdir = x[ir]
# ydir = y[ir]
# φ = atan.(ydir .- mean(y), xdir .- mean(x)) .- π/2
# scatter!(ax, xdir, ydir; rotation = φ, marker = Makie.Polygon(Point2f[(-1, -1), (2, 0), (-1, 1)]), markersize = 10)


fig
