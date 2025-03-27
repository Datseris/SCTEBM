# Here we use the limit cycle set up and apply seasonal forcing
# we vary some parameters to obtain various different behaviors
# of the forced system

using DrWatson
@quickactivate
using DynamicalSystems
using OrdinaryDiffEq
using ConceptualClimateModels
using CairoMakie
using OrderedCollections: OrderedDict
import ClimateBase
import Dates

# include(srcdir("theme.jl"))
include(srcdir("CloudToppedMixedLayerModel", "CloudToppedMixedLayerModel.jl"))
include(srcdir("ctmlm_setups.jl"))

# This setup, with the parameters set below, results in a limit cycle
w_m = true
ftrgrad = :weak
cfversion = :Singer2023
co2_sst = false
co2_ΔT = true
radiation = :detailed # detailed or merged; for merged we need smaller drag coefficient
q_x = radiation == :detailed
input = @strdict(w_m, q_x, cfversion, co2_sst, radiation, co2_ΔT, ftrgrad)

starting_parameters = Dict(:U => 6.68, :T_FTR_0 => 289.0)

function plot_ds(ds, title = "")
    X, t = trajectory(ds, 1000.0)
    t = t ./ 30 # make it in months
    fig, axs = axesgrid(2, 1;
        size = (1000, 400), xlabels = "time (month)", ylabels = ["SST", "C"], sharex = true,
    )

    lines!(axs[1], t, X[:, 1])
    lines!(axs[2], t, X[:, 2], color = Cycled(2))

    ticks = 1:5:36
    names = @. Dates.monthname(mod1(ticks, 12))

    for ax in axs; ax.xticks = (ticks, names); end

    ax3 = Axis(fig[:, 2];
        xlabel = "C", ylabel = "SST", limits = (0, 1, 270, 310)
    )
    r = length(X)÷2:length(X)
    lines!(ax3, X[r, 2], X[r, 1]; color = Cycled(3))
    colsize!(fig.layout, 1, Relative(0.75))
    figuretitle!(fig, title)
    display(fig)
    fig
end


# %% case 1: just limit cycle
ds, eqs = ctmlm_setup(input; starting_parameters)

fig = plot_ds(ds, "limit cycle (no seasonal forcing)")

# wsave(plotsdir("ctbbl", "seasonal", "limitcycle.png"), fig)

# %% Turning on seasonal forcing

# this additional equation will make the insolation follow the seasonal cycle
@register_symbolic ClimateBase.insolation(t, phi)
extra_eqs = [
    CTMLM.S ~ ClimateBase.insolation(CTMLM.t, -28),
]
ds, eqs = ctmlm_setup(input; starting_parameters, extra_eqs)

fig = plot_ds(ds, "seasonal forcing = on, bursting")

# wsave(plotsdir("ctbbl", "seasonal", "default_seasonal.png"), fig)

# %% Playground for the other parameters
quasiperiodic_parameters = Dict(
    :U => 6.68,
    :T_FTR_0 => 289.0,
    :D => 2.24e-6,
    :e_e => 1.51,
    :τ_SST => 57.3,
    :τ_C => 1.95,
)

ds, eqs = ctmlm_setup(input; starting_parameters = chaos_parameters, extra_eqs)
λ = lyapunov(ds, 20000; Ttr = 100.0)
λs = lyapunovspectrum(ds, 20000; Ttr = 100.0)
λ1, λ2 = λs

fig = plot_ds(ds, "seasonal forcing, quasiperiodic-like state")
wsave(plotsdir("ctbbl", "seasonal", "quasiperiodic_like.png"), fig)



# %%

fig = plot_ds(ds, "seasonal forcing, periodic state reversals")

wsave(plotsdir("ctbbl", "seasonal", "periodic_reversals.png"), fig)
