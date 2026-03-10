# Given a SCTEBM instance, plot the timeseries of the
# time derivatives of the observables

using DrWatson
@quickactivate
using DynamicalSystems
using ConceptualClimateModels
include(srcdir("sctebm_setups.jl"))
include(srcdir("simulations_run.jl"))
include(srcdir("theme.jl"))

ds, eqs = sctebm_setup(;
    Ld = :three_layer,
    ΔF_s = :ctrc,
    entrain = :Stevens2006,
    ftrgrad = :weak,
    cloud_rad = :const,
    cloud_limit = true,
)

set_parameter!(ds, :U, 9)
set_parameter!(ds, :D, 4.5e-6)
set_parameter!(ds, :RH₊, 0.2)
set_parameter!(ds, :δ_Δ₊T, 5)
set_parameter!(ds, :δ_FTR, 5)

u0 = Dict(:SST => 290.0, :C => 1.0)

X, t = trajectory(ds, 150.0, u0; Δt = 0.05)

f = dynamic_rule(ds)

Xderiv = StateSpaceSet(map(u -> f(u, current_parameters(ds), 0.0), X); names = X.names)

# %% plot

fig, axs = axesgrid(3, 2; titles = ["variables", "derivatives"], xlabels = "time (days)", ylabels = ["SST (K)", "C", "q_b (g/kg)"])
for (i, name) in enumerate((:SST, :C, :q_b))
    lines!(axs[i, 1], t, X[:, name])
    if i == 1 # different units for sst
        units = 50 * 4 * 10^6 / 86400
    else
        units = 1
    end
    lines!(axs[i, 2], t, Xderiv[:, name] .* units)
end

axs[1, 2].ylabel = "W/m^2"
axs[2, 2].ylabel = "1/day"
axs[3, 2].ylabel = "g/kg/day"

fig
