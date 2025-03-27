# This script performs the continuations (climate change scenario simulations)
# and saves data for a given model variant. Another script loads and plots the data later.

using DrWatson
@quickactivate
using DynamicalSystems
using ConceptualClimateModels
using Statistics
using Distributions
using Statistics: mean
include(srcdir("ctmlm_setups.jl"))
include(srcdir("simulations_run.jl"))
include(srcdir("simulations_process.jl"))
include(srcdir("theme.jl"))
using ProgressMeter

###########################################################################################
# Input
###########################################################################################
# Input 1: equations configuration
cooling = :q_x
invfix = :temperature
Ld = :three_layer # top bottom may be nonsensical; I should study it in the GUI
ΔF = :ctrc
co2 = 2
entrain = :Stevens2006
ftrgrad = :none
input = @dict(cooling, invfix, Ld, ΔF, co2, entrain, ftrgrad)

# Input 2: distributions of parameters
# It's a dictionary mapping named parameters (symbols) to distributions
distributions = Dict(
    :D => Uniform(1e-6, 5e-6),
    :U => Uniform(5.0, 9.0), # the exaggerated range of U takes into account d variability
    :δ_Δ₊T => Uniform(1, 10),
    :RH₊ => Uniform(0.1, 0.4),
    :δ_FTR => Uniform(0, 10),
)

# Input 3: change of environmental conditions with "time"
time = 0:14
rates = Dict(:CO2 => 100, :D => -0.1e-6, :U => -0.2)
starts = Dict(:CO2 => 400, :D => 3e-6, :U => 8.0)

# Input 4: which environmental conditions will change with time
used = sort([:CO2]) # always sort this!!!
N = 100 # number of randomized param simulations

###########################################################################################
# Code runs
###########################################################################################
# Create parameter curve for the continuation
pcurve = [Dict{Symbol, Float64}() for t in time]
for (p, t) in zip(pcurve, time)
    for k in used
        p[k] = starts[k] + rates[k]*t
    end
end

# init stuff
ebm, eqs = ctmlm_setup(input)

collapse = fill(false, N)
has_LC = fill(false, N)
no_Sc_time = fill(NaN, N)
only_Cu_time = fill(NaN, N)
continuations = Vector{Any}(undef, N)
param_values = Vector{Any}(undef, N)

progress = ProgressMeter.Progress(N)

# for/while loop for different random parameter samplings
i = 1
while i ≤ N
    sample_parameters!(ebm, distributions, i)

    for k in used
        set_parameter!(ebm, k, starts[k])
    end

    # sometimes it can fail for a particular combination of parameters
    # (e.g., way too slow convergence / late tipping); we just skip these
    attractors = nothing
    try
        fractions, labels, convergence, attractors = multistability_analysis_auto(ebm; density = 11)
    catch err
        continue
    end

    # only start the continuation if we start with a stratocumulus state
    if !has_stratocumulus(attractors)
        continue
    end

    # perform continuation with our special matching
    fractions_cont = attractors_cont = nothing
    try
        fractions_cont, attractors_cont = continuation_analysis(ebm, pcurve; density = 5, extra = 2)
    catch err
        continue
    end

    continuations[i] = attractors_cont

    # record collapse times, if any.
    j = findlast(!has_stratocumulus, attractors_cont)
    if !isnothing(j)
        no_Sc_time[i] = time[j]
    end
    p = findfirst(a -> has_cumulus(a) && !has_stratocumulus(a) && !has_limitcycle(a), attractors_cont)
    if !isnothing(p)
        only_Cu_time[i] = time[p]
    end

    has_LC[i] = any(has_limitcycle, attractors_cont)

    param_values[i] = named_current_parameters(ebm)

    i += 1
    next!(progress)
end

foldergroup = ["sims", "continuations"]
prefix = "used="*join(string.(used), "+")
name = savename(input)*"_"*prefix
data = @strdict(continuations, only_Cu_time, no_Sc_time, has_LC, param_values)
wsave(datadir(foldergroup..., name)*".jld2", data)
