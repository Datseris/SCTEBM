# sensitivity analysis of a particular variant.
# like a pre-step to the main aggregate analysis.
using DrWatson
@quickactivate
using DynamicalSystems
using ConceptualClimateModels
using Statistics
using Distributions
include(srcdir("sctebm_setups.jl"))
include(srcdir("simulations_run.jl"))
include(srcdir("simulations_process.jl"))

# Run the simulations and save data to disk
variant = Dict(
    :cooling => :q_x,
    :entrain => :Stevens2006,
    :co2 => 1,
    :invfix => :difference,
    :ftrgrad => :weak,
    :Ld => :three_layer,
    :ΔF_s => :three_layer,
    :cloud_rad => :const,
)

ebm, eqs = sctebm_setup(variant)

distributions = Dict(
    :D => Uniform(1e-6, 5e-6),
    :U => Uniform(5.0, 9.0), # the exaggerated range of U takes into account d variability
    :δ_Δ₊T => Uniform(1, 10),
    :RH₊ => Uniform(0.1, 0.4),
    :δ_FTR => Uniform(0, 10),
    :CO2 => Uniform(200, 1600), # doesn't matter that it's not LOG, we use Mutual Info!
)

foldergroup = ["sims", "single_variant"]

observables_to_obtain = [:C, :SST, :q_b, :z_b]

multiparameter_multistability_analysis(ebm, distributions;
    filename = savename(variant), foldergroup, observables_to_obtain,
    density = 5, extra = 2, N = 500
)

# %% Load and analyze
output = multiparameter_multistability_process(variant;
    foldergroup, delete_files = false, observables_to_obtain
)
