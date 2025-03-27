# Script for running SCTEBM simulations given some input parameter distributions
# and potential model aspects.
# It outputs a bunch of information regarding the system attractors,
# basin fractions, global stability, as well as recording a bunch of observables.
# You can run this with or without CO2: simply comment out the CO2 distribution
# in the parameters and the CO2-related aspects in the `dict_list` call.
# To run only a specific combination just call the final line in the script with `single_input`
using DrWatson
@quickactivate
using DynamicalSystems
using ConceptualClimateModels
using Statistics
using Distributions
include(srcdir("CloudToppedMixedLayerModel", "CloudToppedMixedLayerModel.jl"))
include(srcdir("ctmlm_setups.jl"))
include(srcdir("simulations_run.jl"))
include(srcdir("simulations_process.jl"))

###########################################################################################
# Inputs
###########################################################################################
# Input 1: variants configuration

# this is used for a single variant
single_input = Dict(
    :cooling => :q_x,
    :entrain => :Gesso2014,
    :co2 => 1,
    :invfix => :difference,
    :ftrgrad => :weak,
    :Ld => :three_layer,
    :ΔF => :ctrc,
)

# this is used for multiple variants
many_inputs = dict_list(Dict(
    :cooling => :q_x,
    :ftrgrad => :weak,
    :entrain => :Gesso2014,
    :Ld => [:three_layer, :fixed],
    :ΔF => [:Gesso2014],
    :invfix => [:difference, :temperature],
    :co2 => [1, 2, 3],
))

# Input 2: distributions of parameters
# It's a dictionary mapping named parameters (symbols) to distributions
distributions = Dict(
    :D => Uniform(1e-6, 5e-6),
    :U => Uniform(5.0, 9.0), # the exaggerated range of U takes into account d variability
    :δ_Δ₊T => Uniform(1, 10),
    :RH₊ => Uniform(0.1, 0.4),
    :δ_FTR => Uniform(0, 10),
    :CO2 => Uniform(200, 1600), # doesn't matter that it's not LOG, we use Mutual Info!
)

# Input 3: folder name to save to
foldergroup = ["sims", "random_params"]

# The rest is the automated running code that will do the whole
# pipeline of running the model variant through the random multi-parameter
# steadystate (multistability) analysis and then estimate the rMIs.
function run_sim(input)
    process_multiparam_multistability_analysis(
        input, distributions, foldergroup;
        force = false, # whether to overwrite existing (aggregate) data
        # multistability multiparameter keywords
        kw_analysis = (density = 5, extra = 2, N = 1000),
        # keywords for processing
        # make plots true if you are running a single aspect
        kw_process = (plot_rmi = true, plot_density = true, delete_files = true,),
    )
end

# run_sim(single_input) # uncomment for single variant
map(input -> run_sim(input), many_inputs) # uncomment for multiple variants