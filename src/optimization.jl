export integer_parameter_index, integer_state_index, optimize_parameters

# Optimization code interface:
# I provide a function that takes in a dynamical system and outputs a
# vector of numbers. This is the objective.
# The surrounding code makes a wrapper function that given parameters it
# reinits the dynamical system at these parameters and calls the user function.

# TODO: Port this to Optimization.jl to not be tied down to BlackBoxOptim.jl

import BlackBoxOptim
# Alternatives: Hyperopt.jl, Metaheuristics.jl, Optim.jl (NelderMead)
# SciML/PSOGPU.jl, .... I should desparately move to another library
# because BlackBoxOptim.jl is one of the most poorly documented libraries I've seen.
# Like, the keywords of the central function are not even listed.

# is there any scientific comparison (limitations, stability, etc.)?
# "kind of" because then it's very dependent on the hyperparmaeters
# you can tweak the hyperparameters of a PSO and get something that searches widely
# or something that searches narrowly, so it's hard to make a very general statement

# see https://discourse.julialang.org/t/comparison-of-derivative-free-black-box-optimization-algorithms-in-julia/109614

"""
    optimize_parameters(objective, ds::DynamicalSystem, parameter_ranges::AbstractDict)

Find the parameter values of `ds` that minimize the given `objective(ds)` function
using global optimization. The parameters are provided as a dictionary mapping
parameter indices to ranges (min, max) to search in.

Return the optimal parameters as a dictionary and the fitness of the fit,
i.e., the minimum value of `objective` achieved by the returned parameters.
"""
function optimize_parameters(objective, ds::DynamicalSystem, parameter_ranges::AbstractDict)
    ranges = collect(values(parameter_ranges))
    psymbols = collect(keys(parameter_ranges))
    function f(pvalues) # function given to the the optimizer
        for i in eachindex(pvalues)
            set_parameter!(ds, psymbols[i], pvalues[i])
        end
        return objective(ds)
    end
    bbres = BlackBoxOptim.bboptimize(f; SearchRange = ranges, TraceMode = :silent, MaxTime = 10.0)
    fitted_params = BlackBoxOptim.best_candidate(bbres)
    fitness = BlackBoxOptim.best_fitness(bbres)
    params_dict = OrderedDict(psymbols .=> fitted_params)
    return params_dict, fitness
end

# this is necessary because of requiring the parameters to
# be changed in the optimization library.

"""
    integer_parameter_index(symbol, ds::DynamicalSystem) → i::Int

Convert the given symbolic variable representing a parameter to its integer
index in the parameter container ([`current_parameters`](@ref)).
Return `nothing` if `ds` doesn't reference a symbolic model
or the model does not have the given `symbol` as parameter.
"""
integer_parameter_index(s, ds::DynamicalSystem) = integer_parameter_index(s, referrenced_sciml_model(ds))
function integer_parameter_index(symbol, model)
    isnothing(model) && return nothing
    findfirst(isequal(symbol), parameters(model))
end

"""
    integer_state_index(symbol, ds::DynamicalSystem) → i::Int

Convert the given symbolic variable representing a state (dynamic) variable to its integer
index in the state vector ([`current_state`](@ref)).
Return `nothing` if `ds` doesn't reference a symbolic model
or the model does not have the given `symbol` as a state variable.
"""
integer_state_index(s, ds::DynamicalSystem) = integer_state_index(s, referrenced_sciml_model(ds))
function integer_state_index(symbol, model)
    isnothing(model) && return nothing
    findfirst(isequal(symbol), states(model))
end
