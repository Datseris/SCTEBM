include("theme.jl") # to enable plotting during processing

function has_bottleneck(attractors, convergence)
    length(attractors) == 1 || return false # specific to my system: only if 1 attractor exists there is a bottleneck
    if any(≤(0), convergence)
        # TODO: need to fix this in Attractors.jl because it should never happen!
        convergence = filter(>(0), convergence)
    end
    timescale = 50.0 # specialized to my system
    distr = log.(timescale, convergence)
    return count(distr .> 1.1)/length(distr) > 0.25 # want at least 25% more than this threshold
end

has_limitcycle(attractors) = any(A -> length(A) > 2, values(attractors))

function has_fixedpoint(attractors, Cidx = 2, condition = x -> true)
    for A in values(attractors)
        length(A) > 4 && continue # skip limit cycles
        condition(A[end][Cidx]) && return true
    end
    return false
end

has_stratocumulus(attractors, Cidx = 2) = has_fixedpoint(attractors, Cidx, x -> x > 0.5)
has_cumulus(attractors, Cidx = 2) = has_fixedpoint(attractors, Cidx, x -> x < 0.5)

using DynamicalSystems
entropy = DynamicalSystems.entropy
using Random: shuffle!
function rel_mut_info(x, y, bins = 20; trials = 10_000)
    est = ValueHistogram(bins)
    x = copy(x)
    Hx = entropy(est, x)
    Hy = entropy(est, y)
    Hxy = entropy(est, StateSpaceSet(x, y))
    m = Hx + Hy - Hxy
    null = zeros(trials)
    for i in eachindex(null)
        shuffle!(x)
        Hxy = entropy(est, StateSpaceSet(x, y))
        null[i] = Hx + Hy - Hxy
    end
    μ, σ = mean(null), std(null)
    return abs(m - μ)/3σ
end

using ProgressMeter
using OrderedCollections: OrderedDict

# Perform the analysis of multiple saved files, eaching being
# a multistability analysis performed and saved for a particular parameter combination.
function multiparameter_multistability_process(input;
        foldergroup = ["sims", "random_params_sims"],
        parameters_to_obtain = [:U, :D, :RH₊, :δ_Δ₊T, :δ_FTR, :CO2],
        observables_to_obtain = [
            :C, :SST, :q_b, :z_b, :s_b, :CTRC, :Ld, :Lnet, :LHF, :CLT, :SHF, :T₊, :ASW,
        ],
        delete_files = false, plot_density = true, plot_rmi = true, trials = 1000,
        xMI = parameters_to_obtain, yMI = [:SST, :C], estimate_rmi = true,
    )
    # init
    filesfolder = datadir(foldergroup..., savename(input))
    parameters_values = OrderedDict(k => Float64[] for k in parameters_to_obtain)
    observables_values = OrderedDict(k => Float64[] for k in observables_to_obtain)
    output = Dict{String, Any}(string(k) => v for (k, v) in input)
    N = length(readdir(filesfolder))
    natt = zeros(Int, N) # number of attractors in each simulation
    limitcycle = zeros(Bool, N) # boolean for whether there is a limit cycle
    cumulus = zeros(Bool, N) # boolean for whether there is a Cu fixed point
    stratocumulus = zeros(Bool, N) # boolean for whether there is a St fixed point
    diverged = zeros(Float64, N) # fraction of diverged i.c.

    # load
    @showprogress for (i, file) in enumerate(readdir(filesfolder; join = true))
        # load only the necessary stuff
        observables, params, attractors, fractions = wload(file, "observables", "params", "attractors", "fractions")
        # get basic info
        natt[i] = length(attractors)
        diverged[i] = get(fractions, -1, 0.0)
        limitcycle[i] = has_limitcycle(attractors)
        cumulus[i] = has_cumulus(attractors)
        stratocumulus[i] = has_stratocumulus(attractors)

        # get all observables in the same containers
        for O in observables_to_obtain
            # TODO: If O isn't a key of the observables, I can instantiate the dynamical system
            append!(observables_values[O], observables[O])
        end
        # Also get each parameter in the same containers, and I need to duplicate parameters
        m = length(last(first(observables))) # tells us how many states are overall
        for P in parameters_to_obtain
            # We use `get` here in case the particular parameter doesn't exist
            append!(parameters_values[P], fill(get(params, P, NaN), m))
        end
    end

    # Basic quantities
    output["diverged%"] = 100mean(diverged)
    output["limitcycle%"] = 100count(limitcycle)/N
    output["cumulus%"] = 100count(cumulus)/N
    output["stratocumulus%"] = 100count(stratocumulus)/N
    output["multistable%"] = 100count(>(1), natt)/N
    output["quartiles_SST"] = quantile(observables_values[:SST], [0.25, 0.5, 0.75])
    output["quartiles_C"] = quantile(observables_values[:C], [0.25, 0.5, 0.75])
    output["quartiles_z"] = quantile(observables_values[:z_b], [0.25, 0.5, 0.75])
    output["quartiles_q"] = quantile(observables_values[:q_b], [0.25, 0.5, 0.75])
    output["quartiles_s"] = quantile(observables_values[:s_b], [0.25, 0.5, 0.75])

    # densities
    plot_density && plot_densities(input, observables_values, foldergroup)

    # mutual infos
    if estimate_rmi
        rMI = Dict(y => Dict{Symbol, Float64}() for y in yMI)
        for (i, xname) in enumerate(xMI)
            x = parameters_values[xname]
            for (j, yname) in enumerate(yMI)
                y = observables_values[yname]
                any(isempty, (x, y)) && continue # skip non-existent parameters
                rmi = rel_mut_info(x, y; trials)
                rMI[yname][xname] = rmi
            end
        end
        output["rMI"] = rMI
        plot_rmi && plot_rmi_scatters(input, observables_values, parameters_values, xMI, yMI, rMI, foldergroup)
    end

    # Delete the individual files
    if delete_files
        try
            rm(filesfolder; recursive = true, force = true)
        catch
        end
    end

    # save the final output
    outfile = datadir(foldergroup..., savename(input, "jld2"))
    wsave(outfile, output)

    return output
end

function plot_densities(input, observables_values, foldergroup)
    fig = Figure(size = (1200, 1000))
    ci = 1
    for (k, v) in observables_values
        k == :s_b && continue
        i, j = Tuple(CartesianIndices((3,4))[ci])
        mS = string(round(mean(v); sigdigits = 3))
        ax, den = density(fig[i, j], v;
            color = (COLORS[ci], 0.5), strokecolor = COLORS[ci], strokewidth = 3
        )
        textbox!(ax, mS)
        hideydecorations!(ax)
        ylims!(ax, 0, nothing)
        xlims!(ax, density_limits(k))
        ax.titlefont = :regular
        ax.xlabel = String(k)
        ci += 1
    end

    figuretitle!(fig, savename(input; connector = ", "))
    colgap!(fig.layout, 25)
    display(fig)
    wsave(plotsdir(foldergroup..., savename("distributions", input, "png")), fig)
    return
end

function plot_rmi_scatters(input, observables_values, parameters_values, xMI, yMI, rMI, foldergroup)
    fig, axs = axesgrid(length(yMI), length(xMI);
        size = (1200, 400), sharex = true, sharey = true,
        xlabels = string.(xMI), ylabels = string.(yMI)
    )

    for (i, xname) in enumerate(xMI)
        x = parameters_values[xname]
        if xname == :D # special normalization just for D
            x = x ./ 1e-6
            axs[length(yMI), i].xlabel = "D (× 10⁻⁶)"
        end
        for (j, yname) in enumerate(yMI)
            y = observables_values[yname]
            rmi = rMI[yname][xname]
            kw = (color = (COLORS[i], 0.25), strokecolor = (color = (COLORS[i], 0.1)), strokewidth = 2)
            ax = axs[j, i]
            scatter!(ax, x, y; kw...)
            textbox!(ax, string(round(rmi; sigdigits = 2)))
            # Store result in the `output`
            if xname == :δ_FTR
                ax.xticks = [282, 287, 292]
            end
        end
    end

    figuretitle!(fig, savename(input; connector = ", "))
    display(fig)
    wsave(plotsdir(foldergroup..., savename("rmiscatters", input, "png")), fig)

    return
end


function process_multiparam_multistability_analysis(input, params, foldergroup;
        force = false, delete_files = true,
        starting_parameters = Dict(), kw_process = NamedTuple(), kw_analysis = NamedTuple(),
    )
    # skip if the whole simulation has been performed already
    outfile = datadir(foldergroup..., savename(input, "jld2"))
    if isfile(outfile)
        force || return
    end
    println("Currently running simulation with folder configuration:")
    println(foldergroup)
    println("and model configuration")
    display(input)
    # otherwise, perform the whole analysis
    ebm, eqs = ctmlm_setup(; input..., starting_parameters)
    multiparameter_multistability_analysis(ebm, params;
        filename = savename(input), foldergroup, kw_analysis...
    )
    # then aggregate the individual parameter conbimation simulations
    # into a single file and delete the individual files
    multiparameter_multistability_process(input; foldergroup, delete_files, kw_process...)
end

# Make sure the aspects ordering is like in the paper table
# and make sure all match the options!
const ASPECTS_OPTIONS = OrderedDict(
    "ftrgrad" => [:none, :weak, :strong],
    "invfix" => [:difference, :temperature],
    "co2" => [1, 2, 3],
    "ΔF" => [:ctrc, :three_layer, :Gesso2014],
    "Ld" => [:three_layer, :fixed],
    "entrain" => [:Stevens2006, :Gesso2014],
)

"""
    aspects_to_option_idx(row)

Given a `DataFrame` row, or in general a named tuple
mapping aspects to their options, return the integer
(1, 2, 3) corresponding to that option
"""
function aspects_to_option_idx(row)
    idxs = zeros(Int, length(ASPECTS_OPTIONS))
    # Ensure there is an informative error if any of these is nothing
    for i in eachindex(idxs)
        aspect, options = collect(ASPECTS_OPTIONS)[i]
        idx = findfirst(isequal(row[aspect]), options)
        if isnothing(idx)
            @show row[aspect] options
            error("index nothing")
        end
        idxs[i] = idx
    end
    # Same but in a one liner:
    # [findfirst(isequal(row[aspect]), options) for (aspect, options) in pairs(aspects)]
    return idxs
end
