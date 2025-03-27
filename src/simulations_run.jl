using DynamicalSystems, ProgressMeter, Distributions

function random_initial_conditions(density::Int, extra::Int = 5)
    Cgrid = range(0, 1; length = density)
    Tgrid = range(270, 310; length = density)
    ics = Dict{Symbol, Float64}[]
    for C in Cgrid, T in Tgrid
        for _ in 1:extra # add some small variability to each IC
            r() = 1 + 2e-2randn()
            RH0 = clamp(0.7r(), 0, 1)
            s = (T - 2.0)*r() # in equilibrium s is about 1-2 Kelvin less than SST.
            q = RH0*CTMLM.q_saturation(T)
            z = 1200*r() # this doesn't matter
            u0 = Dict(:SST => T, :C => C, :z_b => z, :q_b => q, :s_b => s)
            push!(ics, u0)
        end
    end
    return ics
end

function sample_parameters!(ds, distributions::Dict{<:Any, <:Distribution}, i)
    for (p, dist) in pairs(distributions)
        try
            set_parameter!(ds, p, rand(dist))
        catch err # it doesn't have parameter
        end
    end
    return ds
end
function sample_parameters!(ds, params::Vector{<:Dict}, i)
    set_parameters!(ds, params[i])
    return ds
end

folders_from_params(::Dict{<:Any, <:Distribution}) = ["random_params_sims",]
folders_from_params(::Vector{<:Dict}) = ["observed_params_sims",]

function multistability_analysis(ds::DynamicalSystem, ics; kw...)
    mapper = recurrences_mapper_ctbbl(ds)
    # this try-catch block ignores parameter values with huge stickiness
    # (right after bifurcation) that requires much larger recurrences thresholds.
    fractions = labels = convergence = nothing
    while isnothing(labels)
        try
            fractions, labels, convergence = convergence_and_basins_fractions(mapper, ics; show_progress = false)
        catch err
        end
    end
    attractors = extract_attractors(mapper)
    return fractions, labels, convergence, attractors
end

function recurrences_mapper_ctbbl(ds; kw...)
    # use finer grid for attractor recurrences algorithm
    # This is a special grid that only resolves things
    # in T and C; the other three variables we don't have to resolve
    # because all states that are different are different in temperature.
    # We could rewrite this as a projected system but in the end we want the
    # attractors in the full 5 dimensional space, so we might as well leave it as is.
    Cgrid = range(0, 1; length = 151)
    Tgrid = range(270, 310; length = 151)
    grid = (Tgrid, Cgrid, (0.0:4000:4000), (250.0:100:350), (1.0:50:60),)

    mapper = AttractorsViaRecurrences(ds, grid;
        sparse = true, consecutive_lost_steps = 100,
        consecutive_recurrences = 2000, attractor_locate_steps = 400,
        Δt = 10*0.05, kw...
    )
    return mapper
end

function multistability_analysis_auto(ds::DynamicalSystem; density = 21, extra = 3)
    ics = random_initial_conditions(density, extra)
    multistability_analysis(ds, ics)
end

function continuation_analysis(ds, pcurve; show_progress = false, density = 11, extra = 3, kw...)
    # Here we augment the standard matching by distance so that limit
    # cycles are not matched; rather they are a different attractor.
    function centroid_and_length(A, B)
        # first check we have a fixed point and limit cycle. We do this by
        # checking if there are different lengths and one of the two is 1 (fixed point)
        if length(A) != length(B) && any(isequal(1), length.((A, B)))
            return Inf
        end
        # otherwise we need a weighted euclidean distance
        # note the weights assume a fixed ordering of the state variables
        # as I haven't implemented yet a named column syntax to the state space set!
        weights = (300.0, 1.0, 1200.0, 300.0, 10.0)
        d = maximum(i -> abs( ( mean(A[:, i]) - mean(B[:, i]) )/weights[i] ), 1:5)
        return d
    end

    # The threshold means a difference of about 20% from the normalizing value
    matcher = MatchBySSSetDistance(; distance = centroid_and_length, threshold = 0.2)
    mapper = recurrences_mapper_ctbbl(ds; kw...)
    ascm = AttractorSeedContinueMatch(mapper, matcher)

    ics = random_initial_conditions(density, extra)
    fractions_cont, attractors_cont = global_continuation(ascm, pcurve, ics; show_progress)
    return fractions_cont, attractors_cont
end

function multiparameter_multistability_analysis(
        ebm::DynamicalSystem, paramsampler;
        foldergroup = folders_from_params(paramsampler), density = 21, extra = 5, N = 1000,
        observables_to_obtain = [
            :C, :SST, :q_b, :z_b, :s_b, :CTRC, :Ld, :Lnet, :LHF, :CLT, :SHF, :T₊, :ASW
        ],
        filename, keep_mean = true,
    )
    # Step 1: input dynamical system and parameters;
    # we make copies to parallelize
    savepath = datadir(foldergroup...)
    dss = [deepcopy(ebm) for _ in 1:Threads.nthreads()]
    mtk = referrenced_sciml_model(ebm)
    paramsampler isa Vector && (N = length(paramsampler))

    @showprogress desc="Computing..." Threads.@threads for _prog in 1:N # run many simulations!
    # for _prog in 1:N # run many simulations!

        ds = dss[Threads.threadid()]
        # Step 2: decide the input parameters based on input distributions
        # input could be probability distributions or vectors of observed data
        sample_parameters!(ds, paramsampler, _prog)

        # Step 3: create some initial conditions variability
        ics = random_initial_conditions(density, extra)

        # Step 4: multistability analysis
        fractions, labels, convergence, attractors = multistability_analysis(ds, ics)

        # Step 5: obtain all observables to save; save equations is probably too complex
        observables = extract_observables(ds, attractors, observables_to_obtain; keep_mean)

        # Step 6: save output; but be conservative
        # TODO: once https://github.com/SciML/ModelingToolkit.jl/issues/3444 is fixed,
        # I can use the simple one-liner
        params_names = Symbol.(ModelingToolkit.parameters(mtk))
        params_values = current_parameter.(ds, params_names)
        params = Dict(params_names .=> params_values)
        # params = Dict(Symbol.(ModelingToolkit.parameters(mtk)) .=> current_parameters(ds))
        output = @strdict fractions attractors params observables
        output = Dict{String, Any}(output)
        output["summary"] = dynamical_system_summary(ds) # with this we can inprinciple recreate everything

        tagsave(joinpath(savepath, filename, string(hash(current_parameters(ds))))*".jld2", output; warn = false)
    end
end

using Statistics: mean

function extract_observables(ds::DynamicalSystem, attractors::Dict, observables; keep_mean = false)
    observables_values = Dict(Symbol(k) => Float64[] for k in observables)
    for A in values(attractors)
        # first, select the points
        if length(A) < 4
            us = A
        else # limit cycle; get states of max SST, max C, and overall mean
            Cs = observe_state.(Ref(ds), :C, vec(A))
            SSTs = observe_state.(Ref(ds), :SST, vec(A))
            j1 = argmin(Cs)
            j2 = argmax(Cs)
            j3 = argmin(SSTs)
            j4 = argmax(SSTs)
            us = A[[j1, j2, j3, j4]]
            keep_mean && push!(us, mean(A))
        end
        # then for each point get the requested observable
        for u in us
            for (p, container) in observables_values
                push!(container, observe_state(ds, p, u))
            end
        end
    end
    return observables_values
end
