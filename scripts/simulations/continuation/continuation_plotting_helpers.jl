using Attractors

# Helper functions
TRANSPARENCY = Observable(0.25) # can update this observable in any analysis script

state_colors = Dict(
    "Sc ✓" => @lift((COLORS[3], $TRANSPARENCY)),
    "Sc ↓" => @lift((COLORS[5], $TRANSPARENCY)),
    "Sc ↶" => @lift((COLORS[6], $TRANSPARENCY)),
    "Cu" => @lift((COLORS[1], $TRANSPARENCY)),
    "LC" => @lift((COLORS[2], $TRANSPARENCY)),
    "X" => @lift((COLORS[4], $TRANSPARENCY))
)

state_markers = Dict(
    "Sc ↓" => MARKERS[3],
    "Sc ✓" => MARKERS[2],
    "Sc ↶" => MARKERS[6],
    "Cu" => MARKERS[1],
    "LC" => MARKERS[4],
    "X" => MARKERS[5]
)

function classify_series(series, collapsed)
    idx = findfirst(!isempty, series)
    isnothing(idx) && return "X"
    A = series[idx]
    if length(A) > 4
        return "LC"
    elseif A[length(A), :C] > 0.5
        # First see if it starts with Sc in case of recovery
        if idx > 1
            return "Sc ↶"
        elseif collapsed
            return "Sc ↓"
        else
            return "Sc ✓"
        end
    elseif A[length(A), :C] ≤ 0.5
        return "Cu"
    else
        return "X"
    end
end

# given a global continuation output for a particular attractor ID,
# classify it according to one of the 5 categories of the climate change figures
# (same as keys of the state markers)
function classify_series(series)
    # we first check where to start, if we even start
    idx_start = findfirst(!isempty, series)
    isnothing(idx_start) && return "X"
    # we then check the cases one by one. First, Cu or LC cases do not care about
    # when they started. So we check immediatelly:
    if length(series[idx_start]) > 4
        return "LC"
    elseif series[idx_start][1, :C] < 0.5
        return "Cu"
    end
    # we then do a sanity check that Sc state actually exists
    if mean(series[idx_start][:, :C]) < 0.5
        return "X"
    end
    # we then check for non-starting Sc
    if idx_start > 1
        return "Sc ↶"
    end
    # we then check if the simulation persists until the end or not
    idx_end = findfirst(isempty, series)
    if isnothing(idx_end)
        return "Sc ✓"
    end
    # then check if there is recovery
    idx_final = findlast(!isempty, series)
    if idx_final == length(series)
        return "Sc ↶"
    else
        return "Sc ↓"
    end
    return "X"
end

# Plotting code
function plot_attractor_series!(axs, observables, ds, attractors_cont)
    attractors_series = continuation_series(attractors_cont, StateSpaceSet{5, Float64}(; names = fill(:x, 5)))

    # use the classification to create consistent coloring
    classification = Dict(k => classify_series(series) for (k, series) in attractors_series)
    colors = Dict(k => state_colors[v] for (k, v) in classification)
    markers = Dict(k => state_markers[v] for (k, v) in classification)
    # project to means
    function observable_mean(ds, obs, A)
        isempty(A) && return NaN
        if isnothing(ds)
            mean(A[:, obs])
        else
            mean(observe_state(ds, obs, u) for u in A)
        end
    end
    # plot continuation
    t = 0:length(attractors_cont)-1
    for (ax, obs) in zip(axs, observables)
        cont = map(attractors_cont) do attractors
            Dict(k => observable_mean(ds, obs, A) for (k, A) in attractors)
        end
        Attractors.plot_continuation_curves!(ax, cont, t; colors, markers, add_legend = false)
    end
    # I can try to connect the dots here by adding collapse arrows but damn its hard!
    # anyways, record whether this particular simulation has a collapse or not and return
    has_collapse_label = "Sc ↓" ∈ values(classification)
    return has_collapse_label
end

# This function can plot continuations for ANY observable of the SCTEBM!
function load_n_plot_continuation!(fig_loc, input, used, observables = [:C, :SST]; ids = 1:100, add_legend = true)
    foldergroup = ["sims", "continuations"]
    prefix = "used="*join(string.(used), "+")
    name = savename(input)*"_"*prefix
    data = wload(datadir(foldergroup..., name)*".jld2")
    @unpack continuations, param_values = data
    # Create dynamical system (to observe states)
    ds, eqs = sctebm_setup(input)

    # TODO: make ds, make the function generic

    axs = axesgrid!(fig_loc, length(observables), 1; xlabels = "time (a.u.)", ylabels = string.(observables), sharex = true)

    for i in ids
        attractors_cont = continuations[i]
        set_parameters!(ds, param_values[i])
        plot_attractor_series!(axs, observables, ds, attractors_cont)
    end
    rowgap!(fig_loc.layout, 4)

    if add_legend
        cloud_transition_legend!(fig_loc[:, 2])
    end

    return axs
end

function cloud_transition_legend!(figloc)
    state_names = sort(collect(keys(state_colors)))
    filter!(≠("X"), state_names)
    elements = [
        [LineElement(color = state_colors[k][][1]),
        MarkerElement(color = state_colors[k][][1], marker = state_markers[k], markersize = 25)]
        for k in state_names
    ]
    Legend(figloc, elements, state_names)
end


# This function only plots cloud fraction and therefore does not need to initialize the
# dynamical system instance to observe a state.
function load_n_plot_continuation_only_c!(axc, input, used; add_collapse = true, ids = 1:100)
    foldergroup = ["sims", "continuations"]
    prefix = "used="*join(string.(used), "+")
    name = savename(input)*"_"*prefix
    data = wload(datadir(foldergroup..., name)*".jld2")
    @unpack continuations, no_Sc_time = data
    has_collapse = falses(ids)
    global acont = 1
    for (i, attractors_cont) in enumerate(continuations)
        i ∉ ids && continue

        # Fix for loading old data before named dimensions in SSSet
        if typeof(attractors_cont[1]) <: DrWatson.JLD2.SerializedDict
            stype = typeof(StateSpaceSet{5, Float64}(; names = fill(:x, 5)))
            acont = Vector{Dict{Int, stype}}()
            names = [:SST, :C, :z_b, :s_b, :q_b]
            for j in 1:length(attractors_cont)
                dict = Dict(key => StateSpaceSet(val.data; names) for (key, val) in attractors_cont[j].kvvec)
                push!(acont, dict)
            end
            attractors_cont = acont
        end

        acont = attractors_cont
        has_collapse[i] = plot_attractor_series!([axc], [:C], nothing, attractors_cont)
    end
    if add_collapse
        percent = 100count(has_collapse)/length(has_collapse)
        textbox!(axc, "$(percent)% Sc ↓";
            valign = 0.2, halign = 0.1,
        )
    end
    return axc
end
