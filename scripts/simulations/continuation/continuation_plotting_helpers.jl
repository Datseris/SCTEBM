# Helper functions
TRANSPARENCY = Observable(0.25)
state_colors = Dict(
    "Sc ✓" => @lift((COLORS[3], $TRANSPARENCY)),
    "Sc ↓" => @lift((COLORS[2], $TRANSPARENCY)),
    "Sc ↶" => @lift((COLORS[6], $TRANSPARENCY)),
    "Cu" => @lift((COLORS[1], $TRANSPARENCY)),
    "LC" => @lift((COLORS[4], $TRANSPARENCY)),
    "X" => @lift((COLORS[5], $TRANSPARENCY))
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
    if length(A) > 1
        return "LC"
    elseif A[end][2] > 0.5 # index 2 = cloud fraction
        # First see if it starts with Sc in case of recovery
        if idx > 1
            return "Sc ↶"
        elseif collapsed
            return "Sc ↓"
        else
            return "Sc ✓"
        end
    else A[end][2] ≤ 0.5
        return "Cu"
    end
end


# Plotting code
function plot_attractor_series!(axs, observables, ds, attractors_cont)
    collapse = !isnothing(findlast(!has_stratocumulus, attractors_cont))

    attractors_series = continuation_series(attractors_cont, StateSpaceSet{5, Float64}())

    # use the classification to create consistent coloring
    classification = Dict(k => classify_series(series, collapse) for (k, series) in attractors_series)
    colors = Dict(k => state_colors[v] for (k, v) in classification)
    markers = Dict(k => state_markers[v] for (k, v) in classification)

    # project to means
    function observable_mean(ds, obs, A)
        isempty(A) && return NaN
        mean(observe_state(ds, obs, u) for u in A)
    end
    # plot continuation
    t = 0:length(attractors_cont)-1
    for (ax, obs) in zip(axs, observables)
        cont = map(attractors_cont) do attractors
            Dict(k => observable_mean(ds, obs, A) for (k, A) in attractors)
        end
        plot_continuation_curves!(ax, cont, t; colors, markers, add_legend = false)
    end
    # I can try to connect the dots here by adding collapse arrows but damn its hard!

    return axs
end

function load_n_plot_continuation!(fig_loc, input, used, observables = [:C, :SST]; ids = 1:100, add_legend = true)
    foldergroup = ["sims", "continuations"]
    prefix = "used="*join(string.(used), "+")
    name = savename(input)*"_"*prefix
    data = wload(datadir(foldergroup..., name)*".jld2")
    @unpack continuations, param_values = data
    # Create dynamical system (to observe states)
    ds, eqs = ctmlm_setup(input)

    # TODO: make ds, make the function generic

    axs = axesgrid!(fig_loc, length(observables), 1; xlabels = "time (a.u.)", ylabels = string.(observables), sharex = true)

    for i in ids
        attractors_cont = continuations[i]
        set_parameters!(ds, param_values[i])
        plot_attractor_series!(axs, observables, ds, attractors_cont)
    end
    rowgap!(fig_loc.layout, 4)

    if add_legend
        state_names = sort(collect(keys(state_colors)))
        elements = [
            [LineElement(color = state_colors[k][][1]),
            MarkerElement(color = state_colors[k][][1], marker = state_markers[k], markersize = 25)]
            for k in state_names
        ]
    end
    Legend(fig_loc[:, 2], elements, state_names)

    return axs
end

function load_n_plot_continuation_only_c!(axc, input, used; add_collapse = true)
    foldergroup = ["sims", "continuations"]
    prefix = "used="*join(string.(used), "+")
    name = savename(input)*"_"*prefix
    data = wload(datadir(foldergroup..., name)*".jld2")
    ds, eqs = ctmlm_setup(input)
    @unpack continuations, only_Cu_time, no_Sc_time = data
    for (i, attractors_cont) in enumerate(continuations)
        plot_attractor_series!([axc], [:C], ds, attractors_cont)
    end
    if add_collapse
        percent = 100count(!isnan, no_Sc_time)/length(no_Sc_time)
        textbox!(axc, "$(percent)% Sc ↓";
            valign = 0.2, halign = 0.1,
        )
    end
    return axc
end
