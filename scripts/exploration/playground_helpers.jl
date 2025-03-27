include(srcdir("theme.jl"))
include(srcdir("stateflow.jl"))
using GLMakie; GLMakie.activate!()
using OrderedCollections: OrderedDict

function ctmlm_gui(ds::DynamicalSystem, GUI_par, GUI_obs = :dynamic; dt = 0.05, u0s = ctbbl_random_ics(ds))

    observables = obtain_GUI_obs(GUI_obs)
    parameter_sliders = OrderedDict(k => allparams[k] for k in GUI_par)
    tail = round(Int, 100÷dt)
    idxs = [:C, :SST]

    fig, dsobs = interactive_trajectory_timeseries(
        ds, observables, u0s;
        Δt = dt,
        tail,
        idxs,
        # lims = Tuple(plausible_limits(ds, idxs)),
        lims = ((0, 1), (280, 310)),
        timelabel = "time (days)",
        parameter_sliders,
        starting_step = tail,
        linekwargs = (linewidth = 2.0,),
        colors = color_from_u0.(u0s),
    )

    # add nullclines if 2D
    if dimension(ds) == 2
        ax = content(fig[1,1][1,1])
        # add nullclines and update them on parameter change
        ux, uy = nullclines!(ax, ds, limits[1], limits[2])
        on(dsobs.param_observable) do params
            nullclines!(ux[], uy[], ds, limits[1], limits[2])
            notify.((ux, uy))
        end
    end

    step!(dsobs, tail)
    display(fig)
    return fig, dsobs
end

function obtain_GUI_obs(GUI_obs)
    GUI_obs isa Vector && return GUI_obs
    observables = [:C, :SST, :z_b, :q_b, :s_b]
    if GUI_obs == :dynamic
        observables = [:C, :SST, :q_b, :z_b, :s_b]
    elseif GUI_obs == :SST
        observables = [:SST, :ASW, :Ld, :L₀, :LHF, :SHF]
    elseif GUI_obs == :height
        w_e╱D(u) = observe_state(ds, :w_e, u)/current_parameter(ds, :D)
        w_m╱D(u) = observe_state(ds, :w_m, u)/current_parameter(ds, :D)
        cltex(s, q, z) = (zb = CTMLM.cloud_base_height_exact(s, q, z); clamp((z - zb)/z, 0, 1))
        CLTex(u) = observe_state(ds, cltex(CTMLM.s_b, CTMLM.q_b, CTMLM.z_b), u)
        cltBolton(s, q, z) = (zb = CTMLM.cloud_base_height_bolton1980(s, q); clamp((z - zb)/z, 0, 1))
        CLTBolton(u) = observe_state(ds, cltBolton(CTMLM.s_b, CTMLM.q_b, CTMLM.z_b), u)
        cltRomps(s, q, z) = (zb = CTMLM.cloud_base_height_romps2017(s, q); clamp((z - zb)/z, 0, 1))
        CLTRomps(u) = observe_state(ds, cltRomps(CTMLM.s_b, CTMLM.q_b, CTMLM.z_b), u)
        observables = [:z_b, w_e╱D, CLTex, CLTBolton, CLTRomps]
    elseif GUI_obs == :clouds
        observables = [:C, :𝒟, :CTRC, :LHF, :CLT]
    elseif GUI_obs == :ctrc
        # ΔLss(u) = observe_state(ds, CTMLM.C*0.9*CTMLM.σ_SB*(CTMLM.T_t^4 - CTMLM.T_FTR^4), u)
        observables = [:CTRC, :L_c, :L_FTR, :ε_FTR, :ε_c]
    elseif GUI_obs == :emissivity
        logq_b(u) = log(observe_state(ds, :q_b, u))
        logq₊(u) = log(observe_state(ds, :q₊, u))
        observables = [:C, :ε_b, :ε_FTR]
    elseif GUI_obs == :radiation
        observables = [:ASW, :L_FTR, :Ld, :Lnet, :CTRC, :ΔF]
    elseif GUI_obs == :temperature
        observables = [:SST, :T_b, :T_t, :T₊, :T_FTR]
    elseif GUI_obs == :inversion
        observables = [:T_t, :T₊, :Δ₊T, :Δ₊s, :Δ₊q]
    elseif GUI_obs == :humidity
        RH_b = CTMLM.q_b/CTMLM.q₀
        observables = [:q₀, :q_b, :q₊, :q_x]
    elseif GUI_obs == :boundarylayer
        RH_b = CTMLM.q_b/CTMLM.q₀
        Δs(u) = ((sp, s0) = observe_state.(ds, (:s₊, :s₀), Ref(u)); sp - s0)
        Δq(u) = ((sp, s0) = observe_state.(ds, (:q₊, :q₀), Ref(u)); sp - s0)
        observables = [:z_b, :CLT, :q_b, :s_b, RH_b]
    end
    return observables
end


allparams = OrderedDict(
    # timescales:
    :τ_SST => 1:0.1:100.0,
    :τ_C => 0:0.01:10.0,
    # Environmental conditions
    :S_0 => (0:1:500),
    :SST_X_0 => -20:20,
    :D => (1:0.01:8) .* 1e-6,
    :RH₊ => 0:0.005:0.5,
    :q₊_0 => 0:0.01:5.0,
    :U => (0:0.01:12.0),
    :CO2 => 100:10:3200,
    :δ_Δ₊T => -5:0.1:10,
    :δ_FTR => -10:0.1:10,
    # Atmospheric composition:
    :ε_a => (0:0.01:1),
    :ε_q => (0:0.01:1),
    :ε_CO2 => (0:0.01:1),
    :ε_FTR => (0:0.01:1),
    :ρ => (0.8:0.01:2.0),
    # Clouds:
    :α_C => (0:0.01:1),
    :α_a => (0:0.01:1),
    :ε_C => (0:0.01:1),
    :Cmin => (0:0.01:1),
    :Cmax => (0:0.01:1),
    :𝒟s => (0:0.01:1),
    # Velocity efficiencies:
    :e_e => 0.0:0.01:3.0,
    :e_v => 0:0.01:2.0,
    :e_m => -2:0.01:2.0,
    :c_d => 0.0005:0.00001:0.0015,
    # Free troposphere:
    :ECS_CO2 => (0.0:0.01:10),
    :Δ₊T_C => (0.0:0.01:10),
    :T₊_0 => (0:0.01:20) .+ 280,
    :T_FTR_0 => (0:0.01:50) .+ 260,
    # Auxilary
    :p => 0:0.01:10.0,
    :q_x_rate => 0:0.001:10.0,
    :s_x_0 => -4:0.01:4.0,
)

function ctbbl_random_ics(ds)
    # Some broad spectrum of core initial conditions
    u0s = [Dict(:C => C, :SST => T) for C in [0.05, 0.5, 0.95] for T in [285.0, 295, 305]]

    # Establish generic variability for the rest of the dynamic variables
    if dimension(ds) == 5
        u0s0 = deepcopy(u0s)
        u0s = Dict{Symbol, Real}[]
        for u0 in u0s0
            for _ in 1:3
                # get random variability coefficientt
                r = () -> 1 + 1e-2randn()
                T0 = u0[:SST]
                RH0 = clamp(0.7r(), 0, 1)
                u0[:SST] = T0*r()
                u0[:C] = u0[:C]*r()
                u0[:s_b] = (T0 - 2.0)*r() # in equilibrium s is about 1-2 Kelvin less than SST.
                u0[:q_b] = RH0*CTMLM.q_saturation(T0*r())
                u0[:z_b] = 1200*r()
                push!(u0s, copy(u0))
            end
        end
    end
    return u0s
end

function basic_statespace_plot(ds, dims = [:C, :SST]; total_time = 10.0)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1,1], xlabel = string(dims[1]), ylabel = string(dims[2]))
    u0s = ctbbl_random_ics(ds)
    for u0 in u0s
        reinit!(ds, u0)
        X, t = trajectory(ds, total_time; save_idxs = dims)
        fadelines!(ax, X)
        scatter!(ax, X[end]; color = "red")
    end
    ax.limits = Tuple(plausible_limits(ds, dims))
    return fig
end
