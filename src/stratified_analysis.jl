# Stratify data by various aspects + visualization for it
import LsqFit, Statistics

function visualize_stratified_fits(stratf::Function, X, Y; kw...)
    Xs, titles = stratf(X)
    Ys, titles = stratf(Y)
    xlabel = string(X.name)
    ylabel = string(Y.name)
    visualize_stratified_fits(Xs, Ys; titles, xlabel, ylabel, kw...)
end

# TODO: Docstring here
function visualize_stratified_fits(Xs, Ys;
        fitexprs = [(x, p) -> @. p[1] + p[2]*x], fitstrs = ["a + b*X"], colorid = 1, alpha = 0.02,
        p0s = [rand(5) for _ in fitexprs],
        titles = nothing,
        xlabel = "X", ylabel = "Y", type = :scatter, kw...,
    )

    fig, axs = axesgrid(1, length(Xs);
        titles, xlabels = xlabel, ylabels = ylabel, sharex = true,
        sharey = true, size = (length(Xs)*400, 450),
    )
    ps = Array{Vector{Vector{Float64}}}(undef, size(axs)...)
    # populate parameter container so that we can use `push!`
    for i in eachindex(ps); ps[i] = Vector{Float64}[]; end
    # loop over stratification
    ymin, ymax = Inf, -Inf
    for (j, ax) in enumerate(axs)
        x, y = Xs[j], Ys[j]
        # manual limits for y axis because of labels of fitted functions
        ymin = min(ymin, minimum(y))
        ymax = max(ymax, maximum(y))
        ylims!(ax, ymin, ymax + 0.2(ymax-ymin))
        # plot data
        if type == :scatter
            scatter!(ax, x, y; color = (COLORS[colorid], alpha), kw...)
        elseif type == :hex
            hexbin!(ax, x, y;
                colormap = [(COLORS[colorid], alpha), COLORS[colorid]], bins = 40,
                threshold = 1, kw...
            )
        end
        # loop over fitting functons
        # TODO: use overplot_fits! function from fitting.jl
        k = colorid+1
        for (i, f) in enumerate(fitexprs)
            # perform fit
            fit = LsqFit.curve_fit(f, x, y, p0s[i])
            pfit = LsqFit.coef(fit) # best parameters
            xsort = range(minimum(x), maximum(x); length = 100)
            ysort = f(xsort, pfit)
            # plot fit while reporting error
            e = nrmse(y, f(x, pfit))
            lines!(axs[j], xsort, ysort;
                color = (COLORS[k], 0.9), linewidth = 4, linestyle = LINESTYLES[i],
                label = fitstrs[i]*", nrmse = $(string(round(e; sigdigits = 3)))"
            )
            # push fitted parameters to stored container
            push!(ps[j], pfit)
            k += 1
        end
        axislegend(ax; position = :lt)
    end
    # manual limits for y axis because of labels of fitted functions
    ylims!(axs[1], ymin, ymax + 0.2(ymax-ymin))

    return fig, ps
end

function stratify_by_hemisphere(C)
    descr = ["NH", "SH", "BOTH"]
    Cn, Cs = vec.(hemispheric_functions(C))
    Ca = vcat(Cn, Cs)
    return (Cn, Cs, Ca), descr
end

# note that this function is defined w.r.t. global variables by default!
function stratify_by_sfc_type(C; othreshold = 0.5, ithreshold = 0.5, I = I, O = O)
    descr = ["ice (> $(ithreshold))", "ocean (> $(othreshold))", "land (else)"]
    ice_idxs = findall(I .> ithreshold)
    oce_idxs = findall(O .> othreshold)
    lan_idxs = setdiff(eachindex(C), vcat(ice_idxs, oce_idxs))
    Cice = C[ice_idxs]
    Coce = C[oce_idxs]
    Clan = C[lan_idxs]
    return (Cice, Coce, Clan), descr
end

function stratify_by_cloud_regime(C;
        tropical_ascend = (-5, 15), O = ones(C), othreshold = 0.5,
        midlattitudes_north = (50, 65),
        midlattitudes_south = (-65, -50),
    )
    # subsidence and low clouds (stratocumulus). These boxes come from Fig. 1a of
    # Qu, X., Hall, A., Klein, S.A. et al. On the spread of changes in marine low cloud cover
    # in climate model simulations of the 21st century. Clim Dyn 42, 2603–2626 (2014).
    # https://doi.org/10.1007/s00382-013-1945-z
    # which is cited by the review article on low clouds and controlling factors (Klein, 2017)
    # All boxes are to the west of continents in the subtropics and are 40 by 20 degrees.
    lowcloudboxes = [ # lower-left corner
        (-110+360, -30), # Peru
        (-155+360, 15), # california
        (-25+360, -30), # Namibia
        (-55+360, 10), # Canary
        (75, -38), # Australia
    ]

    lowcloudv = eltype(C)[]
    for box in lowcloudboxes
        s = C[Coord(Lon(Between(box[1], box[1]+40)), Lat(Between(box[2], box[2]+20)))]
        append!(lowcloudv, vec(s))
    end

    # tropical ascend
    accessor = Coord(Lat(Between(tropical_ascend)))
    tropical_idxs = findall(O[accessor] .> othreshold)
    tropical = vec(C[accessor][tropical_idxs])

    # midlatitudes
    # Blanco et al. use between 50-65 degrees for "midlatitudes".
    # one could also use just extra tropics such as 30-polewards
    midlat = eltype(C)[]
    for midlattitudes in (midlattitudes_north, midlattitudes_south)
        accessor = Coord(Lat(Between(midlattitudes)))
        midlat_idxs = findall(O[accessor] .> othreshold)
        m = vec(C[accessor][midlat_idxs])
        append!(midlat, m)
    end

    descr = ["subsidence low clouds", "tropical ascend", "storm tracks"]
    return (lowcloudv, tropical, midlat), descr
end