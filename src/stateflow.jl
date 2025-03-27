# function nullclines_observable!(ax::Axis, dsobs::DynamicalSystemObservable; kw...)
#     ds = dsobs.ds # find the correct non-parallel ds version
#     xlim, ylim = physically_plausible_limits(dsobs.ds) # get correct reference
#     ux, uy = nullclines!(ax, ds, xlim, ylim; kw...)
#     on(dsobs.param_observable) do params
#         nullclines!(ux[], uy[], ds, limits[1], limits[2])
#         notify.((ux, uy))
#     end
#     return ux, uy
# end

"""
    nullclines!(ax::Axis, ds::DynamicalSystem, xlim, ylim;
        n = 1000, nx = n, ny = n, xkw = NamedTuple(), ykw = NamedTuple(),
    )

Plot nullclines.
"""
function nullclines!(ax::Axis, ds::DynamicalSystem, xlim, ylim;
        n = 1000, nx = n, ny = n, xkw = NamedTuple(), ykw = NamedTuple(),
    )
    ux = Observable(zeros(nx, ny))
    uy = Observable(zeros(nx, ny))
    xr = range(xlim[1], xlim[2]; length = nx)
    yr = range(ylim[1], ylim[2]; length = ny)
    nullclines!(ux[], uy[], ds, xlim, ylim)
    contour!(ax, xr, yr, ux; levels = [0.0], color = :black, linewidth = 2.0, linestyle = :dash, xkw...)
    contour!(ax, xr, yr, uy; levels = [0.0], color = :grey, linewidth = 2.0, linestyle = :dash, ykw...)
    return ux, uy
end

function nullclines!(ux::AbstractArray, uy::AbstractArray, ds::DynamicalSystem, xlim, ylim)
    dimension(ds) == 2 || error("only works with 2D systems")
    nx, ny = size(ux)
    xr = range(xlim[1], xlim[2]; length = nx)
    yr = range(ylim[1], ylim[2]; length = ny)
    du = zeros(2)
    for (j, y) in enumerate(yr)
        for (i, x) in enumerate(xr)
            u = SVector(x, y)
            if isinplace(ds)
                dynamic_rule(ds)(du, u, current_parameters(ds), current_time(ds))
            else
                du = dynamic_rule(ds)(u, current_parameters(ds), current_time(ds))
            end
            ux[i, j] = du[1]
            uy[i, j] = du[2]
        end
    end
    return ux, uy
end
