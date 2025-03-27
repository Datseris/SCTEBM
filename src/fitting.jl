import LsqFit

function nrmse(y, yfit)
    ymean = mean(y)
    mse = sum((yfit .- y).^2)
    nse = sum((ymean .- y).^2)
    return sqrt(mse/nse)
end

"""
    overplot_fits!(ax, x, y;
        fitexprs = [(x, p) -> @. p[1] + p[2]*x],
        fitstrs = ["a + b*X"],
    )

Fit and overplot the expressions in `fitexpr`, which is a vector of functions.
Each function is `y = f(x, p)` and is fitted via LsqFit.jl.
"""
function overplot_fits!(ax, x, y;
        fitexprs = [(x, p) -> @. p[1] + p[2]*x],
        fitstrs = ["a + b*X"],
        p0s = [rand(5) for _ in fitexprs], colorid = 3,
        legendkw = (position = :rb, ),
    )
    k = colorid
    ps = Vector{Vector{Float64}}(undef, length(fitexprs))
    errs = zeros(length(fitexprs))
    for (i, f) in enumerate(fitexprs)
        # perform fit
        fit = LsqFit.curve_fit(f, x, y, p0s[i])
        pfit = LsqFit.coef(fit) # best parameters
        xsort = range(minimum(x), maximum(x); length = 100)
        ysort = f(xsort, pfit)
        # plot fit while reporting error
        e = nrmse(y, f(x, pfit))
        lines!(ax, xsort, ysort;
            color = (COLORS[k], 0.9), linewidth = 4, linestyle = LINESTYLES[i],
            label = fitstrs[i]*", nrmse = $(string(round(e; sigdigits = 3)))"
        )
        # push fitted parameters to stored container
        ps[i] = pfit
        errs[i] = e
        k += 1
    end
    axislegend(ax; legendkw...)
    return ps, errs
end
