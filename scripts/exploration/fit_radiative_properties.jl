# For fitting cloud albedo based on expressions of old papers
using DrWatson
@quickactivate
include(srcdir("theme.jl"))

fig, axs = axesgrid(2, 1; xlabels = "LWP (g/m²)", ylabels = ["τ_C", "α_C"], sharex = true)

lwp = 10 .^ range(0, 3; length = 101)
tau = @. lwp^1.7/(10 + lwp)
alpha = @. 0.13tau/(1 + 0.13tau)
lines!(axs[1], lwp, tau)
lines!(axs[2], lwp, alpha)
axs[1].xscale = log10
axs[1].yscale = log10
axs[2].xscale = log10

fig