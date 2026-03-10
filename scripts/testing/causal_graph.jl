using DrWatson
@quickactivate
using ModelingToolkit
using ModelingToolkit: t_nounits as t
import Graphs

@variables x(t) = 1.0 y(t) = 1.0 z(t)

eqs = [
    Differential(t)(x) ~ -2x + z,
    Differential(t)(y) ~ -2y,
    z ~ x*y,
]

sys = ODESystem(eqs, t; name = :test)

ssys = structural_simplify(sys)

go = asgraph(sys; variables = [x, y, z])
gs = asgraph(ssys; variables = [x, y, z])

eqs = equations(sys)

graph = asdigraph(go, sys)

variables = [x, z, y]
varnames = string.(ModelingToolkit.SymbolicIndexingInterface.getname.(variables))
varvardep = varvar_dependencies(asgraph(sys; variables), variable_dependencies(sys; variables))

using CairoMakie, GraphMakie
fig = graphplot(varvardep; nlabels = varnames)
display(fig)

function get_diff_variables(sys)
    only.(get_variables.(getproperty.(diff_equations(sys), :lhs)))
end

diffvars = get_diff_variables(sys)
diffvarindices = [findfirst(isequal(v), variables) for v in diffvars]