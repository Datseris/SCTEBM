using DrWatson
@quickactivate
using PkgCite
PkgCite.get_tool_citation(
    only_direct = true,
    filename = papersdir("pkgbib.tex"),
)