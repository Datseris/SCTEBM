import Dates

function dynamical_system_progress_report(
        ebm::DynamicalSystem, text::String, filename, aux = nothing;
        folder = datadir("simulations", "autoprogressreport"),
    )
    preamble = """
    Report generated at: $(Dates.now())
    With git commit: $(DrWatson.gitdescribe())

    ---

    """
    footer = dynamical_system_summary(ebm)
    if !isnothing(aux)
        if any(a -> a isa Process, aux)
            aux = skipfirstline(sprint(show, MIME"text/plain"(), aux))
            aux = "\n\nProcesses used to create the system:\n"*aux
        else
            @assert aux isa AbstractString
        end
        footer *= aux
    end

    finaltext = preamble*text*"\n---\n\nDynamical system summary:\n\n"*footer

    dir = joinpath(folder, filename)
    write(dir, finaltext)
    return dir
end
