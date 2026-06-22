using Documenter

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const SOURCE_ROOT = joinpath(REPOSITORY_ROOT, "src")
const API_ROOT = joinpath(@__DIR__, "src", "api")

# Load the package directly from this checkout. This keeps local documentation
# builds independent of a docs/Manifest.toml file.
pushfirst!(LOAD_PATH, REPOSITORY_ROOT)
using AiECO

"""Fail when a source script has no matching page under docs/src/api/."""
function check_script_documentation()
    source_files = String[]
    for (directory, _, files) in walkdir(SOURCE_ROOT)
        for file in files
            endswith(file, ".jl") && push!(source_files, joinpath(directory, file))
        end
    end

    missing_pages = String[]
    for source_file in sort(source_files)
        relative_source = relpath(source_file, SOURCE_ROOT)
        relative_page = splitext(relative_source)[1] * ".md"
        page_file = joinpath(API_ROOT, relative_page)
        isfile(page_file) || push!(missing_pages, relative_page)
    end

    isempty(missing_pages) || error(
        "Every source script must have a documentation page. Missing:\n  " *
        join(missing_pages, "\n  "),
    )
end

check_script_documentation()

# The pure-equation files are still standalone scripts and are not yet included
# by AiECO. Load each one in an isolated documentation-only module so that its
# docstrings can be rendered without changing package wiring or causing name
# collisions between scripts.
module ScriptDocs
    module Thermodynamics
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "pure_equations", "Thermodynamics.jl"))
    end

    module Hydraulics
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "pure_equations", "Hydraulics.jl"))
    end

    module FluxLimiters
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "pure_equations", "FluxLimiters.jl"))
    end

    module ConductiveExchange
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "pure_equations", "ConductiveExchange.jl"))
    end

    module PhaseChange
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "pure_equations", "PhaseChange.jl"))
    end

    module Radiation
        const tiny_num = eps(Float64)
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "pure_equations", "Radiation.jl"))
    end

    module SurfaceExchange
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "pure_equations", "SurfaceExchange.jl"))
    end

    module BoundaryFluxes
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "pure_equations", "BoundaryFluxes.jl"))
    end

    module StorageUpdates
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "pure_equations", "StorageUpdates.jl"))
    end

end

const DOCUMENTED_MODULES = [
    AiECO,
    ScriptDocs.Thermodynamics,
    ScriptDocs.Hydraulics,
    ScriptDocs.FluxLimiters,
    ScriptDocs.ConductiveExchange,
    ScriptDocs.PhaseChange,
    ScriptDocs.Radiation,
    ScriptDocs.SurfaceExchange,
]

makedocs(
    sitename = "AiECO.jl",
    modules = DOCUMENTED_MODULES,
    clean = false, # OneDrive/Windows can keep handles open in docs/build.
    checkdocs = :exports,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        mathengine = Documenter.KaTeX(),
    ),
    pages = [
        "Home" => "index.md",
        "Documenting scripts" => "documenting_scripts.md",
        "Script reference" => [
            "AiECO package" => "api/AiECO.md",
            "Soil Water Energy" => [
                "Module entry point" => "api/Soil_Water_Energy/Soil_Water_Energy.md",
                "Pure equations" => [
                    "Thermodynamics" => "api/Soil_Water_Energy/pure_equations/Thermodynamics.md",
                    "Hydraulics" => "api/Soil_Water_Energy/pure_equations/Hydraulics.md",
                    "Flux limiters" => "api/Soil_Water_Energy/pure_equations/FluxLimiters.md",
                    "Conductive exchange" => "api/Soil_Water_Energy/pure_equations/ConductiveExchange.md",
                    "Phase change" => "api/Soil_Water_Energy/pure_equations/PhaseChange.md",
                    "Radiation" => "api/Soil_Water_Energy/pure_equations/Radiation.md",
                    "Surface exchange" => "api/Soil_Water_Energy/pure_equations/SurfaceExchange.md",
                    "Boundary Fluxes" => "api/Soil_Water_Energy/pure_equations/BoundaryFluxes.md",
                    "Storage Updates" => "api/Soil_Water_Energy/pure_equations/StorageUpdates.md"
                ],
            ],
        ],
    ],
)

deploydocs(
    repo = "github.com/hutx2309/AIECO.git",
    devbranch = "main",
)

