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

# The Soil_Water_Energy component files are still standalone scripts and are not
# yet included by AiECO. Load each one in an isolated documentation-only module
# so that its docstrings can be rendered without changing package wiring or
# causing name collisions between scripts.
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

    module ExternalBoundarySurfaceRunoff
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "process_functions", "ExternalBoundarySurfaceRunoff.jl"))
    end

    module BoundaryUnsaturatedSubsurfaceFlow
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "process_functions", "BoundaryUnsaturatedSubsurfaceFlow.jl"))
    end

    module LitterSoilCapillaryExchange
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "process_functions", "LitterSoilCapillaryExchange.jl"))
    end

    module LowerBoundaryConductiveHeat
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "process_functions", "LowerBoundaryConductiveHeat.jl"))
    end

    module SoilLayerPhaseChange
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "process_functions", "SoilLayerPhaseChange.jl"))
    end

    module SoilInterlayerConductiveHeat
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "process_functions", "SoilInterlayerConductiveHeat.jl"))
    end

    module SoilInterlayerMacroporeFlow
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "process_functions", "SoilInterlayerMacroporeFlow.jl"))
    end

    module SoilInterlayerMicroporeFlow
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "process_functions", "SoilInterlayerMicroporeFlow.jl"))
    end

    module SoilInterlayerVaporDiffusion
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "process_functions", "SoilInterlayerVaporDiffusion.jl"))
    end

    module WaterTableBoundaryProcess
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "process_functions", "WaterTableBoundaryProcess.jl"))
    end

    module SoilPoreDomainExchange
        include(joinpath(@__DIR__, "..", "src", "Soil_Water_Energy", "process_functions", "SoilPoreDomainExchange.jl"))
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
    ScriptDocs.BoundaryFluxes,
    ScriptDocs.StorageUpdates,
    ScriptDocs.ExternalBoundarySurfaceRunoff,
    ScriptDocs.BoundaryUnsaturatedSubsurfaceFlow,
    ScriptDocs.LitterSoilCapillaryExchange,
    ScriptDocs.LowerBoundaryConductiveHeat,
    ScriptDocs.SoilLayerPhaseChange,
    ScriptDocs.SoilInterlayerConductiveHeat,
    ScriptDocs.SoilInterlayerMacroporeFlow,
    ScriptDocs.SoilInterlayerMicroporeFlow,
    ScriptDocs.SoilInterlayerVaporDiffusion,
    ScriptDocs.WaterTableBoundaryProcess,
    ScriptDocs.SoilPoreDomainExchange,
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
                "Process functions" => [
                    "External boundary surface runoff" => "api/Soil_Water_Energy/process_functions/ExternalBoundarySurfaceRunoff.md",
                    "Boundary unsaturated subsurface flow" => "api/Soil_Water_Energy/process_functions/BoundaryUnsaturatedSubsurfaceFlow.md",
                    "Litter-soil capillary exchange" => "api/Soil_Water_Energy/process_functions/LitterSoilCapillaryExchange.md",
                    "Lower boundary conductive heat" => "api/Soil_Water_Energy/process_functions/LowerBoundaryConductiveHeat.md",
                    "Soil interlayer conductive heat" => "api/Soil_Water_Energy/process_functions/SoilInterlayerConductiveHeat.md",
                    "Soil interlayer macropore flow" => "api/Soil_Water_Energy/process_functions/SoilInterlayerMacroporeFlow.md",
                    "Soil interlayer micropore flow" => "api/Soil_Water_Energy/process_functions/SoilInterlayerMicroporeFlow.md",
                    "Soil interlayer vapor diffusion" => "api/Soil_Water_Energy/process_functions/SoilInterlayerVaporDiffusion.md",
                    "Soil layer phase change" => "api/Soil_Water_Energy/process_functions/SoilLayerPhaseChange.md",
                    "Water table boundary process" => "api/Soil_Water_Energy/process_functions/WaterTableBoundaryProcess.md",
                    "Soil pore domain exchange" => "api/Soil_Water_Energy/process_functions/SoilPoreDomainExchange.md"
                ],
            ],
        ],
    ],
)

deploydocs(
    repo = "github.com/hutx2309/AIECO.git",
    devbranch = "main",
)
