# AiECO.jl

AiECO is a modular ecosystem modeling platform combining process-based
conservation laws with auditable AI components. The platform is intended to
support independently testable soil, plant, microbe, weather, data, and AI
modules.

This site documents the implementation one Julia script at a time. Each source
file under [`src/`](https://github.com/hutx2309/AIECO/tree/main/src) has a
matching page in the **Script reference** section.

## Current implementation

- `AiECO.jl` defines the package entry point.
- `Soil_Water_Energy.jl` is the developing soil water-energy module entry point.
- `Soil_Water_Energy/pure_equations/` contains non-mutating equations for
  thermodynamics, hydraulics, flux limiting, conductive exchange, phase change,
  radiation, and surface exchange.
- `Soil_Water_Energy/process_functions/` contains process-level functions that
  compose the lower-level equations into boundary, phase-change, runoff, and
  pore-domain exchange calculations.

See [Documenting scripts](@ref) for the required page structure and build
commands.
