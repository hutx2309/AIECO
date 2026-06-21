# Soil Water Energy module entry point

Source: [`src/Soil_Water_Energy/Soil_Water_Energy.jl`](https://github.com/hutx2309/AIECO/blob/main/src/Soil_Water_Energy/Soil_Water_Energy.jl)

This script is the developing entry point for the coupled soil water-energy
module. It currently establishes the module boundary; the pure-equation scripts
listed below have not yet been included into it.

## Integration checklist

- Include each completed component in dependency order.
- Define the public exports explicitly.
- Keep pure equations independent of mutable model state.
- Add process/state orchestration only after the pure functions are tested.
- Add the module to `src/AiECO.jl` when its public interface is stable.

