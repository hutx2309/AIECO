# =============================================================================
# LowerBoundaryConductiveHeat.jl
#
# Process-level function for conductive heat exchange across the lower soil
# boundary.
#
# Design:
#   - pure function
#   - one lower boundary face at a time
#   - no mutation
#   - no @unpack
#   - no direct access to HFLWL
# =============================================================================


@doc raw"""
    lower_boundary_conductive_heat_process(...)

Calculate conductive heat exchange across the lower soil boundary.

This process corresponds to the WF9 lower-boundary heat branch:

```julia
if (N == 3) && (climateZone != -2)
```

The lower-boundary thermal conductance is

```math
G_b =
\frac{k_b}{z_b - z_s}
```

where `k_b` is the lower-boundary thermal conductivity, `z_b` is the
lower thermal boundary depth, and `z_s` is the depth of the current soil
layer.

The conductive heat flux is

```math
H_b =
G_b A (T_s - T_b) \Delta t
```

implemented by

```julia
conductive_heat_exchange(
    soil_temperature,
    lower_boundary_temperature,
    thermal_conductance,
    boundary_area,
    time_factor
)
```

# Sign convention

The returned `heat_flux` follows the sign convention of
`conductive_heat_exchange(...)`.

# Returns

A named tuple with activation flag, thermal conductance, and heat flux.
"""
function lower_boundary_conductive_heat_process(;
is_lower_boundary,
climate_zone,

soil_temperature,
lower_boundary_temperature,

lower_boundary_thermal_conductivity,
lower_boundary_depth,
soil_layer_depth,

boundary_area,
time_factor,

inactive_climate_zone = -2,
tiny = tiny_num2

)
# -------------------------------------------------------------------------
# 1. Activation gate
# -------------------------------------------------------------------------

is_active =
    is_lower_boundary &&
    climate_zone != inactive_climate_zone

if !is_active
    return (
        is_active = false,
        thermal_conductance = 0.0,
        heat_flux = 0.0
    )
end


# -------------------------------------------------------------------------
# 2. Boundary thermal conductance
# -------------------------------------------------------------------------

boundary_distance =
    lower_boundary_depth - soil_layer_depth

thermal_conductance =
    if boundary_distance > tiny
        lower_boundary_thermal_conductivity / boundary_distance
    else
        0.0
    end


# -------------------------------------------------------------------------
# 3. Conductive heat flux
# -------------------------------------------------------------------------

heat_flux =
    conductive_heat_exchange(
        soil_temperature,
        lower_boundary_temperature,
        thermal_conductance,
        boundary_area,
        time_factor
    )


# -------------------------------------------------------------------------
# 4. Return diagnostics
# -------------------------------------------------------------------------

return (
    is_active =
        true,

    thermal_conductance =
        thermal_conductance,

    heat_flux =
        heat_flux
)
end
