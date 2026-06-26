# =============================================================================
# BoundaryUnsaturatedSubsurfaceFlow.jl
#
# Process-level function for unsaturated subsurface liquid-water flow across
# an external or lower boundary when explicit water-table exchange is inactive.
#
# Design:
#   - pure function
#   - one soil layer and one boundary face at a time
#   - no mutation
#   - no @unpack
#   - no direct access to FLWL / FLWHL / HFLWL arrays
# =============================================================================


@doc raw"""
    boundary_unsaturated_subsurface_flow_process(...)

Calculate unsaturated subsurface water flow across one boundary face when
explicit water-table exchange is inactive.

This process corresponds to the legacy WF9 branch

```julia
if (waterTbl_Flag == 0) || (N == 3)
```

It calculates:

1. micropore liquid-water boundary flow,
2. macropore liquid-water boundary flow,
3. advective heat carried by the combined boundary water flow.

The bounded micropore volumetric water content is

```math
\theta =
\min\left[
    \phi,
    \max\left(\theta_{\min}, \frac{W_\mu}{V_s}\right)
\right]
```

where `W_\mu` is micropore liquid water, `V_s` is the soil reference
volume, `\phi` is porosity, and `\theta_{\min}` is the minimum water
content used for hydraulic lookup.

The hydraulic-conductivity table index is

```math
k = I_K(\theta, \phi)
```

implemented by

```julia
conductivity_table_index(theta, porosity)
```

The boundary driving term is

```math
D_b =
s_d g_h \left(-|S|\right)
```

where `s_d` is the boundary direction sign, `g_h` is the hydrostatic
pressure gradient, and `S` is the slope component.

The raw micropore flux before boundary factors is

```math
F_\mu^* =
D_b K_\mu A
```

and it is first limited by the current removable micropore water storage:

```math
\tilde{F}_\mu =
\min\left[
    W_\mu f_r,
    \max\left(-W_\mu f_r, F_\mu^*\right)
\right]
```

The final micropore boundary flux is

```math
F_\mu =
\tilde{F}_\mu f_b f_a \Delta t
```

where `f_b` is the boundary flow factor, `f_a` is the boundary activity
factor, and `\Delta t` is the process time factor.

Macropore boundary flow follows the same structure:

```math
F_M^* =
D_b K_M A
```

```math
\tilde{F}_M =
\min\left[
    W_M f_r,
    \max\left(-W_M f_r, F_M^*\right)
\right]
```

```math
F_M =
\tilde{F}_M f_b f_a \Delta t
```

The advective heat flux is calculated from the combined liquid-water flux:

```math
H =
c_w T_s (F_\mu + F_M)
```

implemented by

```julia
advective_heat_water(micropore_flux + macropore_flux, soil_temperature, cpw)
```

# Sign convention

The returned fluxes use the same sign convention as legacy WF9:

* positive or negative sign is controlled by `direction_sign`
* `micropore_flux` maps to `FLWL[...]`
* `micropore_mobile_flux` maps to `FLWLX[...]`
* `macropore_flux` maps to `FLWHL[...]`

# Returns

A named tuple with water contents, conductivity index, conductivities,
limited micropore/macropore boundary fluxes, and advective heat flux.
"""
function boundary_unsaturated_subsurface_flow_process(;
# -------------------------------------------------------------------------
# Boundary geometry and direction
# -------------------------------------------------------------------------
direction_sign,
slope_component,
boundary_area,

# -------------------------------------------------------------------------
# Soil water state
# -------------------------------------------------------------------------
micropore_water,
macropore_water,

soil_reference_volume,
porosity,
minimum_water_content,

# -------------------------------------------------------------------------
# Hydraulic conductivity
# -------------------------------------------------------------------------
hydraulic_conductivity_table,
boundary_axis_index,
layer_index,
micropore_conductivity_multiplier,
macropore_hydraulic_conductivity,

# -------------------------------------------------------------------------
# Boundary factors and time
# -------------------------------------------------------------------------
boundary_flow_factor,
boundary_activity_factor,
storage_removal_fraction,
process_time_factor,

# -------------------------------------------------------------------------
# Heat
# -------------------------------------------------------------------------
soil_temperature,
cpw,

# -------------------------------------------------------------------------
# Constants / numerical controls
# -------------------------------------------------------------------------
hydraulic_gradient = 0.0098
)
# -------------------------------------------------------------------------
# 1. Micropore water content and hydraulic conductivity
# -------------------------------------------------------------------------

micropore_water_content =
    bounded_volumetric_water_content(
        micropore_water,
        soil_reference_volume,
        minimum_water_content,
        porosity
    )

conductivity_index =
    conductivity_table_index(
        micropore_water_content,
        porosity
    )

micropore_hydraulic_conductivity =
    hydraulic_conductivity_table[
        boundary_axis_index,
        conductivity_index,
        layer_index
    ] *
    micropore_conductivity_multiplier


# -------------------------------------------------------------------------
# 2. Common boundary driving term
#
# Preserve the legacy WF9 expression:
#
#     XN * 0.0098 * (-abs(SLOPE[N+1]))
#
# -------------------------------------------------------------------------

boundary_drive =
    direction_sign *
    hydraulic_gradient *
    (-abs(slope_component))


# -------------------------------------------------------------------------
# 3. Micropore boundary water flux
# -------------------------------------------------------------------------

removable_micropore_water =
    max(
        0.0,
        micropore_water * storage_removal_fraction
    )

raw_micropore_flux =
    boundary_drive *
    micropore_hydraulic_conductivity *
    boundary_area

storage_limited_micropore_flux =
    min(
        removable_micropore_water,
        max(
            -removable_micropore_water,
            raw_micropore_flux
        )
    )

micropore_flux =
    storage_limited_micropore_flux *
    boundary_flow_factor *
    boundary_activity_factor *
    process_time_factor


# -------------------------------------------------------------------------
# 4. Macropore boundary water flux
# -------------------------------------------------------------------------

removable_macropore_water =
    max(
        0.0,
        macropore_water * storage_removal_fraction
    )

raw_macropore_flux =
    boundary_drive *
    macropore_hydraulic_conductivity *
    boundary_area

storage_limited_macropore_flux =
    min(
        removable_macropore_water,
        max(
            -removable_macropore_water,
            raw_macropore_flux
        )
    )

macropore_flux =
    storage_limited_macropore_flux *
    boundary_flow_factor *
    boundary_activity_factor *
    process_time_factor


# -------------------------------------------------------------------------
# 5. Heat carried by combined boundary liquid-water flux
# -------------------------------------------------------------------------

heat_flux =
    advective_heat_water(
        micropore_flux + macropore_flux,
        soil_temperature,
        cpw
    )


# -------------------------------------------------------------------------
# 6. Return process diagnostics
# -------------------------------------------------------------------------

return (
    micropore_water_content =
        micropore_water_content,

    conductivity_index =
        conductivity_index,

    micropore_hydraulic_conductivity =
        micropore_hydraulic_conductivity,

    macropore_hydraulic_conductivity =
        macropore_hydraulic_conductivity,

    boundary_drive =
        boundary_drive,

    removable_micropore_water =
        removable_micropore_water,

    removable_macropore_water =
        removable_macropore_water,

    raw_micropore_flux =
        raw_micropore_flux,

    raw_macropore_flux =
        raw_macropore_flux,

    storage_limited_micropore_flux =
        storage_limited_micropore_flux,

    storage_limited_macropore_flux =
        storage_limited_macropore_flux,

    micropore_flux =
        micropore_flux,

    micropore_mobile_flux =
        micropore_flux,

    macropore_flux =
        macropore_flux,

    heat_flux =
        heat_flux
)

end
