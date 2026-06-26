# =============================================================================
# ExternalBoundarySurfaceRunoff.jl
#
# Process-level function for surface runoff across an external grid boundary.
#
# Design:
#   - pure function
#   - one boundary face at a time
#   - no mutation
#   - no @unpack
#   - no direct access to QR1, HQR1, QR, HQR, QRMN, or IFLBM arrays
# =============================================================================


@doc raw"""
    external_boundary_surface_runoff_process(...)

Calculate surface water runoff across one external horizontal boundary face.

This process combines:

1. Surface-runoff activation/gating.
2. Boundary runoff calculation.
3. Advective heat carried by runoff water.
4. Updated available runoff state.

The surface runoff depth above free-water retention is

```math
d_w =
\frac{W_{excess}}{A}
```

where `W_{excess}` is excess surface water/ice storage and `A` is
surface area.

The downstream reference depth is

```math
d_r =
\frac{W_{ret}}{A}
```

where `W_{ret}` is ground water retention capacity.

The water-surface elevations are

```math
z_1 = z_s + d_w
```

and

```math
z_2 = z_s + d_r - s_d S L
```

where `z_s` is current surface elevation, `s_d` is boundary direction,
`S` is slope, and `L` is boundary flow length.

When the surface-water elevation exceeds the boundary reference elevation,
runoff is

```math
Q_r = -s_d Q_{avail} f_s f_b
```

where `Q_{avail}` is available runoff, `f_s` is the slope/aspect
runoff fraction, and `f_b` is the boundary runoff condition factor.

If the free surface intersects the water table, the available runoff can be
reset by the water-table storage constraint:

```math
Q_{avail} =
\min\left[
    0,
    \left(z_{wt} - z_{surf} + d_w\right) A
\right] f_r
```

The advective heat flux is

```math
H_r = c_w T_s Q_r
```

implemented by `advective_heat_water`.

# Sign convention

The returned `runoff` uses the same sign convention as legacy `QR1`:

* positive or negative sign is determined by `direction_sign`
* `runoff = 0.0` means no boundary surface runoff

# Returns

A named tuple with runoff flux, heat flux, updated available runoff,
updated runoff velocity, and activation diagnostics.
"""
function external_boundary_surface_runoff_process(;
# Activation / boundary context
is_top_layer,
is_horizontal_boundary,
runoff_allowed_flag,
runoff_condition_factor,


# Surface and soil state
surface_depth,
surface_reference_depth,
initial_surface_elevation,
topsoil_bulk_density,

excess_surface_water_ice,
ground_water_retention_capacity,
current_surface_elevation,
natural_water_table_depth,

# Geometry and direction
direction_sign,
slope_component,
boundary_flow_length,
surface_area,
slope_runoff_fraction,

# Available runoff state
available_runoff,
runoff_velocity,

# Time and heat
storage_removal_fraction,
surface_temperature,
cpw
)
# -------------------------------------------------------------------------
# 1. Gate: only horizontal external boundaries at the surface layer
# -------------------------------------------------------------------------

if !(is_top_layer && is_horizontal_boundary)
    return (
        runoff = 0.0,
        heat_flux = 0.0,
        available_runoff = available_runoff,
        runoff_velocity = runoff_velocity,
        boundary_is_active = false,
        runoff_is_allowed = false
    )
end


# -------------------------------------------------------------------------
# 2. Gate: boundary runoff permission
# -------------------------------------------------------------------------

runoff_is_allowed =
    external_surface_runoff_is_allowed(
        runoff_allowed_flag,
        runoff_condition_factor,
        surface_depth,
        surface_reference_depth,
        initial_surface_elevation,
        topsoil_bulk_density
    )

if !runoff_is_allowed
    return (
        runoff = 0.0,
        heat_flux = 0.0,
        available_runoff = available_runoff,
        runoff_velocity = runoff_velocity,
        boundary_is_active = true,
        runoff_is_allowed = false
    )
end


# -------------------------------------------------------------------------
# 3. Boundary runoff calculation
# -------------------------------------------------------------------------

boundary_runoff =
    boundary_surface_runoff_flux(
        excess_surface_water_ice,
        ground_water_retention_capacity,
        surface_area,
        current_surface_elevation,
        direction_sign,
        slope_component,
        boundary_flow_length,
        surface_depth,
        natural_water_table_depth,
        available_runoff,
        runoff_velocity,
        slope_runoff_fraction,
        runoff_condition_factor,
        storage_removal_fraction
    )

runoff =
    boundary_runoff.runoff

heat_flux =
    advective_heat_water(
        runoff,
        surface_temperature,
        cpw
    )


# -------------------------------------------------------------------------
# 4. Return process result
# -------------------------------------------------------------------------

return (
    runoff =
        runoff,

    heat_flux =
        heat_flux,

    available_runoff =
        boundary_runoff.available_runoff,

    runoff_velocity =
        boundary_runoff.runoff_velocity,

    boundary_is_active =
        true,

    runoff_is_allowed =
        true
)
end
