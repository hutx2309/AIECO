# =============================================================================
# BoundaryFluxes.jl
#
# Pure helper equations for surface runoff, grid-to-grid runoff, and
# water-table / external-boundary diagnostics in the water-energy module.
#
# Scope:
# - no mutation
# - no @unpack
# - no direct access to waterVar_copy
# - no array assignment
# - return fluxes, flags, or diagnostic values only
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Surface excess water / ice above litter holding capacity
# -----------------------------------------------------------------------------

@doc raw"""
partition_surface_excess_liquid_ice(
liquid_water,
ice,
holding_capacity;
tiny = tiny_num2
)

Partition excess surface water/ice into liquid and ice components while
preserving the liquid:ice ratio.

This represents the v2023 WF6 logic:

```math
\begin{aligned}
W_t &= W_l + I \\
W_{ret,l} &= \frac{W_l}{W_t}W_{hold} \\
W_{ret,i} &= \frac{I}{W_t}W_{hold} \\
W_{free,l} &= \max(0, W_l - W_{ret,l}) \\
W_{free,i} &= \max(0, I - W_{ret,i})
\end{aligned}
```

Returns:
excess_total
free_liquid
free_ice
"""
function partition_surface_excess_liquid_ice(
    liquid_water,
    ice,
    holding_capacity;
    tiny = tiny_num2
    )
    total_water_ice =
    liquid_water + ice

    if total_water_ice > tiny
        retained_liquid =
            liquid_water / total_water_ice * holding_capacity

        retained_ice =
            ice / total_water_ice * holding_capacity

        return (
            excess_total = max(0.0, total_water_ice - holding_capacity),
            free_liquid = max(0.0, liquid_water - retained_liquid),
            free_ice = max(0.0, ice - retained_ice)
        )
    else
        return (
            excess_total = 0.0,
            free_liquid = 0.0,
            free_ice = 0.0
        )
    end

end

# -----------------------------------------------------------------------------
# 2. Overland runoff from free surface liquid water
# -----------------------------------------------------------------------------

@doc raw"""
surface_overland_runoff(
free_liquid,
retention_capacity,
surface_area,
flow_area,
slope,
roughness_scale,
runoff_time_factor,
storage_removal_fraction,
available_liquid;
hydraulic_radius_divisor = 2.828
)

Calculate overland runoff from surface liquid water above the retention
capacity.

This captures the WF6 runoff pattern:

```math
\begin{aligned}
V_X &= W_{free,l} - W_{ret} \\
D &= \frac{V_X}{A_s} \\
R_h &= \frac{D}{2.828} \\
u &= \frac{R_h^{0.67}\sqrt{S}}{Z_M} \\
Q &= uD A_f (3.6\times10^3)\Delta t_r \\
Q_r &= \min(Q, V_X f_r, W_{avail,l}f_r)
\end{aligned}
```

Arguments:
free_liquid              : liquid water free to run off
retention_capacity       : surface retention capacity
surface_area             : area used to convert storage to water depth
flow_area                : grid-cell flow area used in runoff equation
slope                    : local slope
roughness_scale          : ZM in the current WF6 code
runoff_time_factor       : XNPHX in the current WF6 code
storage_removal_fraction : XNPXX in the current WF6 code
available_liquid         : liquid storage available for runoff removal

Returns:
runoff
velocity
excess_liquid
water_depth
"""
function surface_overland_runoff(
    free_surface_liquid,
    retention_capacity,
    surface_area,
    flow_area,
    slope,
    roughness_scale,
    runoff_time_factor,
    storage_removal_fraction,
    available_liquid;
    hydraulic_radius_divisor = 2.828
    )
    if free_surface_liquid <= retention_capacity
        return (
            runoff = 0.0,
            velocity = 0.0
        )
    end

    excess_liquid =
        free_surface_liquid - retention_capacity

    water_depth =
        excess_liquid / surface_area

    hydraulic_radius =
        water_depth / hydraulic_radius_divisor

    velocity =
        hydraulic_radius^0.67 *
        sqrt(max(0.0, slope)) /
        roughness_scale

    potential_runoff =
        velocity *
        water_depth *
        flow_area *
        3.6e3 *
        runoff_time_factor

    runoff =
        min(
            potential_runoff,
            excess_liquid * storage_removal_fraction,
            max(0.0, available_liquid) * storage_removal_fraction
        )

    return (
        runoff = max(0.0, runoff),
        velocity = velocity
    )
end
# -----------------------------------------------------------------------------
# 3. Wind redistribution of top snowpack material
# -----------------------------------------------------------------------------

@doc raw"""
wind_redistributed_top_snowpack(
wind_speed,
time_factor,
snow_volume,
liquid_water,
ice_volume
)

Calculate wind redistribution/removal of top snowpack material.

Legacy WF6 pattern:

```math
\begin{aligned}
Q_{SX} &= 10^{-7}u_{wind}\Delta t \\
Q_{SM} &= Q_{SX}S \\
Q_{WM} &= Q_{SX}W \\
Q_{IM} &= Q_{SX}I \\
Q_{ST} &= Q_{SM} + Q_{WM} + Q_{IM}
\end{aligned}
```

Returns:
snow_flux
water_flux
ice_flux
total_flux
"""
function wind_redistributed_top_snowpack(
    wind_speed,
    time_factor,
    snow_volume,
    liquid_water,
    ice_volume
    )
    redistribution_fraction =
    1.0e-7 * wind_speed * time_factor

    snow_flux =
        redistribution_fraction * snow_volume

    water_flux =
        redistribution_fraction * liquid_water

    ice_flux =
        redistribution_fraction * ice_volume

    return (
        snow_flux = snow_flux,
        water_flux = water_flux,
        ice_flux = ice_flux,
        total_flux = snow_flux + water_flux + ice_flux
    )

end

# -----------------------------------------------------------------------------
# 4. Runoff between adjacent grid cells
# -----------------------------------------------------------------------------

@doc raw"""
grid_to_grid_runoff_equilibrium(
source_elevation,
destination_elevation,
source_storage,
destination_storage,
source_area,
destination_area
)

Calculate equilibrium runoff amount from a source grid cell to a neighboring
destination grid cell.

This captures the WF7 expression:

```math
Q_{RQ} = \max\left(0,
\frac{(z_s-z_d)A_sA_d - S_sA_s + S_dA_d}{A_s + A_d}
\right)
```

where `S` is surface water/ice storage and `A` is cell area.

Returns:
nonnegative equilibrium runoff
"""
function grid_to_grid_runoff_equilibrium(
    source_elevation,
    destination_elevation,
    source_surface_storage,
    destination_surface_storage,
    source_area,
    destination_area
)
    if source_elevation <= destination_elevation
        return 0.0
    end

    return max(
        0.0,
        (
            (source_elevation - destination_elevation) *
            source_area *
            destination_area -
            source_surface_storage * source_area +
            destination_surface_storage * destination_area
        ) / (source_area + destination_area)
    )
end

# -----------------------------------------------------------------------------
# 5. Water-table and layer-fraction diagnostics
# -----------------------------------------------------------------------------

@doc raw"""
layer_fraction_below_depth(
layer_bottom_depth,
boundary_depth,
layer_thickness
)

Calculate the fraction of a soil layer below a given boundary depth.

Legacy pattern:

```math
f = \min\left(1, \max\left(0, \frac{z_{bottom} - z_b}{\Delta z}\right)\right)
```

Returns:
fraction in [0, 1]
"""
function layer_fraction_below_depth(
    layer_bottom_depth,
    boundary_depth,
    layer_thickness;
    tiny = tiny_num2
    )
    if layer_thickness > tiny
        return clamp(
            (layer_bottom_depth - boundary_depth) / layer_thickness,
            0.0,
            1.0
        )
    else
        return 0.0
    end
end

@doc raw"""
macropore_water_depth(
layer_bottom_depth,
macropore_water,
macropore_ice,
macropore_volume,
layer_thickness;
tiny = tiny_num2
)

Calculate the depth to the top of water/ice stored in the macropore domain.

Legacy WF9 pattern:

```math
z_M = z_{bottom} - \frac{W_M + I_M}{V_M}\Delta z
```

If macropore volume is unavailable, return layer_bottom_depth.

Returns:
depth
"""
function macropore_water_depth(
    layer_bottom_depth,
    macropore_water,
    macropore_ice,
    macropore_volume,
    layer_thickness;
    tiny = tiny_num2
    )
    if macropore_volume > tiny
        return layer_bottom_depth -
               (macropore_water + macropore_ice) /
               macropore_volume *
               layer_thickness
    else
        return layer_bottom_depth
    end
end

# -----------------------------------------------------------------------------
# 6. External-boundary surface runoff helpers
# -----------------------------------------------------------------------------

@doc raw"""
boundary_surface_runoff_flux(
excess_surface_water_ice,
area,
current_surface_elevation,
retention_capacity,
direction_sign,
slope_component,
slope_length,
natural_water_table_depth,
surface_depth,
available_runoff,
slope_fraction,
runoff_condition
)

Calculate external-boundary surface runoff flux.

This is a compact pure helper for the main WF9 boundary surface-runoff logic.
The sign follows the legacy convention: direction_sign controls the runoff
direction, so the returned runoff may be signed.

Returns:
runoff
"""
function boundary_surface_runoff_flux(
    excess_surface_water_ice,
    retention_capacity,
    area,
    current_surface_elevation,
    direction_sign,
    slope_component,
    slope_length,
    surface_depth,
    natural_water_table_depth,
    available_runoff,
    runoff_velocity,
    slope_fraction,
    runoff_condition,
    storage_removal_fraction
 )
    water_depth =
        excess_surface_water_ice / area

    retention_depth =
        retention_capacity / area

    source_elevation =
        current_surface_elevation + water_depth

    boundary_elevation =
        current_surface_elevation +
        retention_depth -
        direction_sign * slope_component * slope_length

    if source_elevation > boundary_elevation &&
       surface_depth - water_depth < natural_water_table_depth

        runoff =
            -direction_sign *
            available_runoff *
            slope_fraction *
            runoff_condition

        return (
            runoff = runoff,
            available_runoff = available_runoff,
            runoff_velocity = runoff_velocity,
            water_depth = water_depth,
            source_elevation = source_elevation,
            boundary_elevation = boundary_elevation
        )

    elseif surface_depth - water_depth > natural_water_table_depth

        water_table_adjustment =
            min(
                0.0,
                (
                    natural_water_table_depth -
                    surface_depth +
                    water_depth
                ) * area
            )

        adjusted_available_runoff =
            water_table_adjustment * storage_removal_fraction

        runoff =
            -direction_sign *
            adjusted_available_runoff *
            slope_fraction *
            runoff_condition

        return (
            runoff = runoff,
            available_runoff = adjusted_available_runoff,
            runoff_velocity = 0.0,
            water_depth = water_depth,
            source_elevation = source_elevation,
            boundary_elevation = boundary_elevation
        )

    else
        return (
            runoff = 0.0,
            available_runoff = available_runoff,
            runoff_velocity = runoff_velocity,
            water_depth = water_depth,
            source_elevation = source_elevation,
            boundary_elevation = boundary_elevation
        )
    end
end

# -----------------------------------------------------------------------------
# 7. Boundary-condition flags
# -----------------------------------------------------------------------------

@doc raw"""
external_surface_runoff_allowed(
runoff_allowed_flag,
runoff_condition,
surface_depth,
surface_layer_thickness,
initial_elevation,
lower_boundary_bulk_density;
tiny = tiny_num
)

Check whether external-boundary surface runoff is allowed.

This helper separates the Boolean gate from the flux equation. It should be
used only if it exactly matches the relevant WF9 branch.
"""
function external_surface_runoff_is_allowed(
    runoff_allowed_flag,
    runoff_condition,
    surface_depth,
    surface_layer_thickness,
    initial_elevation,
    lower_boundary_bulk_density;
    tiny = tiny_num
    )
    if runoff_allowed_flag == 0
        return false
    end

    if runoff_condition == 0
        return false
    end

    if surface_depth - surface_layer_thickness > initial_elevation &&
       lower_boundary_bulk_density <= tiny
        return false
    end

    return true
end

# -----------------------------------------------------------------------------
# 8. Water-table exchange helpers
# -----------------------------------------------------------------------------

@doc raw"""
hydrostatic_pressure_potential(
layer_depth,
water_table_depth;
hydraulic_gradient = 0.0098
)

Calculate hydrostatic pressure potential relative to a water-table depth.

The default hydraulic gradient is:

```math
g_h = 0.0098\ \mathrm{MPa\ m^{-1}}
```

which corresponds approximately to:

```math
\rho_w g = 1000\ \mathrm{kg\ m^{-3}}\times 9.8\ \mathrm{m\ s^{-2}}
= 9800\ \mathrm{Pa\ m^{-1}}
= 0.0098\ \mathrm{MPa\ m^{-1}}
```

Legacy WF9 pattern:

```math
\psi_h = 0.0098(z_L - z_{wt})
```
Returns:
pressure potential in the same units as soil water potential.
"""
function hydrostatic_pressure_potential(
    layer_depth,
    water_table_depth;
    hydraulic_gradient = 0.0098
    )
    return hydraulic_gradient * (layer_depth - water_table_depth)
end

@doc raw"""
adjusted_water_table_depth(
water_table_depth,
saturation_potential;
hydraulic_gradient = 0.0098
)

Calculate an adjusted/effective water-table depth from saturated water
potential.

Legacy WF9 pattern:

```math
z_{wt,x} = z_{wt} + \frac{\psi_{sat}}{0.0098}
```

Returns:
adjusted water-table depth
"""
function adjusted_water_table_depth(
    water_table_depth,
    saturation_potential;
    hydraulic_gradient = 0.0098
    )
    return water_table_depth + saturation_potential / hydraulic_gradient
end

@doc raw"""
slope_pressure_adjustment(
direction_sign,
slope_component,
flow_length,
water_table_slope;
coefficient = 0.0049
)

Calculate slope-related pressure adjustment for lateral water-table exchange.

```math
\Delta\psi_s = s_d c_s S L(1-S_{wt})
```

`c_s` is 0.005 for micropore discharge above natural/tile water table and
0.0049 for macropore discharge and recharge branches.

Returns:
slope pressure adjustment
"""
function slope_pressure_adjustment(
    direction_sign,
    slope_component,
    flow_length,
    water_table_slope;
    coefficient = 0.0049
    )
    return direction_sign *
           coefficient *
           slope_component *
           flow_length *
           (1.0 - water_table_slope)
end

@doc raw"""
water_table_discharge_potential(
matric_potential,
osmotic_potential,
source_depth,
water_table_depth,
drainage_reference_depth,
slope_adjustment;
include_matric = true,
hydraulic_gradient = 0.0098,
osmotic_factor = 0.03
)

Calculate pressure potential for discharge from soil to water-table or
tile-drain boundary.

This helper represents branches where the potential is constrained to be
nonpositive:

```math
\psi_d = \min(0, \psi_b)
```

For macropore branches, use `include_matric = false` because the legacy
macropore expression usually omits `-PSISA1[L]`.

Arguments:
matric_potential          : soil matric potential, e.g., PSISA1[L]
osmotic_potential         : soil osmotic potential
source_depth              : soil midpoint depth or macropore-water depth
water_table_depth         : natural or artificial water-table depth
drainage_reference_depth  : reference drain depth, e.g., waterTbl_DepZ
slope_adjustment          : slope pressure adjustment, e.g., PSISWD

Returns:
nonpositive discharge potential
"""
function water_table_discharge_potential(
    matric_potential,
    osmotic_potential,
    source_depth,
    water_table_depth,
    drainage_reference_depth,
    slope_adjustment;
    include_matric = true,
    clamp_nonpositive = true,
    hydraulic_gradient = 0.0098,
    osmotic_factor = 0.03
    )
    base_potential =
        -osmotic_factor * osmotic_potential +
        hydrostatic_pressure_potential(
            source_depth,
            water_table_depth;
            hydraulic_gradient = hydraulic_gradient
        ) -
        hydraulic_gradient *
        max(0.0, source_depth - drainage_reference_depth)

    if include_matric
        base_potential -= matric_potential
    end

    discharge_potential =
        clamp_nonpositive ? min(0.0, base_potential) : base_potential

    if discharge_potential < 0.0
        discharge_potential -= slope_adjustment
    end

    return discharge_potential
end

@doc raw"""
water_table_recharge_potential(
matric_potential,
osmotic_potential,
source_depth,
water_table_depth,
slope_adjustment;
include_matric = true,
hydraulic_gradient = 0.0098,
osmotic_factor = 0.03
)

Calculate pressure potential for recharge from water table into soil.

This helper represents branches where the potential is constrained to be
nonnegative:

```math
\psi_r = \max(0, \psi_b)
```

For macropore branches, use `include_matric = false` because the legacy
macropore expression usually omits `-PSISA1[L]`.

Returns:
nonnegative recharge potential
"""
function water_table_recharge_potential(
    matric_potential,
    osmotic_potential,
    source_depth,
    water_table_depth,
    slope_adjustment;
    include_matric = true,
    clamp_nonnegative = true,
    hydraulic_gradient = 0.0098,
    osmotic_factor = 0.03
    )
    base_potential =
        -osmotic_factor * osmotic_potential +
        hydrostatic_pressure_potential(
            source_depth,
            water_table_depth;
            hydraulic_gradient = hydraulic_gradient
        )

    if include_matric
        base_potential -= matric_potential
    end

    recharge_potential =
        clamp_nonnegative ? max(0.0, base_potential) : base_potential

    if recharge_potential > 0.0
        recharge_potential += slope_adjustment
    end

    return recharge_potential
end

@doc raw"""
boundary_water_flux_from_potential(
pressure_potential,
hydraulic_conductivity,
area,
active_fraction,
resistance_factor,
boundary_activity_factor,
time_factor
)

Convert a pressure potential into a boundary water flux.

```math
F_w = \frac{\psi_b K A f}{R + 1}\alpha_b\Delta t
```

This helper can be used for:
- micropore discharge to natural water table
- micropore recharge from natural water table
- macropore discharge to natural water table
- macropore recharge from natural water table
- artificial water table / tile-drain analogs

Returns:
signed water flux
"""
function boundary_water_flux_from_potential(
    pressure_potential,
    hydraulic_conductivity,
    area,
    active_fraction,
    resistance_factor,
    boundary_activity_factor,
    time_factor
)
    return pressure_potential *
           hydraulic_conductivity *
           area *
           active_fraction /
           (resistance_factor + 1.0) *
           boundary_activity_factor *
           time_factor
end

