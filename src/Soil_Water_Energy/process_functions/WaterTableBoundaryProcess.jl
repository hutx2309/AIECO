# =============================================================================
# WaterTableBoundaryProcess.jl
#
# Process-level functions for soil water exchange with natural/artificial
# water-table boundaries and tile-drain boundaries.
#
# Design:
#   - pure functions only
#   - no mutation
#   - no @unpack
#   - no direct access to waterVar_copy
#   - process functions call lower-level pure equations when useful
# =============================================================================
"""
    water_table_layer_boundary_exchange(; kwargs...)

Calculate water exchange between one soil layer and one water-table boundary.

This process function represents the physical sequence

```math
\\Delta \\psi_s =
s_d c_s S L (1 - S_{wt})
```

where `\\Delta \\psi_s` is the slope-induced pressure adjustment,
`s_d` is the directional sign, `c_s` is a slope-pressure coefficient,
`S` is the local slope component, `L` is the lateral flow length, and
`S_{wt}` is the water-table slope.

For discharge from soil to the boundary, the pressure potential is

```math
\\psi_d =
\\min\\left[
0,
-I_m \\psi_m
-c_o \\psi_o
+ g_h (z - z_{wt})
- g_h \\max(0, z - z_d)
\\right]
```

and the slope-adjusted discharge potential is

```math
\\psi_d^* =
\\begin{cases}
\\psi_d - \\Delta \\psi_s, & \\psi_d < 0 \\\\
\\psi_d, & \\psi_d \\ge 0
\\end{cases}
```

For recharge from the water table into soil, the pressure potential is

```math
\\psi_r =
\\max\\left[
0,
-I_m \\psi_m
-c_o \\psi_o
+ g_h (z - z_{wt})
\\right]
```

and the slope-adjusted recharge potential is

```math
\\psi_r^* =
\\begin{cases}
\\psi_r + \\Delta \\psi_s, & \\psi_r > 0 \\\\
\\psi_r, & \\psi_r \\le 0
\\end{cases}
```

The boundary water flux is

```math
F_w =
\\frac{
    \\psi_b K A f
}{
    R + 1
}
\\alpha_b \\Delta t
```

where `\\psi_b` is the boundary pressure potential, `K` is hydraulic
conductivity or conductance, `A` is boundary area, `f` is the active
layer fraction, `R` is a resistance factor, `\\alpha_b` is the boundary
activity factor, and `\\Delta t` is the process time step.

# Arguments
## Geometry and boundary

* `direction_sign`: directional sign, e.g. `XN`
* `slope_component`: terrain slope component, e.g. `SLOPE[N+1]`
* `flow_length`: lateral flow length, e.g. `soil_cube[N, L]`
* `water_table_slope`: water-table slope factor
* `layer_mid_depth`: soil layer midpoint depth
* `macropore_water_depth`: depth to macropore water/ice surface
* `water_table_depth`: natural or artificial water-table depth
* `drainage_reference_depth`: drainage reference depth, e.g. `waterTbl_DepZ`
* `boundary_area`: lateral exchange area

## Potentials

* `matric_potential`: soil matric potential
* `osmotic_potential`: soil osmotic potential

## Conductivities

* `micropore_conductivity`: micropore hydraulic conductivity
* `macropore_conductivity`: macropore hydraulic conductivity or conductance

## Active fractions

* `fraction_above_water_table`: layer fraction above the water table
* `fraction_below_water_table`: layer fraction below the water table

## Resistance/activity

* `micropore_resistance_factor`: resistance factor for micropore exchange
* `macropore_resistance_factor`: resistance factor for macropore exchange
* `boundary_activity_factor`: boundary activity/recharge factor
* `time_factor`: process time factor

## Switches

* `do_micropore_discharge`
* `do_macropore_discharge`
* `do_micropore_recharge`
* `do_macropore_recharge`

# Returns

A named tuple containing pressure potentials, water fluxes, and total water
fluxes.
"""
function water_table_layer_boundary_exchange(;
direction_sign,
slope_component,
flow_length,
water_table_slope,


layer_mid_depth,
macropore_water_depth,
water_table_depth,
drainage_reference_depth,
boundary_area,

matric_potential,
osmotic_potential,

micropore_conductivity,
macropore_conductivity,

fraction_above_water_table,
fraction_below_water_table,

micropore_resistance_factor,
macropore_resistance_factor,
boundary_activity_factor,
time_factor,

do_micropore_discharge = true,
do_macropore_discharge = true,
do_micropore_recharge = true,
do_macropore_recharge = true,

micropore_discharge_slope_coefficient = 0.0049,   # CLAUDE DEBUG: was 0.005 — typo introduced in 2026-06-21 refactor; Fortran uses 0.0049 uniformly at watsub.f90:4036/4077/4115/4153/4192/4236. 2% error compounded explains the observed WTR_TBL 6-9% deeper drift from D=250.
other_slope_coefficient = 0.0049,

hydraulic_gradient = 0.0098,
osmotic_factor = 0.03

)
# -------------------------------------------------------------------------
# 1. Slope pressure adjustments
# -------------------------------------------------------------------------

slope_adjustment_micropore_discharge =
    slope_pressure_adjustment(
        direction_sign,
        slope_component,
        flow_length,
        water_table_slope;
        coefficient = micropore_discharge_slope_coefficient
    )

slope_adjustment_other =
    slope_pressure_adjustment(
        direction_sign,
        slope_component,
        flow_length,
        water_table_slope;
        coefficient = other_slope_coefficient
    )

# -------------------------------------------------------------------------
# 2. Discharge potentials: soil -> water-table boundary
# -------------------------------------------------------------------------

micropore_discharge_potential =
    do_micropore_discharge ?
    water_table_discharge_potential(
        matric_potential,
        osmotic_potential,
        layer_mid_depth,
        water_table_depth,
        drainage_reference_depth,
        slope_adjustment_micropore_discharge;
        include_matric = true,
        clamp_nonpositive = true,
        hydraulic_gradient = hydraulic_gradient,
        osmotic_factor = osmotic_factor
    ) :
    0.0

macropore_discharge_potential =
    do_macropore_discharge ?
    water_table_discharge_potential(
        0.0,
        osmotic_potential,
        macropore_water_depth,
        water_table_depth,
        drainage_reference_depth,
        slope_adjustment_other;
        include_matric = false,
        clamp_nonpositive = false,
        hydraulic_gradient = hydraulic_gradient,
        osmotic_factor = osmotic_factor
    ) :
    0.0

# -------------------------------------------------------------------------
# 3. Recharge potentials: water table -> soil
# -------------------------------------------------------------------------

micropore_recharge_potential =
    do_micropore_recharge ?
    water_table_recharge_potential(
        matric_potential,
        osmotic_potential,
        layer_mid_depth,
        water_table_depth,
        slope_adjustment_other;
        include_matric = true,
        clamp_nonnegative = true,
        hydraulic_gradient = hydraulic_gradient,
        osmotic_factor = osmotic_factor
    ) :
    0.0

macropore_recharge_potential =
    do_macropore_recharge ?
    water_table_recharge_potential(
        0.0,
        osmotic_potential,
        macropore_water_depth,
        water_table_depth,
        slope_adjustment_other;
        include_matric = false,
        clamp_nonnegative = false,
        hydraulic_gradient = hydraulic_gradient,
        osmotic_factor = osmotic_factor
    ) :
    0.0


# -------------------------------------------------------------------------
# 4. Convert pressure potentials to boundary water fluxes
# -------------------------------------------------------------------------

micropore_discharge_flux =
    boundary_water_flux_from_potential(
        micropore_discharge_potential,
        micropore_conductivity,
        boundary_area,
        fraction_above_water_table,
        micropore_resistance_factor,
        boundary_activity_factor,
        time_factor
    )

macropore_discharge_flux =
    boundary_water_flux_from_potential(
        macropore_discharge_potential,
        macropore_conductivity,
        boundary_area,
        fraction_above_water_table,
        macropore_resistance_factor,
        boundary_activity_factor,
        time_factor
    )

micropore_recharge_flux =
    boundary_water_flux_from_potential(
        micropore_recharge_potential,
        micropore_conductivity,
        boundary_area,
        fraction_below_water_table,
        micropore_resistance_factor,
        boundary_activity_factor,
        time_factor
    )

macropore_recharge_flux =
    boundary_water_flux_from_potential(
        macropore_recharge_potential,
        macropore_conductivity,
        boundary_area,
        fraction_below_water_table,
        macropore_resistance_factor,
        boundary_activity_factor,
        time_factor
    )


# -------------------------------------------------------------------------
# 5. Return all process diagnostics
# -------------------------------------------------------------------------

return (
    slope_adjustment_micropore_discharge =
        slope_adjustment_micropore_discharge,

    slope_adjustment_other =
        slope_adjustment_other,

    micropore_discharge_potential =
        micropore_discharge_potential,

    macropore_discharge_potential =
        macropore_discharge_potential,

    micropore_recharge_potential =
        micropore_recharge_potential,

    macropore_recharge_potential =
        macropore_recharge_potential,

    micropore_discharge_flux =
        micropore_discharge_flux,

    macropore_discharge_flux =
        macropore_discharge_flux,

    micropore_recharge_flux =
        micropore_recharge_flux,

    macropore_recharge_flux =
        macropore_recharge_flux,

    total_micropore_flux =
        micropore_discharge_flux + micropore_recharge_flux,

    total_macropore_flux =
        macropore_discharge_flux + macropore_recharge_flux,

    total_water_flux =
        micropore_discharge_flux +
        macropore_discharge_flux +
        micropore_recharge_flux +
        macropore_recharge_flux,

    # CODEX DEBUG BEGIN: no raw-flux heat outputs; heat must use WF9's limited fluxes.
    # CODEX DEBUG END
)
end
