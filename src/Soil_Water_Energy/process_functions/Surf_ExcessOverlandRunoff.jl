# =============================================================================
# SurfaceExcessOverlandRunoff.jl
#
# Process-level function for partitioning excess surface water/ice and
# calculating overland runoff from free surface liquid water.
#
# Design:
#   - pure function
#   - one surface/litter state at a time
#   - no mutation
#   - no @unpack
#   - no direct access to QRM / QRV / XVOLTM / XVOLWM / XVOLIM arrays
# =============================================================================


@doc raw"""
    surface_excess_overland_runoff_process(...)

Partition excess surface litter water/ice above litter holding capacity and
calculate overland runoff from free surface liquid water.

The candidate post-flux surface liquid water is

```math
W_l' =
W_l + F_{surf} + F_{evap} + F_{ft}
```

where ``W_l`` is current litter liquid water, ``F_{surf}`` is the net liquid
flux into surface litter storage, ``F_{evap}`` is vapor-liquid exchange, and
``F_{ft}`` is freeze-thaw liquid-water-equivalent flux.

The candidate post-flux surface ice is

```math
I_l' =
I_l - \frac{F_{ft}}{\rho_i}
```

where ``\rho_i`` is the ice-density conversion factor.

The total water plus ice above litter holding capacity is partitioned into
free liquid and free ice using the legacy ratio-preserving rule:

```math
E =
\max(0, W_l' + I_l' - C_l)
```

```math
W_f =
\max\left(0, W_l' - \frac{W_l'}{W_l' + I_l'} C_l\right)
```

```math
I_f =
\max\left(0, I_l' - \frac{I_l'}{W_l' + I_l'} C_l\right)
```

where ``C_l`` is litter water holding capacity.

Overland runoff occurs only when free surface liquid exceeds ground-water
retention capacity:

```math
W_x =
W_f - W_r
```

If ``W_x > 0``, the water depth and hydraulic radius are

```math
D = \frac{W_x}{A_s}
```

and

```math
R = \frac{D}{2.828}
```

The surface velocity is

```math
v =
\frac{
    R^{0.67} \sqrt{\max(0, S)}
}{
    z_m
}
```

and the potential runoff is

```math
Q^* =
v D A_f 3.6 \times 10^3 \Delta t
```

The final runoff is limited by excess free liquid and available liquid water:

```math
Q =
\min(Q^*, W_x f_r, W_{avail} f_r)
```

# note
 
```julia
available_liquid_for_runoff = max(0.0, surfW_Micpore_L + watFlux_littFT)
```

rather than the full post-flux liquid storage. This function preserves that
behavior by default.

# Returns

A named tuple containing post-flux surface liquid/ice, total excess water/ice,
free liquid, free ice, runoff, runoff velocity, and runoff-limiting diagnostics.
"""
function surface_excess_overland_runoff_process(;
    # -------------------------------------------------------------------------
    # Current surface litter state
    # -------------------------------------------------------------------------
    surface_liquid_water,
    surface_ice,

    # -------------------------------------------------------------------------
    # Fluxes already diagnosed for this substep/hour
    # -------------------------------------------------------------------------
    net_surface_liquid_flux,
    vapor_liquid_flux,
    freeze_thaw_water_flux,

    # -------------------------------------------------------------------------
    # Surface storage properties
    # -------------------------------------------------------------------------
    litter_water_holding_capacity,
    ground_water_retention_capacity,

    # -------------------------------------------------------------------------
    # Runoff geometry and forcing
    # -------------------------------------------------------------------------
    surface_area_for_depth,
    runoff_flow_area,
    slope,
    roughness_scale,

    # -------------------------------------------------------------------------
    # Time and limiting factors
    # -------------------------------------------------------------------------
    runoff_time_factor,
    storage_removal_fraction,

    # -------------------------------------------------------------------------
    # Constants
    # -------------------------------------------------------------------------
    ice_density_factor = soil_iceDensty,

    # Preserve current WF6 runoff availability expression
    use_legacy_available_liquid = true,

    # Optional override if you later want to pass availability explicitly
    available_liquid_override = nothing
)
    # -------------------------------------------------------------------------
    # 1. Candidate post-flux surface liquid and ice
    #
    # Preserve WF6:
    #
    #     VOLW10 = surfW_Micpore_L + FLWRL + surfLitt_EVP + watFlux_littFT
    #     VOLI10 = surfIce_Micpore_L - watFlux_littFT / soil_iceDensty
    # -------------------------------------------------------------------------

    surface_liquid_after_fluxes =
        surface_liquid_water +
        net_surface_liquid_flux +
        vapor_liquid_flux +
        freeze_thaw_water_flux

    surface_ice_after_phase_change =
        surface_ice -
        freeze_thaw_water_flux / ice_density_factor


    # -------------------------------------------------------------------------
    # 2. Partition excess water/ice above litter holding capacity
    # -------------------------------------------------------------------------

    surface_excess =
        partition_surface_excess_liquid_ice(
            surface_liquid_after_fluxes,
            surface_ice_after_phase_change,
            litter_water_holding_capacity
        )

    excess_surface_water_ice =
        surface_excess.excess_total

    free_surface_liquid =
        surface_excess.free_liquid

    free_surface_ice =
        surface_excess.free_ice


    # -------------------------------------------------------------------------
    # 3. Available liquid for runoff limiting
    #
    # Preserve WF6:
    #
    #     VOLW1X = max(0.0, surfW_Micpore_L + watFlux_littFT)
    #
    # This intentionally does not include FLWRL or surfLitt_EVP.
    # -------------------------------------------------------------------------

    available_liquid_for_runoff =
        if available_liquid_override === nothing
            if use_legacy_available_liquid
                max(
                    0.0,
                    surface_liquid_water + freeze_thaw_water_flux
                )
            else
                max(
                    0.0,
                    surface_liquid_after_fluxes
                )
            end
        else
            max(0.0, available_liquid_override)
        end


    # -------------------------------------------------------------------------
    # 4. Surface overland runoff
    # -------------------------------------------------------------------------

    runoff =
        surface_overland_runoff(
            free_surface_liquid,
            ground_water_retention_capacity,
            surface_area_for_depth,
            runoff_flow_area,
            slope,
            roughness_scale,
            runoff_time_factor,
            storage_removal_fraction,
            available_liquid_for_runoff
        )


    # -------------------------------------------------------------------------
    # 5. Return diagnostics
    # -------------------------------------------------------------------------

    return (
        surface_liquid_after_fluxes =
            surface_liquid_after_fluxes,

        surface_ice_after_phase_change =
            surface_ice_after_phase_change,

        excess_surface_water_ice =
            excess_surface_water_ice,

        free_surface_liquid =
            free_surface_liquid,

        free_surface_ice =
            free_surface_ice,

        available_liquid_for_runoff =
            available_liquid_for_runoff,

        runoff =
            runoff.runoff,

        runoff_velocity =
            runoff.velocity
    )
end