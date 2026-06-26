# =============================================================================
# SoilLayerPhaseChange.jl
#
# Process-level function for local soil-layer phase change:
#   1. vapor <-> liquid exchange in soil micropores
#   2. liquid <-> ice exchange in soil micropores
#   3. liquid <-> ice exchange in soil macropores
#
# Design:
#   - pure function
#   - one soil layer at a time
#   - no mutation
#   - no @unpack
#   - no direct access to waterVar_copy
# =============================================================================


@doc raw"""
    soil_layer_phase_change_process(...)

Calculate local phase-change fluxes in one soil layer.

This process includes:

1. Vapor--liquid exchange in the soil micropore air space.
2. Liquid--ice phase change in soil micropores.
3. Liquid--ice phase change in soil macropores.

The soil water potential used for freezing-point depression is

```math
\psi_f = \psi_m + \psi_o
```

where `\psi_m` is matric potential and `\psi_o` is osmotic potential.

The freezing temperature is

```math
T_f = T_f(\psi_f)
```

implemented by

```julia
freezing_temperature_from_water_potential(freezing_water_potential)
```

The equilibrium vapor concentration in micropore air is

```math
C_v^{eq} = C_v(T, \psi_v)
```

implemented by

```julia
VP_at_psi(temperature, vapor_exchange_potential)
```

The vapor--liquid flux is limited by available liquid water when evaporation
is negative:

```math
F_v =
\max \left(
    F_v^*,
    - W_l f_r
\right)
```

where `F_v^*` is the potential vapor--liquid flux, `W_l` is available
liquid water, and `f_r` is the removal fraction for the current substep.

Micropore freeze--thaw is driven by a heat potential

```math
H_{ft,\mu}^* =
C_\mu (T_f - T) f_r
```

implemented by

```julia
freeze_thaw_heat_potential(...)
```

and limited by available liquid water or ice through

```math
F_{ft,\mu} =
\frac{H_{ft,\mu}}{L_f}
```

implemented by

```julia
freeze_thaw_limited_fluxes(...)
```

Macropore freeze--thaw follows the same structure, but uses macropore heat
capacity

```math
C_M = c_w W_M + c_i I_M
```

where `W_M` and `I_M` are macropore water and ice storage.

# Sign convention

For returned water fluxes:

* positive freeze--thaw water flux means freezing, liquid water becomes ice
* negative freeze--thaw water flux means thawing, ice becomes liquid
* positive vapor-liquid flux follows the sign convention of
  `pore_vapor_condensation_potential(...)`

# Returns

A named tuple containing vapor-liquid fluxes, freeze-thaw fluxes, latent heat
fluxes, freezing diagnostics, and sanitized storages used by the calculation.
"""
function soil_layer_phase_change_process(;
temperature,
vapor_exchange_potential,
matric_potential,
osmotic_potential,

vapor_storage,
liquid_water_micropore,
ice_micropore,
air_micropore,

liquid_water_macropore,
ice_macropore,

layer_volume,
layer_heat_capacity,

removal_fraction,

cpw,
cpi,

latent_heat_vaporization = H_vap,
latent_heat_fusion = 333.0,
ice_density_factor = soil_iceDensty,
reference_ice_temperature = TFice,

# CODEX DEBUG BEGIN: keep thermal and water-potential denominator coefficients independent.
thermal_slope_coeff = 6.2913e-3,
micropore_water_potential_coeff = 0.0,
macropore_water_potential_coeff = 0.10
# CODEX DEBUG END
)
# -------------------------------------------------------------------------
# 1. Sanitize storages used by phase-change limiters
# -------------------------------------------------------------------------


liquid_micropore =
    max(0.0, liquid_water_micropore)

vapor_micropore =
    max(0.0, vapor_storage)

ice_micro =
    max(0.0, ice_micropore)

liquid_macropore =
    max(0.0, liquid_water_macropore)

ice_macro =
    max(0.0, ice_macropore)


# -------------------------------------------------------------------------
# 2. Vapor <-> liquid exchange in soil micropores
# -------------------------------------------------------------------------

if air_micropore > 0.0
    equilibrium_vapor_concentration =
        VP_at_psi(
            temperature,
            vapor_exchange_potential
        )

    potential_vapor_liquid_flux =
        pore_vapor_condensation_potential(
            vapor_micropore,
            equilibrium_vapor_concentration,
            air_micropore
        )

    vapor_liquid_flux =
        max(
            potential_vapor_liquid_flux,
            -liquid_micropore * removal_fraction
        )

    vapor_liquid_latent_heat =
        latent_heat_evaporation(
            vapor_liquid_flux,
            latent_heat_vaporization
        )
else
    equilibrium_vapor_concentration =
        0.0

    potential_vapor_liquid_flux =
        0.0

    vapor_liquid_flux =
        0.0

    vapor_liquid_latent_heat =
        0.0
end


# -------------------------------------------------------------------------
# 3. Freezing-point depression from soil water potential
# -------------------------------------------------------------------------

freezing_water_potential =
    vapor_water_potential(
        matric_potential,
        osmotic_potential
    )

freezing_temperature =
    freezing_temperature_from_water_potential(
        freezing_water_potential
    )


# -------------------------------------------------------------------------
# 4. Micropore liquid <-> ice phase change
# -------------------------------------------------------------------------

if freeze_thaw_is_active(
    temperature,
    freezing_temperature,
    liquid_micropore,
    ice_micro,
    layer_volume
)
    micropore_heat_potential =
        freeze_thaw_heat_potential(
            layer_heat_capacity,
            temperature,
            freezing_temperature,
            freezing_water_potential,
            removal_fraction;
            thermal_slope_coeff = thermal_slope_coeff,
            water_potential_coeff = micropore_water_potential_coeff
        )

    micropore_phase_change =
        freeze_thaw_limited_fluxes(
            micropore_heat_potential,
            liquid_micropore,
            ice_micro,
            removal_fraction;
            latent_heat_fusion = latent_heat_fusion,
            ice_density_factor = ice_density_factor
        )

    micropore_freeze_thaw_heat =
        micropore_phase_change.latent_heat_flux

    micropore_freeze_thaw_water =
        micropore_phase_change.water_flux
else
    micropore_heat_potential =
        0.0

    micropore_freeze_thaw_heat =
        0.0

    micropore_freeze_thaw_water =
        0.0
end


# -------------------------------------------------------------------------
# 5. Macropore liquid <-> ice phase change
# -------------------------------------------------------------------------

macropore_heat_capacity =
    cpw * liquid_macropore +
    cpi * ice_macro

if freeze_thaw_is_active(
    temperature,
    reference_ice_temperature,
    liquid_macropore,
    ice_macro,
    layer_volume
)
    macropore_heat_potential =
        freeze_thaw_heat_potential(
            macropore_heat_capacity,
            temperature,
            freezing_temperature,
            freezing_water_potential,
            removal_fraction;
            thermal_slope_coeff = thermal_slope_coeff,
            water_potential_coeff = macropore_water_potential_coeff
        )

    macropore_phase_change =
        freeze_thaw_limited_fluxes(
            macropore_heat_potential,
            liquid_macropore,
            ice_macro,
            removal_fraction;
            latent_heat_fusion = latent_heat_fusion,
            ice_density_factor = ice_density_factor
        )

    macropore_freeze_thaw_heat =
        macropore_phase_change.latent_heat_flux

    macropore_freeze_thaw_water =
        macropore_phase_change.water_flux
else
    macropore_heat_potential =
        0.0

    macropore_freeze_thaw_heat =
        0.0

    macropore_freeze_thaw_water =
        0.0
end


# -------------------------------------------------------------------------
# 6. Return process fluxes and diagnostics
# -------------------------------------------------------------------------

return (
    equilibrium_vapor_concentration =
        equilibrium_vapor_concentration,

    potential_vapor_liquid_flux =
        potential_vapor_liquid_flux,

    vapor_liquid_flux =
        vapor_liquid_flux,

    vapor_liquid_latent_heat =
        vapor_liquid_latent_heat,

    freezing_water_potential =
        freezing_water_potential,

    freezing_temperature =
        freezing_temperature,

    micropore_heat_potential =
        micropore_heat_potential,

    micropore_freeze_thaw_water =
        micropore_freeze_thaw_water,

    micropore_freeze_thaw_heat =
        micropore_freeze_thaw_heat,

    macropore_heat_capacity =
        macropore_heat_capacity,

    macropore_heat_potential =
        macropore_heat_potential,

    macropore_freeze_thaw_water =
        macropore_freeze_thaw_water,

    macropore_freeze_thaw_heat =
        macropore_freeze_thaw_heat,

    total_freeze_thaw_heat =
        micropore_freeze_thaw_heat +
        macropore_freeze_thaw_heat,

    liquid_micropore =
        liquid_micropore,

    vapor_micropore =
        vapor_micropore,

    ice_micropore =
        ice_micro,

    liquid_macropore =
        liquid_macropore,

    ice_macropore =
        ice_macro
)
end

