# =============================================================================
# PhaseChange.jl

# Pure phase-change helper equations for water-energy exchange.

# Design rules:

# - no mutation
# - no direct access to waterVar_copy
# - no process-level state update

# - functions return algebraic phase-change quantities only

# Vapor / evaporation:
# positive vapor exchange  -> condensation/deposition to surface
# negative vapor exchange  -> evaporation/sublimation from surface

# Freeze-thaw heat:
# positive latent heat flux -> freezing, heat released to storage
# negative latent heat flux -> thawing/melting, heat consumed from storage

# Freeze-thaw water flux:
# water_flux = -latent_heat_flux / latent_heat_fusion

# positive water_flux -> ice melts to liquid water
# negative water_flux -> liquid water freezes to ice
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Vapor exchange limits
# -----------------------------------------------------------------------------

"""
limit_evaporation_by_liquid(
potential_vapor_exchange,
liquid_available,
removal_fraction
)

Limit vapor exchange by available liquid water.

Physical meaning:

```
potential_vapor_exchange > 0:
    condensation to surface; no liquid-water availability limit needed.

potential_vapor_exchange < 0:
    evaporation from surface; cannot evaporate more than available liquid water.
```

Legacy form:

```
EVAP = max(potential, -max(0, liquid_available * removal_fraction))
```

where `removal_fraction` is usually XNPAX, XNPXX, etc.
"""
function limit_evaporation_by_liquid(
    potential_vapor_exchange,
    liquid_available,
    removal_fraction
    )
    return max(
    potential_vapor_exchange,
    -max(0.0, liquid_available * removal_fraction)
    )
end

"""
partition_vapor_exchange_liquid_solid(
potential_vapor_exchange,
liquid_available,
solid_available,
removal_fraction
)

Partition vapor exchange between liquid water and solid snow/ice.

This is useful for snowpack vapor exchange, where evaporation/sublimation is
first limited by available liquid water and then by available snow/ice.

Legacy structure:

```
EVAPW = max(EVAPT, -max(0, liquid_available * removal_fraction))
residual = min(0, EVAPT - EVAPW)
EVAPS = max(residual, -max(0, solid_available * removal_fraction))
```

Physical meaning:
positive potential exchange:
condensation/deposition is assigned to liquid_vapor_exchange.

```
negative potential exchange:
    evaporation uses liquid water first;
    remaining vapor demand is supplied by sublimation from snow/ice.
```

Returns a NamedTuple:
liquid_vapor_exchange
solid_vapor_exchange
residual_unmet_exchange
"""
function partition_vapor_exchange_liquid_solid(
    potential_vapor_exchange,
    liquid_available,
    solid_available,
    removal_fraction
    )
    liquid_vapor_exchange =
    max(
    potential_vapor_exchange,
    -max(0.0, liquid_available * removal_fraction)
    )

    # CODEX DEBUG BEGIN: remove stray Markdown fence from executable function body.
    residual =
        min(0.0, potential_vapor_exchange - liquid_vapor_exchange)

    solid_vapor_exchange =
        max(
            residual,
            -max(0.0, solid_available * removal_fraction)
        )

    residual_unmet_exchange =
        potential_vapor_exchange -
        liquid_vapor_exchange -
        solid_vapor_exchange

    return (
        liquid_vapor_exchange = liquid_vapor_exchange,
        solid_vapor_exchange = solid_vapor_exchange,
        residual_unmet_exchange = residual_unmet_exchange
    )
    # CODEX DEBUG END

end

"""
pore_vapor_condensation_potential(
vapor_storage,
equilibrium_vapor_concentration,
air_filled_volume
)

Potential vapor exchange inside an air-filled pore space.

Physical form:

```
E_potential = V_vapor_current - C_eq * V_air
```

where

```
V_vapor_current = current vapor storage
C_eq            = equilibrium vapor concentration
V_air           = air-filled pore volume
```

If positive, vapor storage exceeds equilibrium and condensation occurs.
If negative, vapor storage is below equilibrium and evaporation/sublimation occurs.
"""
function pore_vapor_condensation_potential(
    vapor_storage,
    equilibrium_vapor_concentration,
    air_filled_volume
    )
    return vapor_storage - equilibrium_vapor_concentration * air_filled_volume
end

# -----------------------------------------------------------------------------

# 3. Freezing-point depression

# -----------------------------------------------------------------------------

"""
freezing_temperature_from_water_potential(
water_potential;
latent_heat_fusion = 333.0
)

Freezing temperature modified by water potential.

Legacy form used for litter and soil:

```
TFREEZ = -9.0959E+04 / (ψ - 333.0)
```

where ψ is water potential and 333.0 is the latent heat of fusion in the
legacy model units.

Physical meaning:
more negative water potential depresses the freezing temperature below
273.15 K, so water in dry/salty/tightly bound pores freezes at a lower
temperature than free water.
"""
function freezing_temperature_from_water_potential(
    water_potential;
    latent_heat_fusion = 333.0
    )
    return -9.0959e4 / (water_potential - latent_heat_fusion)
end

# -----------------------------------------------------------------------------

# 4. Freeze-thaw activation

# -----------------------------------------------------------------------------

"""
freeze_thaw_is_active(
temperature,
freezing_temperature,
liquid_water,
ice_amount,
reference_volume;
tiny = tiny_num
)

Determine whether freeze-thaw should occur.

Legacy logic:

```
freezing can occur when:
    T < T_freeze and liquid water is available

thawing can occur when:
    T > T_freeze and ice is available
```

The reference volume is used only to define a small numerical threshold.
"""
function freeze_thaw_is_active(
    temperature,
    freezing_temperature,
    liquid_water,
    ice_amount,
    reference_volume;
    tiny = tiny_num
    )
    liquid_threshold = tiny * reference_volume
    ice_threshold    = tiny * reference_volume

    # CODEX DEBUG BEGIN: remove stray Markdown fence from executable function body.
    return (
        (temperature < freezing_temperature && liquid_water > liquid_threshold) ||
        (temperature > freezing_temperature && ice_amount > ice_threshold)
    )
    # CODEX DEBUG END

end

# -----------------------------------------------------------------------------

# 5. Freeze-thaw latent heat demand

# -----------------------------------------------------------------------------

"""
freeze_thaw_heat_potential(
heat_capacity,
temperature,
freezing_temperature,
water_potential,
time_factor;
thermal_slope_coeff = 6.2913e-3,
water_potential_coeff = 0.10
)

Potential latent heat available for freeze-thaw.

Legacy soil/litter form:

```
H_potential​=C * (Tfreeze​−T) /(1+6.2913*10−3*Tfreeze​)(1−aψ)* ​Δt
```

where

```
C          = heat capacity of the water/ice-containing domain
T_freeze   = effective freezing temperature
T          = current temperature
ψ          = water potential
a          = water-potential coefficient, either 0.10 or 0.00 in your code
Δt         = substep time factor
```

Sign convention:
positive H_potential:
    temperature below freezing point → freezing tendency

negative H_potential:
    temperature above freezing point → thawing/melting tendency

Notes:
- Surface litter and soil use water_potential_coeff = 0.10 in parts of
the legacy code.
- Some micropore soil code uses water_potential_coeff = 0.00.
Pass water_potential_coeff = 0.00 when preserving that behavior.
"""
function freeze_thaw_heat_potential(
    heat_capacity,
    temperature,
    freezing_temperature,
    water_potential,
    time_factor;
    thermal_slope_coeff = 6.2913e-3,
    water_potential_coeff = 0.10
    )
    denominator =
    (1.0 + thermal_slope_coeff * freezing_temperature) *
    (1.0 - water_potential_coeff * water_potential)

    # CODEX DEBUG BEGIN: remove stray Markdown fence from executable function body.
    return heat_capacity *
        (freezing_temperature - temperature) /
        denominator *
        time_factor
    # CODEX DEBUG END

end

"""
snow_freeze_thaw_heat_potential(
heat_capacity,
temperature,
time_factor;
freezing_temperature = 273.15,
snow_damping_factor = 2.7185
)

Potential latent heat available for snowpack freeze-thaw.

```
H_potential = C * (273.15 - T) / 2.7185 * time_factor
```

Sign convention:
H_potential > 0:
    snow temperature below 273.15 K
    liquid water tends to freeze

H_potential < 0:
    snow temperature above 273.15 K
    snow/ice tends to melt
"""
function snow_freeze_thaw_heat_potential(
    heat_capacity,
    temperature,
    time_factor;
    freezing_temperature = 273.15,
    snow_damping_factor = 2.7185
    )
    return heat_capacity *
    (freezing_temperature - temperature) /
    snow_damping_factor *
    time_factor
end

# -----------------------------------------------------------------------------

# 6. Freeze-thaw limiting by available water or ice

# -----------------------------------------------------------------------------

"""
limit_freeze_thaw_heat(
potential_latent_heat,
liquid_water,
ice_amount,
time_factor;
latent_heat_fusion = 333.0,
ice_density_factor = 1.0
)

Limit freeze-thaw latent heat by available liquid water or ice.

Physical meaning:

```
potential_latent_heat > 0:
    freezing tendency; cannot freeze more liquid water than available.

potential_latent_heat < 0:
    thawing/melting tendency; cannot melt more ice than available.
```

Legacy forms:

```
freezing:
    H_limited = min( L_fusion * liquid_water * time_factor,
                     H_potential )

thawing:
    H_limited = max(-L_fusion * ice_density_factor * ice_amount * time_factor,
                     H_potential )
```

Sign convention:
positive H_limited -> freezing, canot freeze more liquid water than available
negative H_limited -> thawing/melting, cannot melt more ice than available 
"""
function limit_freeze_thaw_heat(
    potential_latent_heat,
    liquid_water,
    ice_amount,
    time_factor;
    latent_heat_fusion = 333.0,
    ice_density_factor = 1.0
    )
    if potential_latent_heat < 0.0
        return max(
        -latent_heat_fusion * ice_density_factor * ice_amount * time_factor,
        potential_latent_heat
        )
    else
        return min(
        latent_heat_fusion * liquid_water * time_factor,
        potential_latent_heat
        )
    end
end

"""
freeze_thaw_water_flux(
latent_heat_flux;
latent_heat_fusion = 333.0
)

Convert freeze-thaw latent heat flux to liquid-water-equivalent flux.

Legacy form:

```
water_flux = -H_freeze_thaw / L_fusion
```

Sign convention:
positive water_flux -> melting/thawing adds liquid water
negative water_flux -> freezing removes liquid water
"""
freeze_thaw_water_flux(
    latent_heat_flux;
    latent_heat_fusion = 333.0
    ) = -latent_heat_flux / latent_heat_fusion

"""
freeze_thaw_limited_fluxes(
potential_latent_heat,
liquid_water,
ice_amount,
time_factor;
latent_heat_fusion = 333.0,
ice_density_factor = 1.0
)

Return both limited latent heat and liquid-water-equivalent freeze-thaw flux.

potential freeze-thaw heat
→ limited by available water or ice
→ converted to liquid-water-equivalent flux

This is suitable for surface litter and soil freeze-thaw updates.
"""
function freeze_thaw_limited_fluxes(
    potential_latent_heat,
    liquid_water,
    ice_amount,
    time_factor;
    latent_heat_fusion = 333.0,
    ice_density_factor = 1.0
    )
    latent_heat_flux = limit_freeze_thaw_heat(
    potential_latent_heat,
    liquid_water,
    ice_amount,
    time_factor;
    latent_heat_fusion = latent_heat_fusion,
    ice_density_factor = ice_density_factor
    )

    # CODEX DEBUG BEGIN: remove stray Markdown fence from executable function body.
    water_flux = freeze_thaw_water_flux(
        latent_heat_flux;
        latent_heat_fusion = latent_heat_fusion
    )

    return (
        latent_heat_flux = latent_heat_flux,
        water_flux = water_flux
    )
    # CODEX DEBUG END

end

# -----------------------------------------------------------------------------

# 7. Snow-specific partition of thawing between snow and ice pools

# -----------------------------------------------------------------------------

"""
snow_solid_water_fractions(
snow_volume,
ice_volume;
ice_density_factor = 1.0,
tiny = tiny_num
)

Calculate the fractional contribution of snow and ice to total solid water.

Legacy snow thawing uses:

```
total_solid_water = snow_volume + ice_volume * ice_density_factor

snow_fraction = snow_volume / total_solid_water
ice_fraction  = ice_volume * ice_density_factor / total_solid_water
```

These fractions are then used to partition thawing between snow and ice pools.
"""
function snow_solid_water_fractions(
    snow_volume,
    ice_volume;
    ice_density_factor = 1.0,
    tiny = tiny_num
    )
    total_solid_water =
    snow_volume + ice_volume * ice_density_factor

    # CODEX DEBUG BEGIN: remove stray Markdown fence from executable function body.
    if total_solid_water > tiny
        return (
            snow_fraction = snow_volume / total_solid_water,
            ice_fraction = ice_volume * ice_density_factor / total_solid_water,
            total_solid_water = total_solid_water
        )
    else
        return (
            snow_fraction = 0.0,
            ice_fraction = 0.0,
            total_solid_water = total_solid_water
        )
    end
    # CODEX DEBUG END

end

"""
    snow_freeze_thaw_limited_fluxes(
        heat_capacity,
        temperature,
        liquid_water,
        snow_volume,
        ice_volume,
        time_factor;
        freezing_temperature = 273.15,
        snow_damping_factor = 2.7185,
        latent_heat_fusion = 333.0,
        ice_density_factor = 1.0,
        tiny = tiny_num
    )

Calculate snowpack freeze-thaw latent heat and separate water-equivalent
fluxes for snow and ice pools.

This preserves the legacy snow freeze-thaw logic:

    H_potential = C * (T_freeze - T) / damping * Δt

If H_potential < 0:
    thawing/melting occurs.
    Melt demand is limited by total solid water:
        snow_volume + ice_volume * ice_density_factor
    The melt flux is partitioned between snow and ice pools.

If H_potential > 0:
    freezing occurs.
    Freezing is limited by available liquid water.
    Frozen water is added to the ice pool.

Returns:
    latent_heat_flux
    snow_water_flux
    ice_water_flux
"""
function snow_freeze_thaw_limited_fluxes(
    heat_capacity,
    temperature,
    liquid_water,
    snow_volume,
    ice_volume,
    time_factor;
    freezing_temperature = 273.15,
    snow_damping_factor = 2.7185,
    latent_heat_fusion = 333.0,
    ice_density_factor = 1.0,
    tiny = tiny_num
)
    potential_latent_heat =
        heat_capacity *
        (freezing_temperature - temperature) /
        snow_damping_factor *
        time_factor

    if potential_latent_heat < 0.0
        # CODEX DEBUG BEGIN: calculate snow/ice fractions with shared solid-water helper.
        solid_parts = snow_solid_water_fractions(
            snow_volume,
            ice_volume;
            ice_density_factor = ice_density_factor,
            tiny = tiny,
        )
        total_solid_water = solid_parts.total_solid_water
        snow_fraction = solid_parts.snow_fraction
        ice_fraction = solid_parts.ice_fraction
        # CODEX DEBUG END

        latent_heat_flux = max(
            -latent_heat_fusion * total_solid_water * time_factor,
            potential_latent_heat
        )

        snow_water_flux =
            -latent_heat_flux * snow_fraction / latent_heat_fusion

        ice_water_flux =
            -latent_heat_flux * ice_fraction / latent_heat_fusion

    else
        latent_heat_flux = min(
            latent_heat_fusion * liquid_water * time_factor,
            potential_latent_heat
        )

        snow_water_flux = 0.0
        ice_water_flux  = -latent_heat_flux / latent_heat_fusion
    end

    return (
        latent_heat_flux = latent_heat_flux,
        snow_water_flux = snow_water_flux,
        ice_water_flux = ice_water_flux,
        potential_latent_heat = potential_latent_heat
    )
end
