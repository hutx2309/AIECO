# =============================================================================
# StorageUpdates.jl
#

# Pure storage/state-update algebra for the water-energy module.
#
# Scope:
# - convert net fluxes into updated snow, litter, and soil storages
# - calculate air-filled volume after water/ice updates
# - calculate diagnostic water/ice/air fractions
# - update temperature from old energy + net heat
#

# Design rules:
# - no mutation
# - no @unpack
# - no direct access to waterVar_copy
# - no direct model-state assignment
# - functions return updated values as NamedTuples
#
# Sign conventions:
# - positive water/vapor/ice flux into a storage increases that storage
# - positive freeze-thaw water flux means melting/thawing adds liquid water
# - negative freeze-thaw water flux means freezing removes liquid water
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Generic storage and fraction helpers
# -----------------------------------------------------------------------------

@doc raw"""
air_storage_from_pore_volume(pore_volume, water_storage, ice_storage)

Calculate air-filled storage from pore volume after water and ice occupy space.

Legacy pattern:

```math
A_{air} = \max(0, V_p - W - I)
```

Inputs:
pore_volume    : available pore volume
water_storage  : liquid water storage
ice_storage    : ice storage

Output:
nonnegative air-filled storage
"""
function air_storage_from_pore_volume(
    pore_volume,
    water_storage,
    ice_storage
    )
    return max(0.0, pore_volume - water_storage - ice_storage)
end

@doc raw"""
storage_fraction(storage, reference_volume; fallback = 0.0, tiny = tiny_num2)

Convert a storage amount to a nonnegative fraction using a reference volume.

Legacy pattern:

```math
f =
\begin{cases}
\max\left(0, \frac{S}{V_{ref}}\right), & V_{ref} > \epsilon \\
f_{fallback}, & V_{ref} \le \epsilon
\end{cases}
```

This helper should be used for diagnostic water/ice/air fractions, not to
silently alter conserved storage amounts.
"""
function storage_fraction(
    storage,
    reference_volume;
    fallback = 0.0,
    tiny = tiny_num2
    )
    if reference_volume > tiny
    return max(0.0, storage / reference_volume)
    else
    return fallback
    end
end

@doc raw"""
storage_fractions_water_ice_air(
water_storage,
ice_storage,
air_storage,
reference_volume;
fallback_water = 0.0,
fallback_ice = 0.0,
fallback_air = 1.0,
air_modifier = 1.0,
tiny = tiny_num2
)

Calculate diagnostic water, ice, and air fractions from storages.

The optional `air_modifier` is useful for surface litter, where the apparent
air-filled fraction is reduced by excess surface water/ice cover.

Returns:
water_fraction
ice_fraction
air_fraction
"""
function storage_fractions_water_ice_air(
    water_storage,
    ice_storage,
    air_storage,
    reference_volume;
    fallback_water = 0.0,
    fallback_ice = 0.0,
    fallback_air = 1.0,
    air_modifier = 1.0,
    tiny = tiny_num2
    )
    if reference_volume > tiny
    return (
    water_fraction = max(0.0, water_storage / reference_volume),
    ice_fraction   = max(0.0, ice_storage / reference_volume),
    air_fraction   = max(0.0, air_storage / reference_volume) *
    max(0.0, air_modifier)
    )
    else
    return (
    water_fraction = fallback_water,
    ice_fraction   = fallback_ice,
    air_fraction   = fallback_air
    )
    end
end
# CLAUDE DEBUG BEGIN: define missing helper. `litter_storage_fractions` (below) calls
# `surface_air_modifier_from_excess` but it was never defined → UndefVarError at runtime.
# Formula per the litter_storage_fractions docstring: air_fraction *= max(0, 1 - excess/capacity).
# Guard capacity > tiny to avoid div-by-zero; degenerate capacity returns 1.0 (no modifier).
@doc raw"""
surface_air_modifier_from_excess(
    excess_surface_water_ice,
    surface_water_holding_capacity;
    tiny = tiny_num2
)

Legacy surface-litter air-fraction reduction factor:

```math
m_{air} = \max\left(0, 1 - \frac{W_{excess}}{W_{hold}}\right)
```

Returns 1.0 (no reduction) when `holding_capacity <= tiny`.
"""
function surface_air_modifier_from_excess(
    excess_surface_water_ice,
    surface_water_holding_capacity;
    tiny = tiny_num2
    )
    if surface_water_holding_capacity > tiny
    return max(0.0, 1.0 - excess_surface_water_ice / surface_water_holding_capacity)
    else
    return 1.0
    end
end
# CLAUDE DEBUG END

@doc raw"""
litter_storage_fractions(
water,
ice,
air,
dry_litter_volume,
excess_surface_water_ice,
surface_water_holding_capacity;
tiny = tiny_num2
)

Calculate surface litter water, ice, and air volume fractions.

The air fraction includes the legacy excess-surface-water reduction factor:

```math
f_{air} \leftarrow f_{air}\max\left(0, 1 - \frac{W_{excess}}{W_{hold}}\right)
```

"""
function litter_storage_fractions(
    water,
    ice,
    air,
    dry_litter_volume,
    excess_surface_water_ice,
    surface_water_holding_capacity;
    tiny = tiny_num2
)
    air_modifier =
        surface_air_modifier_from_excess(
            excess_surface_water_ice,
            surface_water_holding_capacity;
            tiny = tiny
        )

    return storage_fractions_water_ice_air(
        water,
        ice,
        air,
        dry_litter_volume;
        fallback_water = 0.0,
        fallback_ice = 0.0,
        fallback_air = 1.0,
        air_modifier = air_modifier,
        tiny = tiny
    )
end
@doc raw"""
soil_storage_fractions(
micropore_water,
micropore_ice,
macropore_water,
macropore_ice,
air_micropore,
air_macropore,
soil_reference_volume,
macropore_volume;
tiny = tiny_num2
)

Calculate diagnostic total soil water, ice, and air fractions.

Legacy denominator:

```math
V_{total} = V_{soil} + V_M
```

"""
function soil_storage_fractions(
    micropore_water,
    micropore_ice,
    macropore_water,
    macropore_ice,
    air_micropore,
    air_macropore,
    soil_reference_volume,
    macropore_volume;
    fallback_water = 0.0,
    fallback_ice = 0.0,
    fallback_air = 0.0,
    tiny = tiny_num2
)
    total_volume =
        soil_reference_volume + macropore_volume

    return storage_fractions_water_ice_air(
        micropore_water + macropore_water,
        micropore_ice + macropore_ice,
        air_micropore + air_macropore,
        total_volume;
        fallback_water = fallback_water,
        fallback_ice = fallback_ice,
        fallback_air = fallback_air,
        tiny = tiny
    )
end

# -----------------------------------------------------------------------------
# 3. Snow storage update
# -----------------------------------------------------------------------------

@doc raw"""
update_snow_layer_storage(
snow_volume,
liquid_water,
ice_volume,
vapor_volume,
net_snow_flux,
net_water_flux,
net_ice_flux,
net_vapor_flux,
freeze_thaw_snow_water_flux,
freeze_thaw_ice_water_flux,
vapor_flux_from_snow,
vapor_flux_from_liquid;
ice_density_factor
)

Update snow-layer snow, liquid water, ice, and vapor storages.

This preserves the WF10 snow update algebra:

```math
\begin{aligned}
S' &= S + F_S - F_{ft,S} + F_{v,S} \\
W' &= W + F_W + F_{ft,S} + F_{ft,I} + F_{v,W} \\
I' &= I + F_I - \frac{F_{ft,I}}{\rho_i} \\
V' &= V + F_V - F_{v,S} - F_{v,W}
\end{aligned}
```

Returns:
snow_volume
liquid_water
ice_volume
vapor_volume
"""
function update_snow_layer_storage(
    snow_volume,
    liquid_water,
    ice_volume,
    vapor_volume,
    net_snow_flux,
    net_water_flux,
    net_ice_flux,
    net_vapor_flux,
    freeze_thaw_snow_water_flux,
    freeze_thaw_ice_water_flux,
    vapor_flux_from_snow,
    vapor_flux_from_liquid;
    ice_density_factor
    )
    new_snow_volume =
    snow_volume +
    net_snow_flux -
    freeze_thaw_snow_water_flux +
    vapor_flux_from_snow

    new_liquid_water =
        liquid_water +
        net_water_flux +
        freeze_thaw_snow_water_flux +
        freeze_thaw_ice_water_flux +
        vapor_flux_from_liquid

    new_ice_volume =
        ice_volume +
        net_ice_flux -
        freeze_thaw_ice_water_flux / ice_density_factor

    new_vapor_volume =
        vapor_volume +
        net_vapor_flux -
        vapor_flux_from_snow -
        vapor_flux_from_liquid

    return (
        snow_volume = new_snow_volume,
        liquid_water = new_liquid_water,
        ice_volume = new_ice_volume,
        vapor_volume = new_vapor_volume
    )
end

@doc raw"""
update_snow_layer_temperature(
old_heat_capacity,
old_temperature,
snow_volume,
liquid_water,
ice_volume,
vapor_volume,
net_heat,
fallback_temperature,
min_heat_capacity;
cps,
cpw,
cpi
)

Update snow-layer heat capacity and temperature after snow storage changes.

Returns:
heat_capacity
temperature
"""
function update_snow_layer_temperature(
    old_heat_capacity,
    old_temperature,
    snow_volume,
    liquid_water,
    ice_volume,
    vapor_volume,
    net_heat,
    fallback_temperature,
    min_heat_capacity;
    cps,
    cpw,
    cpi
    )

    new_heat_capacity =
        heat_capacity_snow(
            snow_volume,
            liquid_water,
            ice_volume,
            vapor_volume;
            cps = cps,
            cpw = cpw,
            cpi = cpi
        )

    new_temperature =
        temperature_from_energy(
            old_heat_capacity,
            old_temperature,
            net_heat,
            new_heat_capacity,
            fallback_temperature,
            min_heat_capacity
        )

    return (
        heat_capacity = new_heat_capacity,
        temperature = new_temperature
    )

end

# -----------------------------------------------------------------------------
# 5. Surface litter micropore storage update
# -----------------------------------------------------------------------------

@doc raw"""
update_litter_micropore_storage(
water,
vapor,
ice,
pore_volume,
liquid_flux,
vapor_flux,
evaporation_flux,
freeze_thaw_water_flux,
runoff_flux;
ice_density_factor
)

Update surface litter micropore water, vapor, ice, and air storages.

This preserves the WF10 litter storage algebra:

```math
\begin{aligned}
W' &= W + F_l + E + F_{ft} + F_r \\
V' &= V + F_v - E \\
I' &= I - \frac{F_{ft}}{\rho_i} \\
A' &= \max(0, V_p - W' - I')
\end{aligned}
```

Notes:
`evaporation_flux` follows the model sign convention:
positive -> condensation to liquid water
negative -> evaporation from liquid water to vapor
"""
function update_litter_micropore_storage(
    water,
    vapor,
    ice,
    pore_volume,
    liquid_flux,
    vapor_flux,
    evaporation_flux,
    freeze_thaw_water_flux,
    runoff_flux;
    ice_density_factor
    )
    new_water =
    water +
    liquid_flux +
    evaporation_flux +
    freeze_thaw_water_flux +
    runoff_flux


    new_vapor =
        vapor +
        vapor_flux -
        evaporation_flux

    new_ice =
        ice -
        freeze_thaw_water_flux / ice_density_factor

    new_air =
        air_storage_from_pore_volume(
            pore_volume,
            new_water,
            new_ice
        )

    return (
        water = new_water,
        vapor = new_vapor,
        ice = new_ice,
        air = new_air
    )

end

@doc raw"""
update_litter_temperature(
old_heat_capacity,
old_temperature,
organic_mass,
water,
vapor,
ice,
net_heat,
fallback_temperature,
min_heat_capacity;
cpo,
cpw,
cpi
)

Update surface litter heat capacity and temperature.

`organic_mass` is usually surf_SOC + surf_charcoal.
"""
function update_litter_temperature(
    old_heat_capacity,
    old_temperature,
    organic_mass,
    water,
    vapor,
    ice,
    net_heat,
    fallback_temperature,
    min_heat_capacity;
    cpo,
    cpw,
    cpi
    )
    new_heat_capacity =
    heat_capacity_litter(
    organic_mass,
    water,
    vapor,
    ice;
    cpo = cpo,
    cpw = cpw,
    cpi = cpi
    )

    new_temperature =
        temperature_from_energy(
            old_heat_capacity,
            old_temperature,
            net_heat,
            new_heat_capacity,
            fallback_temperature,
            min_heat_capacity
        )

    return (
        heat_capacity = new_heat_capacity,
        temperature = new_temperature
    )
end
# -----------------------------------------------------------------------------
# 6. Soil micropore and macropore storage update
# -----------------------------------------------------------------------------

@doc raw"""
update_soil_micropore_storage(
water,
vapor,
ice,
wet_front_water,
liquid_flux,
vapor_flux,
evaporation_flux,
freeze_thaw_water_flux,
infiltration_flux,
subsurface_water_flux;
ice_density_factor
)

Update soil micropore water, vapor, ice, and wet-front water.

This preserves the WF10 micropore storage algebra:

```math
\begin{aligned}
W' &= W + F_l + E + F_{ft} + F_{inf} + F_{sub} \\
V' &= V + F_v - E \\
W_{front}' &= W_{front} + F_{front} + F_{inf} + F_{ft} + F_{sub} \\
I' &= I - \frac{F_{ft}}{\rho_i}
\end{aligned}
```

The wet-front flux itself is passed separately as `wet_front_liquid_flux`
in `update_soil_wet_front(...)` or can be included in the caller.
"""
function update_soil_micropore_storage(
    water,
    vapor,
    ice,
    wet_front_water,
    liquid_flux,
    vapor_flux,
    evaporation_flux,
    freeze_thaw_water_flux,
    infiltration_flux,
    subsurface_water_flux;
    wet_front_liquid_flux = 0.0,
    ice_density_factor
    )
    new_water =
    water +
    liquid_flux +
    evaporation_flux +
    freeze_thaw_water_flux +
    infiltration_flux +
    subsurface_water_flux


    new_vapor =
        vapor +
        vapor_flux -
        evaporation_flux

    new_wet_front_water =
        wet_front_water +
        wet_front_liquid_flux +
        infiltration_flux +
        freeze_thaw_water_flux +
        subsurface_water_flux

    new_wet_front_water =
        min(new_water, new_wet_front_water)

    new_ice =
        ice -
        freeze_thaw_water_flux / ice_density_factor

    return (
        water = new_water,
        vapor = new_vapor,
        ice = new_ice,
        wet_front_water = new_wet_front_water
    )


end

@doc raw"""
update_soil_macropore_storage(
water,
ice,
liquid_flux,
infiltration_to_micropore,
freeze_thaw_water_flux;
ice_density_factor
)

Update soil macropore water and ice.

Legacy algebra:

```math
\begin{aligned}
W_M' &= W_M + F_l - F_{inf,\mu} + F_{ft} \\
I_M' &= I_M - \frac{F_{ft}}{\rho_i}
\end{aligned}
```

"""
function update_soil_macropore_storage(
water,
ice,
liquid_flux,
infiltration_to_micropore,
freeze_thaw_water_flux;
ice_density_factor
)
new_water =
water +
liquid_flux -
infiltration_to_micropore +
freeze_thaw_water_flux

# CODEX DEBUG BEGIN: remove stray Markdown fence from executable function body.
new_ice =
    ice -
    freeze_thaw_water_flux / ice_density_factor

return (
    water = new_water,
    ice = new_ice
)
# CODEX DEBUG END

end

@doc raw"""
soil_air_storages(
micropore_volume,
macropore_volume,
micropore_water,
micropore_ice,
macropore_water,
macropore_ice
)

Calculate excess storage and air-filled storage for soil micro- and macropores.

Returns:
excess_micropore
air_micropore
excess_macropore
air_macropore
"""
function soil_air_storages(
    micropore_volume,
    macropore_volume,
    micropore_water,
    micropore_ice,
    macropore_water,
    macropore_ice
)
    excess_micropore =
        micropore_volume -
        micropore_water -
        micropore_ice

    air_micropore =
        max(0.0, excess_micropore)

    excess_macropore =
        macropore_volume -
        macropore_water -
        macropore_ice

    air_macropore =
        max(0.0, excess_macropore)

    return (
        excess_micropore = excess_micropore,
        air_micropore = air_micropore,
        excess_macropore = excess_macropore,
        air_macropore = air_macropore
    )
end


@doc raw"""
soil_air_porosity(
air_micropore,
air_macropore,
micropore_volume,
macropore_volume;
tiny = tiny_num2
)

Calculate air-filled porosity relative to total pore + macropore volume.
"""
function soil_air_porosity(
    air_micropore,
    air_macropore,
    micropore_volume,
    macropore_volume;
    tiny = tiny_num2
    )
    pore_volume =
    micropore_volume + macropore_volume

    if pore_volume > tiny
        return max(0.0, (air_micropore + air_macropore) / pore_volume)
    else
        return 0.0
    end
end

@doc raw"""
dynamic_macropore_volume(
base_macropore_volume,
clay_concentration,
micropore_water_fraction,
wilting_point,
soil_layer_volume;
shrink_swell_coeff = FVOLAH
)

Calculate dynamic macropore volume after shrink-swell adjustment.

Legacy pattern:

```math
V_M = \max\left(0, V_{M,0} - F_{VOLAH}C_{clay}(\theta_w - WP)V_{layer}\right)
```

If `shrink_swell_coeff = 0`, this returns the base macropore volume.
"""
function dynamic_macropore_volume(
    base_macropore_volume,
    clay_concentration,
    micropore_water_fraction,
    wilting_point,
    soil_layer_volume;
    shrink_swell_coeff = FVOLAH
    )
    return max(
    0.0,
    base_macropore_volume -
    shrink_swell_coeff *
    clay_concentration *
    (micropore_water_fraction - wilting_point) *
    soil_layer_volume
    )
end

@doc raw"""
macropore_fraction_and_conductivity(
macropore_volume,
base_macropore_volume,
base_macropore_fraction,
base_macropore_conductivity;
tiny = tiny_num2
)

Calculate updated macropore fraction and hydraulic conductivity after dynamic
macropore volume changes.

Legacy patterns:

```math
\begin{aligned}
f_M &= f_{M,0}\frac{V_M}{V_{M,0}} \\
K_M &= K_{M,0}\left(\frac{V_M}{V_{M,0}}\right)^2
\end{aligned}
```

"""
function macropore_fraction_and_conductivity(
    current_macropore_volume,
    reference_macropore_volume,
    reference_macropore_fraction,
    reference_macropore_conductivity;
    tiny = tiny_num2
    )
    if current_macropore_volume > tiny && reference_macropore_volume > tiny
        volume_ratio =
            current_macropore_volume / reference_macropore_volume

        frac_macropore =
            reference_macropore_fraction * volume_ratio

        K_macropore =
            reference_macropore_conductivity * volume_ratio^2

        return (
            frac_macropore = frac_macropore,
            frac_micropore = 1.0 - frac_macropore,
            K_macropore = K_macropore
        )
    else
        return (
            frac_macropore = 0.0,
            frac_micropore = 1.0,
            K_macropore = 0.0
        )
    end
end

@doc raw"""
update_soil_temperature(
old_heat_capacity,
old_temperature,
dry_heat_capacity,
micropore_water,
vapor,
micropore_ice,
macropore_water,
macropore_ice,
net_heat,
fallback_temperature;
cpw,
cpi,
min_heat_capacity = tiny_num
)

Update total soil heat capacity and temperature.

Returns:
heat_capacity
micropore_heat_capacity
macropore_heat_capacity
temperature
"""
function update_soil_temperature(
    old_heat_capacity,
    old_temperature,
    dry_heat_capacity,
    micropore_water,
    vapor,
    micropore_ice,
    macropore_water,
    macropore_ice,
    net_heat,
    fallback_temperature;
    cpw,
    cpi,
    min_heat_capacity = tiny_num
    )
    new_heat_capacity =
        heat_capacity_soil(
            dry_heat_capacity,
            micropore_water,
            vapor,
            micropore_ice,
            macropore_water,
            macropore_ice;
            cpw = cpw,
            cpi = cpi
        )

    micropore_heat_capacity =
        dry_heat_capacity +
        cpw * (micropore_water + vapor) +
        cpi * micropore_ice

    macropore_heat_capacity =
        cpw * macropore_water +
        cpi * macropore_ice

    new_temperature =
        temperature_from_energy(
            old_heat_capacity,
            old_temperature,
            net_heat,
            new_heat_capacity,
            fallback_temperature,
            min_heat_capacity
        )

    return (
        heat_capacity = new_heat_capacity,
        micropore_heat_capacity = micropore_heat_capacity,
        macropore_heat_capacity = macropore_heat_capacity,
        temperature = new_temperature
    )
end
