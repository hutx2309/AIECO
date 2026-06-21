# =============================================================================

# Radiation.jl

# Pure radiation helper equations for water-energy exchange.

# Design rules:

# - no mutation

# - no direct access to waterVar_copy

# - no process-level state update

# - functions return algebraic radiation quantities only

# =============================================================================

# -----------------------------------------------------------------------------
# 1. Radiation partitioning
# -----------------------------------------------------------------------------

"""
partition_ground_radiation(
incoming_radiation,
surface_fraction;
multiplier = 1.0
)

Partition incoming radiation to a surface fraction.

Physical form:

```
R_surface = R_in * f_surface * m
```

where

```
R_surface = radiation assigned to the surface
R_in      = incoming shortwave or longwave radiation
f_surface = active surface fraction
m         = optional multiplier, e.g. dt_snow or dt_litt
```

This helper is appropriate for partitioning incoming shortwave or incoming
longwave radiation among snow, bare soil, and litter/residue surfaces.
"""
function partition_ground_radiation(
    shortwave_ground,
    longwave_ground,
    frac_snowCover,
    frac_snowFree,
    frac_bareGrnd,
    frac_littCover;
    snow_multiplier = 1.0,
    litter_multiplier = 1.0
    )
    snow_fraction = frac_snowCover
    soil_fraction = frac_snowFree * frac_bareGrnd
    litt_fraction = frac_snowFree * frac_littCover

    return (
        Rad_SWsnow = shortwave_ground * snow_fraction * snow_multiplier,
        Rad_SWsoil = shortwave_ground * soil_fraction,
        Rad_SWlitt = shortwave_ground * litt_fraction * litter_multiplier,
        Rad_LW2snow = longwave_ground * snow_fraction * snow_multiplier,
        Rad_LW2soil = longwave_ground * soil_fraction,
        Rad_LW2litt = longwave_ground * litt_fraction * litter_multiplier
    )
end

# -----------------------------------------------------------------------------

# 2. Longwave emission coefficient

# -----------------------------------------------------------------------------

"""
    surface_longwave_coefficients(
        area,
        frac_snowCover,
        frac_snowFree,
        frac_bareGrnd,
        frac_littCover,
        snow_time_factor,
        soil_time_factor,
        litter_time_factor,
        frac_grndRad,
        snow_emissivity,
        soil_emissivity,
        litter_emissivity;
        stefan_boltzmann_time = 2.04e-10
    )

Calculate area- and time-integrated longwave radiation coefficients for
snow, soil, and litter surfaces.

Physical form:

    LW_out = ε * σ* * A * f_surface * Δt * T_surface^4

The returned `*_L` coefficients include `frac_grndRad`, representing the
open-ground/sky radiation fraction.

The returned `*Canpy_L` coefficients do not include `frac_grndRad`; they are
used later for canopy/standing-dead longwave exchange:

    LW_exchange = K * (T_canopy^4 - T_surface^4) * view_fraction
"""
function surface_longwave_coefficients(
    area,
    frac_snowCover,
    frac_snowFree,
    frac_bareGrnd,
    frac_littCover,
    snow_time_factor,
    soil_time_factor,
    litter_time_factor,
    frac_grndRad,
    snow_emissivity,
    soil_emissivity,
    litter_emissivity;
    stefan_boltzmann_time = 2.04e-10
    )
    snow_factor = area * frac_snowCover * snow_time_factor
    soil_factor = area * frac_snowFree * frac_bareGrnd * soil_time_factor
    litt_factor = area * frac_snowFree * frac_littCover * litter_time_factor

    return (
        RadLW_fromSnow_L =
            snow_emissivity * stefan_boltzmann_time *
            snow_factor * frac_grndRad,

        RadLW_fromSnowCanpy_L =
            snow_emissivity * stefan_boltzmann_time *
            snow_factor,

        RadLW_fromSoil_L =
            soil_emissivity * stefan_boltzmann_time *
            soil_factor * frac_grndRad,

        RadLW_fromSoilCanpy_L =
            soil_emissivity * stefan_boltzmann_time *
            soil_factor,

        RadLW_fromLitt_L =
            litter_emissivity * stefan_boltzmann_time *
            litt_factor * frac_grndRad,

        RadLW_fromLittCanpy_L =
            litter_emissivity * stefan_boltzmann_time *
            litt_factor
    )
end

# -----------------------------------------------------------------------------

# 3. Longwave temperature exchange

# -----------------------------------------------------------------------------

"""
longwave_temperature_exchange(
coefficient,
source_temperature,
surface_temperature,
view_fraction
)

Longwave radiation exchange between a radiating source and a receiving surface.

Physical form:

```
LW_exchange = K_LW * (T_source^4 - T_surface^4) * f_view
```

where

```
LW_exchange = net longwave energy received by the surface
K_LW        = longwave exchange coefficient
T_source    = source temperature in K, e.g. canopy or standing dead
T_surface   = receiving surface temperature in K
f_view      = view/radiation partition fraction
```

Positive value means the receiving surface gains longwave energy.
Negative value means the receiving surface loses longwave energy to the source.
"""
function longwave_temperature_exchange(
    coefficient,
    source_temperature,
    surface_temperature,
    view_fraction
    )
    return coefficient *
    (source_temperature^4 - surface_temperature^4) *
    view_fraction
end

# -----------------------------------------------------------------------------

# 5. Surface albedo

# -----------------------------------------------------------------------------
"""
    weighted_surface_albedo(
        component_amounts,
        component_albedos;
        fallback_albedo,
        tiny = tiny_num
    )

Weighted-average albedo for a mixed surface.

Physical/empirical form:

    α_surface = Σ(αᵢ * amountᵢ) / Σ(amountᵢ)

where each component can be dry soil, dry litter, snow, ice, or liquid water.

If the total amount is too small, return `fallback_albedo`.
"""
function weighted_surface_albedo(
    component_amounts,
    component_albedos;
    fallback_albedo,
    tiny = tiny_num
)
    total_amount = sum(component_amounts)

    if total_amount > tiny
        return sum(component_amounts .* component_albedos) / total_amount
    else
        return fallback_albedo
    end
end

# CODEX DEBUG BEGIN: route snow, soil, and litter albedo calculations through the shared weighted helper.
snow_surface_albedo(
    snow_volume,
    ice_volume,
    water_volume;
    tiny = tiny_num
) = weighted_surface_albedo(
    [snow_volume, ice_volume, water_volume],
    [0.85, 0.30, 0.06];
    fallback_albedo = 0.30,
    tiny = tiny,
)

soil_surface_albedo(
    dry_soil_albedo,
    dry_soil_amount,
    water_volume,
    ice_volume;
    tiny = tiny_num
) = weighted_surface_albedo(
    [dry_soil_amount, water_volume, ice_volume],
    [dry_soil_albedo, 0.06, 0.30];
    fallback_albedo = dry_soil_albedo,
    tiny = tiny,
)

litter_surface_albedo(
    dry_litter_albedo,
    dry_litter_amount,
    water_volume,
    ice_volume;
    tiny = tiny_num
) = weighted_surface_albedo(
    [dry_litter_amount, water_volume, ice_volume],
    [dry_litter_albedo, 0.06, 0.30];
    fallback_albedo = dry_litter_albedo,
    tiny = tiny,
)
# CODEX DEBUG END
