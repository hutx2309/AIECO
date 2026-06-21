# =============================================================================
# SurfaceExchange.jl
#
# Pure helper equations for surface cover fractions, surface resistances,
# and exchange conductances used by the water-energy module.
#
# These functions should not mutate model state.
# They should not unpack model structs.
# They should only express small physical / numerical equations.
# =============================================================================


# -----------------------------------------------------------------------------
# Surface cover fractions
# -----------------------------------------------------------------------------

"""
    snow_cover_fraction(snow_depth, reference_depth)

Calculate fractional snow cover from snow depth.

Legacy equation:
    f_snow = min(1, sqrt(max(0, snow_depth) / reference_depth))

Inputs:
    snow_depth      : snow depth
    reference_depth : reference snow depth for full cover, e.g. DPTHSX

Output:
    snow cover fraction in [0, 1]
"""
function snow_cover_fraction(snow_depth, reference_depth)
    if reference_depth > 0.0
        return clamp(sqrt(max(0.0, snow_depth) / reference_depth), 0.0, 1.0)
    else
        return snow_depth > 0.0 ? 1.0 : 0.0
    end
end


"""
    snow_free_fraction(frac_snow_cover)

Calculate snow-free surface fraction.
"""
snow_free_fraction(frac_snow_cover) =
    1.0 - clamp(frac_snow_cover, 0.0, 1.0)


"""
    bare_ground_fraction(surface_litter_carbon, surface_area; coefficient=0.5e-2)

Calculate bare-ground fraction from surface litter organic C + charcoal.

Legacy equation:
    frac_bare = min(1, max(0, exp(-0.5E-02 * litter_mass / area)))

Inputs:
    surface_litter_carbon : surf_SOC + surf_charcoal
    surface_area          : surface area used by legacy code
    coefficient           : empirical extinction coefficient

Output:
    bare-ground fraction in [0, 1]
"""
function bare_ground_fraction(
    surface_litter_carbon,
    surface_area;
    coefficient = 0.5e-2
)
    if surface_area > 0.0
        return clamp(exp(-coefficient * surface_litter_carbon / surface_area), 0.0, 1.0)
    else
        return 1.0
    end
end


"""
    water_cover_fraction_from_excess(excess_water_ice, holding_capacity; tiny=tiny_num)

Convert excess surface water/ice into a bounded surface-cover fraction.

The raw ratio:

    excess_water_ice / holding_capacity

can be larger than 1.0, but the surface-cover fraction cannot exceed 1.0.

Inputs:
    excess_water_ice : surface water + ice above litter holding capacity
    holding_capacity : surface/litter water holding capacity

Output:
    water-cover fraction contribution in [0, 1]
"""
function water_cover_fraction_from_excess(
    excess_water_ice,
    holding_capacity;
    tiny = tiny_num
)
    if holding_capacity > tiny
        return clamp(excess_water_ice / holding_capacity, 0.0, 1.0)
    else
        return 0.0
    end
end


"""
    water_free_fraction(frac_bare_ground, water_cover_from_excess)

Calculate water-free bare-ground fraction.

Legacy equation:
    frac_waterFree = max(0, frac_bareGrnd - water_cover_ratio)
"""
function water_free_fraction(frac_bare_ground, water_cover_from_excess)
    return max(0.0, frac_bare_ground - clamp(water_cover_from_excess, 0.0, 1.0))
end


"""
    surface_fraction_partition(
        snow_depth,
        reference_snow_depth,
        frac_bare_ground,
        excess_surface_water_ice,
        surface_water_holding_capacity;
        tiny=tiny_num
    )

Diagnose snow-covered, snow-free, water-free, and water-covered surface fractions.

This groups the repeated WF1/WF2 logic:

    frac_snowCover = snow_cover_fraction(...)
    frac_snowFree  = 1 - frac_snowCover

    water_cover_from_excess = clamp(excess / holding_capacity, 0, 1)
    frac_waterFree  = max(0, frac_bare_ground - water_cover_from_excess)
    frac_waterCover = 1 - frac_waterFree

Inputs:
    snow_depth                     : current snow depth
    reference_snow_depth           : reference depth for full snow cover, e.g. DPTHSX
    frac_bare_ground               : bare-ground fraction from litter/residue cover
    excess_surface_water_ice       : excess water + ice on surface
    surface_water_holding_capacity : surface/litter water holding capacity

Output:
    NamedTuple with:
        frac_snowCover
        frac_snowFree
        frac_waterFree
        frac_waterCover
        water_cover_from_excess
"""
function surface_fraction_partition(
    snow_depth,
    reference_snow_depth,
    frac_bare_ground,
    excess_surface_water_ice,
    surface_water_holding_capacity;
    tiny = tiny_num
)
    frac_snowCover =
        snow_cover_fraction(snow_depth, reference_snow_depth)

    frac_snowFree =
        1.0 - frac_snowCover

    water_cover_from_excess =
        water_cover_fraction_from_excess(
            excess_surface_water_ice,
            surface_water_holding_capacity;
            tiny = tiny
        )

    frac_waterFree =
        water_free_fraction(frac_bare_ground, water_cover_from_excess)

    frac_waterCover =
        1.0 - frac_waterFree

    return (
        frac_snowCover = frac_snowCover,
        frac_snowFree = frac_snowFree,
        frac_waterFree = frac_waterFree,
        frac_waterCover = frac_waterCover,
        water_cover_from_excess = water_cover_from_excess
    )
end

# -----------------------------------------------------------------------------
# Resistance equations
# -----------------------------------------------------------------------------

"""
    bounded_resistance(raw_resistance, lower_bound, upper_bound)

Constrain a resistance value between lower and upper bounds.

Useful for aerodynamic resistance limits such as:
    min(r_max, max(r_min, raw_resistance))
"""
function bounded_resistance(raw_resistance, lower_bound, upper_bound)
    return clamp(raw_resistance, lower_bound, upper_bound)
end


"""
    stability_corrected_resistance(base_resistance, richardson_factor; tiny=tiny_num)

Apply a Richardson/stability correction to a resistance.

Legacy-style equation:
    resistance = base_resistance / R

where R is a Richardson/stability correction factor.
"""
function stability_corrected_resistance(
    base_resistance,
    richardson_factor;
    tiny = tiny_num
)
    if abs(richardson_factor) > tiny
        return base_resistance / richardson_factor
    else
        return base_resistance / tiny
    end
end

"""
    soil_surface_resistance(frac_water_free, r_soil, frac_water_cover, r_litter_water; tiny=tiny_num)

Effective soil/litter surface resistance used for snow-free ground exchange.

Legacy WF2 equation:
    r_effective_surf = 1 / (frac_waterFree/r_grndAir2surf +
                frac_waterCover/(r_grndAir2surf + r_littBondryLimit_L))

Inputs:
    frac_water_free            : water-free bare-ground fraction
    soil_boundary_resistance   : soil-side boundary resistance, e.g. r_grndAir2surf
    frac_water_cover           : water/litter-covered fraction
    litter_limited_resistance  : water/litter pathway resistance, e.g. r_grndAir2surf + r_littBondryLimit_L

Output:
    effective surface resistance r_effective_surf
"""
function soil_surface_resistance(frac_water_free, soil_boundary_resistance,
                                 frac_water_cover, litter_limited_resistance;
                                 tiny = tiny_num)
    conductance =
        frac_water_free / soil_boundary_resistance +
        frac_water_cover / (soil_boundary_resistance + litter_limited_resistance)

    return conductance > tiny ? 1.0 / conductance : 1.0 / tiny
end


# -----------------------------------------------------------------------------
# Exchange conductance equations
# -----------------------------------------------------------------------------
"""
    vapor_exchange_coefficient(
        area,
        time_factor,
        cover_fraction,
        resistance;
        multiplier = 1.0,
        tiny = tiny_num
    )

Integrated vapor-exchange coefficient for one model substep.

Physical form:

    E = (A * Δt * f * m / r_v) * (C_air - C_surface)

where

    E          = vapor/water exchange amount during the substep
    A          = exchange area
    Δt         = substep time factor
    f          = active surface fraction
    m          = optional multiplier, e.g. dt_snow or dt_litt
    r_v        = vapor-transfer resistance
    C_air      = vapor concentration in air
    C_surface  = equilibrium vapor concentration at the surface

Therefore this function returns only the coefficient:

    K_v = A * Δt * f * m / r_v

and the final vapor exchange is:

    E = K_v * (C_air - C_surface)

Sign convention:
    positive E means vapor moves from air to surface;
    negative E means evaporation/sublimation from surface to air.
"""
function vapor_exchange_coefficient(
    area,
    time_factor,
    cover_fraction,
    resistance;
    multiplier = 1.0,
    tiny = tiny_num
)
    if resistance > tiny
        return area * time_factor * cover_fraction * multiplier / resistance
    else
        return 0.0
    end
end
 
"""
    sensible_heat_exchange_coefficient(
        area,
        time_factor,
        cover_fraction,
        resistance;
        air_volumetric_heat_capacity = 1.25e-3,
        multiplier = 1.0,
        tiny = tiny_num
    )

Integrated sensible-heat exchange coefficient for one model substep.

Physical form:

    H = (A * Δt * f * m * ρ_air * c_p_air / r_h) *
        (T_air - T_surface)

where

    H          = sensible heat exchanged during the substep
    A          = exchange area
    Δt         = substep time factor
    f          = active surface fraction
    m          = optional multiplier, e.g. dt_snow or dt_litt
    ρ_air c_p  = volumetric heat capacity of air
    r_h        = heat-transfer resistance
    T_air      = air temperature
    T_surface  = surface temperature

Therefore this function returns only the coefficient:

    K_h = A * Δt * f * m * ρ_air * c_p_air / r_h

and the final sensible heat exchange is:

    H = K_h * (T_air - T_surface)

The default air_volumetric_heat_capacity = 1.25e-3 is the legacy model value,
approximately 0.00125 MJ m⁻³ K⁻¹.
"""
function sensible_heat_exchange_coefficient(
    area,
    time_factor,
    cover_fraction,
    resistance;
    air_volumetric_heat_capacity = 1.25e-3,
    multiplier = 1.0,
    tiny = tiny_num
)
    if resistance > tiny
        return area *
               time_factor *
               cover_fraction *
               multiplier *
               air_volumetric_heat_capacity / resistance
    else
        return 0.0
    end
end

 