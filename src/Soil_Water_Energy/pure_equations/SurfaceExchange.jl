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

@doc raw"""
    snow_cover_fraction(snow_depth, reference_depth)

Calculate fractional snow cover from snow depth.

Legacy equation:

```math
f_{snow} = \min\left(1, \sqrt{\frac{\max(0, d_{snow})}{d_{ref}}}\right)
```

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


@doc raw"""
    snow_free_fraction(frac_snow_cover)

Calculate snow-free surface fraction.
"""
snow_free_fraction(frac_snow_cover) =
    1.0 - clamp(frac_snow_cover, 0.0, 1.0)


@doc raw"""
    bare_ground_fraction(surface_litter_carbon, surface_area; coefficient=0.5e-2)

Calculate bare-ground fraction from surface litter organic C + charcoal.

Legacy equation:

```math
f_{bare} = \min\left(1, \max\left(0, \exp\left[-c\frac{M_{litter}}{A}\right]\right)\right)
```

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


@doc raw"""
    water_cover_fraction_from_excess(excess_water_ice, holding_capacity; tiny=tiny_num)

Convert excess surface water/ice into a bounded surface-cover fraction.

The raw ratio is

```math
r_w = \frac{W_{excess}}{W_{hold}}
```

which can be larger than 1.0, but the surface-cover fraction cannot exceed
1.0.

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


@doc raw"""
    water_free_fraction(frac_bare_ground, water_cover_from_excess)

Calculate water-free bare-ground fraction.

Legacy equation:

```math
f_{waterfree} = \max(0, f_{bare} - f_{water})
```
"""
function water_free_fraction(frac_bare_ground, water_cover_from_excess)
    return max(0.0, frac_bare_ground - clamp(water_cover_from_excess, 0.0, 1.0))
end


@doc raw"""
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

```math
\begin{aligned}
f_{snowfree} &= 1 - f_{snow} \\
f_{water} &= \operatorname{clamp}\left(\frac{W_{excess}}{W_{hold}}, 0, 1\right) \\
f_{waterfree} &= \max(0, f_{bare} - f_{water}) \\
f_{watercover} &= 1 - f_{waterfree}
\end{aligned}
```

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

@doc raw"""
    bounded_resistance(raw_resistance, lower_bound, upper_bound)

Constrain a resistance value between lower and upper bounds.

Useful for aerodynamic resistance limits such as:

```math
r = \min(r_{max}, \max(r_{min}, r_{raw}))
```
"""
function bounded_resistance(raw_resistance, lower_bound, upper_bound)
    return clamp(raw_resistance, lower_bound, upper_bound)
end


@doc raw"""
    stability_corrected_resistance(base_resistance, richardson_factor; tiny=tiny_num)

Apply a Richardson/stability correction to a resistance.

Legacy-style equation:

```math
r = \frac{r_0}{R}
```

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

@doc raw"""
    soil_surface_resistance(frac_water_free, r_soil, frac_water_cover, r_litter_water; tiny=tiny_num)

Effective soil/litter surface resistance used for snow-free ground exchange.

Legacy WF2 equation:

```math
r_{eff} = \left(\frac{f_{waterfree}}{r_s} +
\frac{f_{watercover}}{r_s + r_l}\right)^{-1}
```

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
@doc raw"""
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

```math
E = \frac{A\Delta t\,f\,m}{r_v}(C_{air} - C_{surface})
```

where `E` is vapor or water exchange during the substep, `A` is exchange area,
`\Delta t` is the substep time factor, `f` is active surface fraction, `m` is
an optional multiplier, `r_v` is vapor-transfer resistance, and `C_{air} -
C_{surface}` is the vapor concentration gradient.

Therefore this function returns only the coefficient:

```math
K_v = \frac{A\Delta t\,f\,m}{r_v}
```

and the final vapor exchange is:

```math
E = K_v(C_{air} - C_{surface})
```

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

@doc raw"""
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

```math
H = \frac{A\Delta t\,f\,m\,\rho_{air}c_{p,air}}{r_h}(T_{air} - T_{surface})
```

where `H` is sensible heat exchanged during the substep, `A` is exchange area,
`\Delta t` is the substep time factor, `f` is active surface fraction, `m` is
an optional multiplier, `\rho_{air}c_{p,air}` is volumetric heat capacity of
air, and `r_h` is heat-transfer resistance.

Therefore this function returns only the coefficient:

```math
K_h = \frac{A\Delta t\,f\,m\,\rho_{air}c_{p,air}}{r_h}
```

and the final sensible heat exchange is:

```math
H = K_h(T_{air} - T_{surface})
```

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


