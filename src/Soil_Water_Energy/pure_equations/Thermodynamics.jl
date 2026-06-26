# ============================================================
# Thermodynamics.jl
#
# Pure thermodynamic helper functions for water-energy exchange.
# These functions should:
#   - not mutate model state
#   - not access water variables directly
#   - not unpack structs
#   - return values only
# ============================================================
@doc raw"""
    cal_SatVP_conc(T)

Calculate saturated vapor concentration from temperature.

```math
VP_{sat}(T) =
\frac{2.173\times10^{-3}}{T}\,0.61\,
\exp\left[5360\left(3.661\times10^{-3}-\frac{1}{T}\right)\right]
```

`VP_{sat}` is the saturation vapor concentration in water-equivalent volume
per air volume. The coefficient `0.61` is the saturation vapor pressure in kPa
at 0 deg C, and the exponential term is the Clausius-Clapeyron temperature
response.
"""


@inline function cal_SatVP_conc(T)
    VP_sat=2.173E-03/T*0.61*exp(5360.0*(3.661E-03-1.0/T))
    return VP_sat
end




@doc raw"""
    VP_at_psi(TK, ψ)

Calculate equilibrium vapor concentration near a water/soil/litter surface,
reduced by water potential.

Equation:

```math
VP(T, \psi) = VP_{sat}(T)\,\exp\left(\frac{18\psi}{RT}\right)
```

Inputs:
    TK : temperature [K]
    ψ  : water potential [MPa]
    R : gas constant [8.3143 J/(mol*K)]
    18: J/mol/MPa conversion factor
Output:
    vapor concentration (m3 water-equivalent vapor /m3 air = g/m3)
"""
@inline function VP_at_psi(TK::Real, ψ::Real)
    return cal_SatVP_conc(TK) * exp(18.0 * ψ / (8.3143 * TK))
end


@doc raw"""
    latent_heat_evaporation(water_flux, latent_heat)

Convert evaporation/condensation water flux to latent heat flux.
Equation:

```math
H_{latent} = W L_v
```
Inputs:
    W: evaporation (positive) or condensation (negative) flux [m3]
    Lv: latent heat of vaporization [MJ/m3]

Positive or negative sign follows the sign of `water_flux`.
"""
@inline function latent_heat_evaporation(water_flux::Real, latent_heat::Real)
    return water_flux * latent_heat
end


@doc raw"""
    advective_heat_water(water_flux, TK, cpw)

Heat carried by liquid water flux.

Equation:

```math
H = c_{pw} T W
```
"""
@inline function advective_heat_water(water_flux::Real, TK::Real, cpw::Real)
    return cpw * TK * water_flux
end


@doc raw"""
    advective_heat_vapor(vapor_flux, TK, cpw)

Heat carried by moving vapor flux.

Inputs:
    vapor_flux : vapor water-equivalent volume/depth for current substep
    TK         : temperature [K]
    cpw        : water-equivalent heat capacity [MJ m⁻³ K⁻¹]

Output:
    heat carried by vapor flux [MJ] or [MJ m⁻²], consistent with vapor_flux

Note:
    This is not latent heat of evaporation. It is the sensible/enthalpy
    content transported by vapor movement.
"""
@inline function advective_heat_vapor(vapor_flux::Real, TK::Real, cpw::Real)
    return cpw * TK * vapor_flux
end


@doc raw"""
    heat_capacity_snow(snow, water, ice, vapor; cps, cpw, cpi)

Snowpack heat capacity from snow, liquid water, ice, and vapor.

```math
C_{snow} = c_{ps}S + c_{pw}(W + V) + c_{pi}I
```
"""
@inline function heat_capacity_snow(
    snow::Real,
    water::Real,
    ice::Real,
    vapor::Real;
    cps::Real,
    cpw::Real,
    cpi::Real
)
    return cps * snow + cpw * (water + vapor) + cpi * ice
end


@doc raw"""
    heat_capacity_litter(organic_mass, water, vapor, ice; cpo, cpw, cpi)

Surface litter/residue heat capacity.

```math
C_{litter} = c_{po}M_o + c_{pw}(W + V) + c_{pi}I
```

`organic_mass` can be surf_SOC + surf_charcoal or another dry organic pool.
"""
@inline function heat_capacity_litter(
    organic_mass::Real,
    water::Real,
    vapor::Real,
    ice::Real;
    cpo::Real,
    cpw::Real,
    cpi::Real
)
    return cpo * organic_mass + cpw * (water + vapor) + cpi * ice
end


@doc raw"""
    heat_capacity_soil(dry_heat_capacity, water, vapor, ice, water_macro, ice_macro; cpw, cpi)

Soil heat capacity including dry soil, micropore water/vapor/ice,
and macropore water/ice.

```math
C_{soil} = C_{dry} + c_{pw}(W_\mu + V_\mu + W_M) + c_{pi}(I_\mu + I_M)
```
"""
@inline function heat_capacity_soil(
    dry_heat_capacity::Real,
    water::Real,
    vapor::Real,
    ice::Real,
    water_macro::Real,
    ice_macro::Real;
    cpw::Real,
    cpi::Real
)
    return dry_heat_capacity +
           cpw * (water + vapor + water_macro) +
           cpi * (ice + ice_macro)
end


@doc raw"""
    temperature_from_energy(old_heat_capacity, old_temp, net_heat, new_heat_capacity, fallback_temp, min_heat_capacity)

Update temperature from energy balance.

Equation:

```math
T_{new} = \frac{C_{old}T_{old} + H_{net}}{C_{new}}
```

If new heat capacity is too small, return fallback temperature.
"""
@inline function temperature_from_energy(
    old_heat_capacity::Real,
    old_temp::Real,
    net_heat::Real,
    new_heat_capacity::Real,
    fallback_temp::Real,
    min_heat_capacity::Real
)
    if new_heat_capacity > min_heat_capacity
        return (old_heat_capacity * old_temp + net_heat) / new_heat_capacity
    else
        return fallback_temp
    end
end


@doc raw"""
    equilibrium_temperature(C1, T1, C2, T2)

Energy-weighted equilibrium temperature between two connected stores.

```math
T_{eq} = \frac{C_1T_1 + C_2T_2}{C_1 + C_2}
```
"""
@inline function equilibrium_temperature(C1::Real, T1::Real, C2::Real, T2::Real)
    denom = C1 + C2
    return denom > 0.0 ? (C1 * T1 + C2 * T2) / denom : 0.5 * (T1 + T2)
end


@doc raw"""
    sensible_heat_limited(T_source, T_equil, heat_capacity, substep_fraction)

Maximum sensible heat exchange needed to move source toward equilibrium.

```math
H_{sens} = (T_{source} - T_{eq}) C\, f_{step}
```
"""
@inline function sensible_heat_limited(
    T_source::Real,
    T_equil::Real,
    heat_capacity::Real,
    substep_fraction::Real
)
    return (T_source - T_equil) * heat_capacity * substep_fraction
end


@doc raw"""
    conductive_heat_flux(conductance, T_source, T_dest)

Simple temperature-gradient heat flux.

```math
H = G(T_{source} - T_{dest})
```

Positive value means heat moves from source to destination if `T_source > T_dest`.
"""
@inline function conductive_heat_flux(
    conductance::Real,
    T_source::Real,
    T_dest::Real
)
    return conductance * (T_source - T_dest)
end
