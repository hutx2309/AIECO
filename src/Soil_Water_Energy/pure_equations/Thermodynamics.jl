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

"""
    VP_at_psi(TK, ψ)

Calculate equilibrium vapor concentration near a water/soil/litter surface,
reduced by water potential.

Equation:
    VP = VP_sat(T) * exp(18 * ψ / (R * T))

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


"""
    latent_heat_evaporation(water_flux, latent_heat)
   
Convert evaporation/condensation water flux to latent heat flux.
Equation:
    H_latent = W * Lv
Inputs:
    W: evaporation (positive) or condensation (negative) flux [m3]
    Lv: latent heat of vaporization [MJ/m3]

Positive or negative sign follows the sign of `water_flux`.
"""
@inline function latent_heat_evaporation(water_flux::Real, latent_heat::Real)
    return water_flux * latent_heat
end


"""
    convective_heat_water(water_flux, TK, cpw)

Heat carried by liquid water flux.

Equation:
    H = cpw * T * W
"""
@inline function convective_heat_water(water_flux::Real, TK::Real, cpw::Real)
    return cpw * TK * water_flux
end


"""
    convective_heat_vapor(vapor_flux, TK, cpw)

Heat carried by vapor flux.

"""
@inline function convective_heat_vapor(vapor_flux::Real, TK::Real, cpw::Real)
    return cpw * TK * vapor_flux
end


"""
    heat_capacity_snow(snow, water, ice, vapor; cps, cpw, cpi)

Snowpack heat capacity from snow, liquid water, ice, and vapor.
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


"""
    heat_capacity_litter(organic_mass, water, vapor, ice; cpo, cpw, cpi)

Surface litter/residue heat capacity.

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


"""
    heat_capacity_soil(dry_heat_capacity, water, vapor, ice, water_macro, ice_macro; cpw, cpi)

Soil heat capacity including dry soil, micropore water/vapor/ice,
and macropore water/ice.
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


"""
    temperature_from_energy(old_heat_capacity, old_temp, net_heat, new_heat_capacity, fallback_temp, min_heat_capacity)

Update temperature from energy balance.

Equation:
    T_new = (C_old * T_old + H_net) / C_new

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


"""
    equilibrium_temperature(C1, T1, C2, T2)

Energy-weighted equilibrium temperature between two connected stores.
"""
@inline function equilibrium_temperature(C1::Real, T1::Real, C2::Real, T2::Real)
    denom = C1 + C2
    return denom > 0.0 ? (C1 * T1 + C2 * T2) / denom : 0.5 * (T1 + T2)
end


"""
    sensible_heat_limited(T_source, T_equil, heat_capacity, substep_fraction)

Maximum sensible heat exchange needed to move source toward equilibrium.
"""
@inline function sensible_heat_limited(
    T_source::Real,
    T_equil::Real,
    heat_capacity::Real,
    substep_fraction::Real
)
    return (T_source - T_equil) * heat_capacity * substep_fraction
end


"""
    conductive_heat_flux(conductance, T_source, T_dest)

Simple temperature-gradient heat flux.

Positive value means heat moves from source to destination if T_source > T_dest.
"""
@inline function conductive_heat_flux(
    conductance::Real,
    T_source::Real,
    T_dest::Real
)
    return conductance * (T_source - T_dest)
end