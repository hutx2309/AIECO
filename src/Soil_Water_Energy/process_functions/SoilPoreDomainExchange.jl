# =============================================================================
# SoilPoreDomainExchange.jl
#
# Process-level function for local water exchange between soil macropores
# and micropores.
#
# Design:
#   - pure function
#   - one soil layer at a time
#   - no mutation
#   - no @unpack
#   - no direct access to waterVar_copy
# =============================================================================


@doc raw"""
    soil_macropore_micropore_exchange_process(...)

Calculate water exchange between macropores and micropores in one soil layer.

This process represents radial water exchange from macropores into the
surrounding micropore/matrix domain.

The potential exchange flux is

```math
F^* =
\frac{
    2\pi K A
    \left(\psi_{sat} - \psi_m\right)
}{
    \ln(r_p / r_m)
}
\Delta t
```

where

* `F^*` is the potential macropore--micropore water exchange flux,
* `K` is hydraulic conductivity,
* `A` is the exchange area,
* `\psi_{sat}` is saturated soil water potential,
* `\psi_m` is matric potential,
* `r_p` is the representative macropore path length,
* `r_m` is the representative macropore radius,
* `\Delta t` is the process time factor.

The legacy implementation uses `6.283` as an approximation to `2\pi`:

```math
2\pi \approx 6.283
```

Positive flux means water moves from macropores into micropores:

```math
F =
\min(F^*, W_M, P_\mu),
\quad F^* > 0
```

where `W_M` is available macropore water and `P_\mu` is available
micropore pore capacity.

Negative flux means water moves from micropores into macropores:

```math
F =
\max(F^*, -P_M, -W_\mu),
\quad F^* < 0
```

where `P_M` is available macropore pore capacity and `W_\mu` is
available micropore water.

# Arguments

* `macropore_water`: current macropore liquid water storage
* `micropore_water`: current micropore liquid water storage
* `macropore_ice`: current macropore ice storage
* `micropore_ice`: current micropore ice storage
* `micropore_water_flux`: net external/internal micropore water flux already calculated
* `macropore_water_flux`: net external/internal macropore water flux already calculated
* `subsurface_water_input`: direct water input to micropores
* `micropore_volume`: micropore pore volume
* `macropore_volume`: macropore pore volume
* `hydraulic_conductivity`: hydraulic conductivity controlling radial exchange
* `exchange_area`: exchange area
* `saturated_potential`: saturated soil water potential
* `matric_potential`: soil matric potential
* `macropore_path_length`: representative macropore path length
* `macropore_radius`: representative macropore radius
* `time_factor`: process time factor

# Returns

A named tuple containing the potential exchange flux, the limited exchange
flux, updated effective water storages/capacities used by the limiter, and
the geometry denominator.
"""
function soil_macropore_micropore_exchange_process(;
macropore_water,
micropore_water,
macropore_ice,
micropore_ice,


micropore_water_flux,
macropore_water_flux,
subsurface_water_input,

micropore_volume,
macropore_volume,

hydraulic_conductivity,
exchange_area,
saturated_potential,
matric_potential,

macropore_path_length,
macropore_radius,

time_factor,

circular_geometry_factor = 6.283,
tiny = tiny_num2


)
# -------------------------------------------------------------------------
# 1. Inactive macropore domain
# -------------------------------------------------------------------------

if macropore_water <= tiny
    return (
        potential_flux = 0.0,
        limited_flux = 0.0,
        micropore_water_after_prior_fluxes = 0.0,
        macropore_water_after_prior_fluxes = 0.0,
        micropore_capacity = 0.0,
        macropore_capacity = 0.0,
        geometry_denominator = 0.0
    )
end


# -------------------------------------------------------------------------
# 2. Potential radial exchange flux
# -------------------------------------------------------------------------

geometry_denominator =
    log(macropore_path_length / macropore_radius)

potential_flux =
    circular_geometry_factor *
    hydraulic_conductivity *
    exchange_area *
    (saturated_potential - matric_potential) /
    geometry_denominator *
    time_factor


# -------------------------------------------------------------------------
# 3. Effective storages after previously calculated water fluxes
# -------------------------------------------------------------------------

micropore_water_after_prior_fluxes =
    micropore_water +
    micropore_water_flux +
    subsurface_water_input

micropore_capacity =
    max(
        0.0,
        micropore_volume -
        micropore_water_after_prior_fluxes -
        micropore_ice
    )

macropore_water_after_prior_fluxes =
    macropore_water +
    macropore_water_flux

macropore_capacity =
    max(
        0.0,
        macropore_volume -
        macropore_water_after_prior_fluxes -
        macropore_ice
    )


# -------------------------------------------------------------------------
# 4. Storage-limited exchange flux
# -------------------------------------------------------------------------

if potential_flux > 0.0
    limited_flux =
        max(
            0.0,
            min(
                potential_flux,
                macropore_water_after_prior_fluxes,
                micropore_capacity
            )
        )
else
    limited_flux =
        min(
            0.0,
            max(
                potential_flux,
                -macropore_capacity,
                -micropore_water_after_prior_fluxes
            )
        )
end


# -------------------------------------------------------------------------
# 5. Return process diagnostics
# -------------------------------------------------------------------------

return (
    potential_flux =
        potential_flux,

    limited_flux =
        limited_flux,

    micropore_water_after_prior_fluxes =
        micropore_water_after_prior_fluxes,

    macropore_water_after_prior_fluxes =
        macropore_water_after_prior_fluxes,

    micropore_capacity =
        micropore_capacity,

    macropore_capacity =
        macropore_capacity,

    geometry_denominator =
        geometry_denominator
)

end
