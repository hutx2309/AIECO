# =============================================================================
# SurfaceLitterPhaseChange.jl
#
# Process-level function for vapor-liquid exchange and freeze-thaw in the
# surface litter micropore domain.
#
# Design:
#   - pure function
#   - one surface/litter state at a time
#   - no mutation
#   - no @unpack
#   - no direct access to waterVar_copy or fluxVar_4Acc
# =============================================================================


@doc raw"""
    surface_litter_phase_change_process(...)

Calculate vapor-liquid exchange and freeze-thaw in the surface litter
micropore domain.
 
The nonnegative surface litter storages are

```math
W_l = \max(0, W)
```

```math
V_l = \max(0, V)
```

```math
I_l = \max(0, I)
```

where ``W`` is liquid water, ``V`` is vapor storage, and ``I`` is ice.

The free litter pore volume is

```math
P_l =
\max(0, P - W_l - I_l)
```

where ``P`` is total surface litter micropore volume.

When free pore volume is available, the equilibrium vapor concentration is

```math
C_v^* =
C_v(T, \psi)
```

implemented by

```julia
VP_at_psi(surface_temperature, surface_matric_potential)
```

The potential vapor-liquid exchange is

```math
F_v^* =
V_l - C_v^* P_l
```

where positive flux means condensation to liquid water and negative flux means
evaporation from liquid water.

The vapor-liquid exchange is limited by available liquid water:

```math
F_v =
\max(F_v^*, -W_l f_r)
```

The associated latent heat flux is

```math
H_v =
L_v F_v
```

The litter freezing temperature is calculated from water potential:

```math
T_f =
\frac{-9.0959 \times 10^4}{\psi - L_f}
```

Freeze-thaw is active only when the temperature and available phase storage
allow freezing or thawing.

The potential freeze-thaw latent heat is

```math
H_f^* =
C_l
\frac{
    T_f - T_l
}{
    (1 + a_T T_f)(1 - a_\psi \psi)
}
f_r
```

The limited freeze-thaw heat and liquid-water-equivalent flux are calculated by

```julia
freeze_thaw_limited_fluxes(...)
```

# Sign convention

For vapor-liquid exchange:

- `vapor_liquid_flux > 0`: condensation to litter liquid water.
- `vapor_liquid_flux < 0`: evaporation from litter liquid water.

For freeze-thaw water flux:

- `freeze_thaw_water_flux > 0`: ice melts to liquid water.
- `freeze_thaw_water_flux < 0`: liquid water freezes to ice.

# Legacy behavior note

Current WF6 zeroes both litter evaporation and litter freeze-thaw fluxes when
the freeze-thaw activation condition is false. This is physically surprising,
because evaporation-condensation and freeze-thaw could be independent. However,
to preserve the current WF6 behavior, this function defaults to

```julia
zero_vapor_exchange_when_freeze_thaw_inactive = true
```

Set this to `false` only if you intentionally want to decouple litter vapor
exchange from freeze-thaw activation.

# Returns

A named tuple with litter storages, free pore volume, equilibrium vapor
concentration, vapor-liquid flux, latent heat, freezing temperature,
freeze-thaw flux, and freeze-thaw heat.
"""
function surface_litter_phase_change_process(;
    # -------------------------------------------------------------------------
    # Surface litter state
    # -------------------------------------------------------------------------
    liquid_water,
    vapor_storage,
    ice_storage,
    micropore_volume,

    surface_temperature,
    surface_matric_potential,
    surface_heat_capacity,
    surface_layer_volume,

    # -------------------------------------------------------------------------
    # Time and thermodynamic constants
    # -------------------------------------------------------------------------
    removal_fraction,
    latent_heat_vaporization,
    latent_heat_fusion = 333.0,
    ice_density_factor = soil_iceDensty,

    # -------------------------------------------------------------------------
    # Freeze-thaw coefficients
    # -------------------------------------------------------------------------
    thermal_slope_coeff = 6.2913e-3,
    water_potential_coeff = 0.10,

    # -------------------------------------------------------------------------
    # Legacy behavior switch
    # -------------------------------------------------------------------------
    zero_vapor_exchange_when_freeze_thaw_inactive = true,

    # -------------------------------------------------------------------------
    # Numerical controls
    # -------------------------------------------------------------------------
    tiny = tiny_num2
)
    # -------------------------------------------------------------------------
    # 1. Nonnegative litter storages and free pore volume
    # -------------------------------------------------------------------------

    liquid_water_litter =
        max(0.0, liquid_water)

    vapor_storage_litter =
        max(0.0, vapor_storage)

    ice_storage_litter =
        max(0.0, ice_storage)

    free_pore_volume =
        max(
            0.0,
            micropore_volume -
            liquid_water_litter -
            ice_storage_litter
        )


    # -------------------------------------------------------------------------
    # 2. Vapor-liquid exchange in surface litter
    # -------------------------------------------------------------------------

    if free_pore_volume > tiny
        equilibrium_vapor_concentration =
            VP_at_psi(
                surface_temperature,
                surface_matric_potential
            )

        potential_vapor_liquid_flux =
            pore_vapor_condensation_potential(
                vapor_storage_litter,
                equilibrium_vapor_concentration,
                free_pore_volume
            )

        vapor_liquid_flux =
            limit_evaporation_by_liquid(
                potential_vapor_liquid_flux,
                liquid_water_litter,
                removal_fraction
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
    # 3. Freezing temperature from litter water potential
    # -------------------------------------------------------------------------

    freezing_temperature =
        freezing_temperature_from_water_potential(
            surface_matric_potential;
            latent_heat_fusion = latent_heat_fusion
        )


    # -------------------------------------------------------------------------
    # 4. Freeze-thaw activation
    # -------------------------------------------------------------------------

    freeze_thaw_active =
        freeze_thaw_is_active(
            surface_temperature,
            freezing_temperature,
            liquid_water_litter,
            ice_storage_litter,
            surface_layer_volume
        )


    # -------------------------------------------------------------------------
    # 5. Freeze-thaw heat and water fluxes
    # -------------------------------------------------------------------------

    if freeze_thaw_active
        potential_freeze_thaw_heat =
            freeze_thaw_heat_potential(
                surface_heat_capacity,
                surface_temperature,
                freezing_temperature,
                surface_matric_potential,
                removal_fraction;
                thermal_slope_coeff = thermal_slope_coeff,
                water_potential_coeff = water_potential_coeff
            )

        freeze_thaw =
            freeze_thaw_limited_fluxes(
                potential_freeze_thaw_heat,
                liquid_water_litter,
                ice_storage_litter,
                removal_fraction;
                latent_heat_fusion = latent_heat_fusion,
                ice_density_factor = ice_density_factor
            )

        freeze_thaw_heat_flux =
            freeze_thaw.latent_heat_flux

        freeze_thaw_water_flux =
            freeze_thaw.water_flux

    else
        potential_freeze_thaw_heat =
            0.0

        freeze_thaw_heat_flux =
            0.0

        freeze_thaw_water_flux =
            0.0

        if zero_vapor_exchange_when_freeze_thaw_inactive
            vapor_liquid_flux =
                0.0

            vapor_liquid_latent_heat =
                0.0
        end
    end


    # -------------------------------------------------------------------------
    # 6. Return diagnostics
    # -------------------------------------------------------------------------

    return (
        liquid_water_litter =
            liquid_water_litter,

        vapor_storage_litter =
            vapor_storage_litter,

        ice_storage_litter =
            ice_storage_litter,

        free_pore_volume =
            free_pore_volume,

        equilibrium_vapor_concentration =
            equilibrium_vapor_concentration,

        potential_vapor_liquid_flux =
            potential_vapor_liquid_flux,

        vapor_liquid_flux =
            vapor_liquid_flux,

        vapor_liquid_latent_heat =
            vapor_liquid_latent_heat,

        freezing_temperature =
            freezing_temperature,

        freeze_thaw_active =
            freeze_thaw_active,

        potential_freeze_thaw_heat =
            potential_freeze_thaw_heat,

        freeze_thaw_water_flux =
            freeze_thaw_water_flux,

        freeze_thaw_heat_flux =
            freeze_thaw_heat_flux
    )
end