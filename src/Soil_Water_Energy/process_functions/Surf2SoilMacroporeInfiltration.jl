# =============================================================================
# SurfaceToSoilMacroporeInfiltration.jl
#
# Process-level function for infiltration of free surface water into the
# top-soil macropore domain.
#
# Design:
#   - pure function
#   - one surface/top-soil state at a time
#   - no mutation
#   - no @unpack
#   - no direct access to FLWHL / HFLWL / FLWRL / HFLWRL arrays
# =============================================================================


@doc raw"""
    surface_to_soil_macropore_infiltration_process(...)

Calculate infiltration of free surface water into the top-soil macropore
domain.

The process is active when both top-soil macropore air capacity and free
surface water are available:

```math
P_M > \epsilon
\quad \text{and} \quad
W_f > \epsilon
```

where ``P_M`` is top-soil macropore air capacity and ``W_f`` is free surface
water.

The potential infiltrating water is limited by removable free surface water
and available macropore air capacity:

```math
F_M =
\min(W_f f_r, P_M)
```

where ``f_r`` is the substep removal fraction.

The associated advective heat flux is

```math
H_M =
c_w T_s F_M
```

implemented by

```julia
advective_heat_water(F_M, surface_temperature, cpw)
```

# Sign convention

- `infiltration_flux > 0`: free surface water enters the top-soil macropore
  domain.
- This process does not allow reverse flow.

# Returns

A named tuple with activation status, infiltration water flux, heat flux,
remaining free surface water, and remaining macropore air capacity.
"""
function surface_to_soil_macropore_infiltration_process(;
    # -------------------------------------------------------------------------
    # Surface/top-soil state
    # -------------------------------------------------------------------------
    free_surface_water,
    topsoil_macropore_air_capacity,

    # -------------------------------------------------------------------------
    # Time and heat
    # -------------------------------------------------------------------------
    storage_removal_fraction,
    surface_temperature,
    cpw,

    # -------------------------------------------------------------------------
    # Numerical controls
    # -------------------------------------------------------------------------
    tiny = tiny_num2
)
    # -------------------------------------------------------------------------
    # 1. Activation gate
    # -------------------------------------------------------------------------

    is_active =
        topsoil_macropore_air_capacity > tiny &&
        free_surface_water > tiny

    if !is_active
        return (
            is_active =
                false,

            infiltration_flux =
                0.0,

            heat_flux =
                0.0,

            remaining_free_surface_water =
                free_surface_water,

            remaining_topsoil_macropore_air_capacity =
                topsoil_macropore_air_capacity
        )
    end


    # -------------------------------------------------------------------------
    # 2. Storage-limited infiltration into top-soil macropores
    #
    # Preserve WF6:
    #
    #     FLQHR = min(water_freeSurfX * XNPXX, soilAir_Macpore_L[NUM])
    #
    # -------------------------------------------------------------------------

    infiltration_flux =
        min(
            free_surface_water * storage_removal_fraction,
            topsoil_macropore_air_capacity
        )


    # -------------------------------------------------------------------------
    # 3. Advective heat carried by infiltrating surface water
    # -------------------------------------------------------------------------

    heat_flux =
        advective_heat_water(
            infiltration_flux,
            surface_temperature,
            cpw
        )


    # -------------------------------------------------------------------------
    # 4. Return updated local storage diagnostics
    # -------------------------------------------------------------------------

    return (
        is_active =
            true,

        infiltration_flux =
            infiltration_flux,

        heat_flux =
            heat_flux,

        remaining_free_surface_water =
            free_surface_water - infiltration_flux,

        remaining_topsoil_macropore_air_capacity =
            topsoil_macropore_air_capacity - infiltration_flux
    )
end