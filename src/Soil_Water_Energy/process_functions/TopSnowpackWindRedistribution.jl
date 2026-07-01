# =============================================================================
# TopSnowpackWindRedistribution.jl
#
# Process-level function for wind redistribution/removal of material from the
# top snowpack layer.
#
# Design:
#   - pure function
#   - one surface/top-snowpack state at a time
#   - no mutation
#   - no @unpack
#   - no direct access to QSM / QWM / QIM / QST arrays
# =============================================================================


@doc raw"""
    top_snowpack_wind_redistribution_process(...)

Calculate wind redistribution/removal of snow, liquid water, and ice from the
top snowpack layer.

The process is active only when snow depth is positive:

```math
D_s > 0
```

where ``D_s`` is snow depth.

The redistribution fraction follows the legacy WF6 expression:

```math
f_d =
10^{-7} U_a \Delta t
```

where ``U_a`` is wind speed and ``\Delta t`` is the process time factor.

The redistributed snow, liquid water, and ice fluxes are

```math
Q_s = f_d S
```

```math
Q_w = f_d W
```

```math
Q_i = f_d I
```

where ``S`` is top-layer snow volume, ``W`` is top-layer liquid-water volume,
and ``I`` is top-layer ice volume.

The total redistributed snowpack material is

```math
Q_t = Q_s + Q_w + Q_i
```

# Sign convention

All returned fluxes are nonnegative redistribution/removal fluxes from the
top snowpack layer.

# Returns

A named tuple with activation status, redistribution fraction, snow flux,
liquid-water flux, ice flux, and total flux.
"""
function top_snowpack_wind_redistribution_process(;
    # -------------------------------------------------------------------------
    # Activation state
    # -------------------------------------------------------------------------
    snow_depth,

    # -------------------------------------------------------------------------
    # Wind and time
    # -------------------------------------------------------------------------
    wind_speed,
    time_factor,

    # -------------------------------------------------------------------------
    # Top snowpack storage
    # -------------------------------------------------------------------------
    top_snow_volume,
    top_liquid_water,
    top_ice_volume,

    # -------------------------------------------------------------------------
    # Numerical / legacy constants
    # -------------------------------------------------------------------------
    snow_depth_threshold = 0.0,
    redistribution_coefficient = 1.0e-7
)
    # -------------------------------------------------------------------------
    # 1. Activation gate
    # -------------------------------------------------------------------------

    is_active =
        snow_depth > snow_depth_threshold

    if !is_active
        return (
            is_active =
                false,

            redistribution_fraction =
                0.0,

            snow_flux =
                0.0,

            water_flux =
                0.0,

            ice_flux =
                0.0,

            total_flux =
                0.0
        )
    end


    # -------------------------------------------------------------------------
    # 2. Legacy wind redistribution fraction
    #
    # Preserve WF6:
    #
    #     QSX = 1.0e-7 * UA * XNPHX
    #
    # -------------------------------------------------------------------------

    redistribution_fraction =
        redistribution_coefficient *
        wind_speed *
        time_factor


    # -------------------------------------------------------------------------
    # 3. Redistributed top-snowpack material
    #
    # Use the pure helper for consistency with BoundaryFluxes.jl.
    # -------------------------------------------------------------------------

    snow_drift =
        wind_redistributed_top_snowpack(
            wind_speed,
            time_factor,
            top_snow_volume,
            top_liquid_water,
            top_ice_volume
        )


    # -------------------------------------------------------------------------
    # 4. Return diagnostics
    # -------------------------------------------------------------------------

    return (
        is_active =
            true,

        redistribution_fraction =
            redistribution_fraction,

        snow_flux =
            snow_drift.snow_flux,

        water_flux =
            snow_drift.water_flux,

        ice_flux =
            snow_drift.ice_flux,

        total_flux =
            snow_drift.total_flux
    )
end