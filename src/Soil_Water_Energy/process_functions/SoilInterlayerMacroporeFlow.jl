# =============================================================================
# SoilInterlayerMacroporeFlow.jl
#
# Process-level function for liquid-water exchange between adjacent soil
# macropore domains.
#
# Design:
#   - pure function
#   - one source-destination soil interface at a time
#   - no mutation
#   - no @unpack
#   - no direct access to FLWHL / FLWHM / HFLWL arrays
# =============================================================================


@doc raw"""
    soil_interlayer_macropore_flow_process(...)

Calculate liquid-water flow between adjacent soil macropore domains.

```julia
# MACROPORE FLOW FROM POISEUILLE FLOW if MACROPORES PRESENT
```

Macropore water potential is approximated as gravitational potential plus
a storage-depth correction:

```math
\Psi_M =
\psi_g +
g_h L_z
\left[
    \min\left(1, \max\left(0, \frac{W_M}{V_M}\right)\right)
    - 0.5
\right]
```

where ``\psi_g`` is gravitational potential, ``g_h`` is the hydrostatic
gradient, ``L_z`` is vertical layer thickness, ``W_M`` is macropore liquid
water, and ``V_M`` is macropore volume.

The potential macropore flux is

```math
F_M^* =
C_M A (\Psi_{M,1} - \Psi_{M,2}) \Delta t
```

where ``C_M`` is macropore conductance, ``A`` is exchange area, and
``\Delta t`` is the process time factor.

For lateral interfaces, flux can move in either direction:

```math
F_M =
\begin{cases}
\min(F_M^*, W_{M,1} f_r, P_{M,2} f_r), & \Psi_{M,1} > \Psi_{M,2} \\
\max(F_M^*, -W_{M,2} f_r, -P_{M,1} f_r), & \Psi_{M,1} < \Psi_{M,2} \\
0, & \Psi_{M,1} = \Psi_{M,2}
\end{cases}
```

For vertical interfaces, legacy WF8 allows only positive/downward macropore
flow:

```math
F_M =
\max\left[
    0,
    \min\left(
        W_{M,1} f_r + F_{M,in},
        P_{M,2} f_r,
        F_M^*
    \right)
\right]
```

where ``F_{M,in}`` is the already accumulated vertical macropore inflow into
the source layer.

During freezing, vertical transfer is additionally corrected by excess
macropore water+ice in the destination layer:

```math
F_M \leftarrow F_M + \min(0, E_{M,2})
```

The advective heat flux is calculated from the temperature of the compartment
where the moving water originates:

```math
H_M =
\begin{cases}
c_w T_1 F_M, & F_M > 0 \\
c_w T_2 F_M, & F_M < 0 \\
0, & F_M = 0
\end{cases}
```

implemented by

```julia
advective_heat_by_water_flux(F_M, T_1, T_2, cpw)
```

# Sign convention

- `limited_flux > 0`: macropore water moves from source to destination.
- `limited_flux < 0`: macropore water moves from destination to source.

# Returns

A named tuple with activation status, macropore potentials, potential flux,
limited flux, heat flux, and updated macropore block flag.
"""
function soil_interlayer_macropore_flow_process(;
    # -------------------------------------------------------------------------
    # Activation state
    # -------------------------------------------------------------------------
    source_macropore_volume,
    destination_macropore_volume,
    macropore_block_flag,

    # -------------------------------------------------------------------------
    # Macropore water and air-capacity state
    # -------------------------------------------------------------------------
    source_macropore_water,
    destination_macropore_water,
    source_macropore_air_capacity,
    destination_macropore_air_capacity,
    destination_excess_macropore_water_ice,

    # -------------------------------------------------------------------------
    # Potentials and geometry
    # -------------------------------------------------------------------------
    source_gravity_potential,
    destination_gravity_potential,
    source_vertical_length,
    destination_vertical_length,

    macropore_conductance,
    exchange_area,

    # -------------------------------------------------------------------------
    # Direction / topology
    # -------------------------------------------------------------------------
    is_vertical_boundary,
    source_incoming_macropore_flux,

    # -------------------------------------------------------------------------
    # Time and storage factors
    # -------------------------------------------------------------------------
    process_time_factor,
    storage_removal_fraction,

    # -------------------------------------------------------------------------
    # Heat
    # -------------------------------------------------------------------------
    source_temperature,
    destination_temperature,
    cpw,

    # -------------------------------------------------------------------------
    # Numerical constants
    # -------------------------------------------------------------------------
    hydraulic_gradient = 0.0098,
    tiny = tiny_num2
)
    # -------------------------------------------------------------------------
    # 1. Activation gate
    # -------------------------------------------------------------------------

    is_active =
        source_macropore_volume > tiny &&
        destination_macropore_volume > tiny &&
        macropore_block_flag == 0

    if !is_active
        updated_macropore_block_flag =
            destination_macropore_air_capacity <= 0.0 ? 1 : macropore_block_flag

        return (
            is_active = false,

            source_macropore_potential = 0.0,
            destination_macropore_potential = 0.0,

            source_macropore_fraction = 0.0,
            destination_macropore_fraction = 0.0,

            potential_flux = 0.0,
            limited_flux = 0.0,
            heat_flux = 0.0,

            excess_macropore_correction = 0.0,

            updated_macropore_block_flag = updated_macropore_block_flag
        )
    end


    # -------------------------------------------------------------------------
    # 2. Macropore water fractions
    # -------------------------------------------------------------------------

    source_macropore_fraction =
        min(
            1.0,
            max(
                0.0,
                source_macropore_water / source_macropore_volume
            )
        )

    destination_macropore_fraction =
        min(
            1.0,
            max(
                0.0,
                destination_macropore_water / destination_macropore_volume
            )
        )


    # -------------------------------------------------------------------------
    # 3. Macropore water potentials
    # -------------------------------------------------------------------------

    source_macropore_potential =
        source_gravity_potential +
        hydraulic_gradient *
        source_vertical_length *
        (source_macropore_fraction - 0.5)

    destination_macropore_potential =
        destination_gravity_potential +
        hydraulic_gradient *
        destination_vertical_length *
        (destination_macropore_fraction - 0.5)


    # -------------------------------------------------------------------------
    # 4. Potential macropore water flux
    # -------------------------------------------------------------------------

    potential_flux =
        water_flux_from_potential(
            source_macropore_potential,
            destination_macropore_potential,
            macropore_conductance,
            exchange_area,
            process_time_factor
        )


    # -------------------------------------------------------------------------
    # 5. Storage-limited macropore flux
    # -------------------------------------------------------------------------

    if !is_vertical_boundary
        if source_macropore_potential > destination_macropore_potential
            limited_flux =
                min(
                    potential_flux,
                    source_macropore_water * storage_removal_fraction,
                    destination_macropore_air_capacity * storage_removal_fraction
                )

        elseif source_macropore_potential < destination_macropore_potential
            limited_flux =
                max(
                    potential_flux,
                    -destination_macropore_water * storage_removal_fraction,
                    -source_macropore_air_capacity * storage_removal_fraction
                )

        else
            limited_flux =
                0.0
        end

    else
        limited_flux =
            max(
                0.0,
                min(
                    min(
                        source_macropore_water * storage_removal_fraction +
                        source_incoming_macropore_flux,

                        destination_macropore_air_capacity *
                        storage_removal_fraction
                    ),
                    potential_flux
                )
            )
    end


    # -------------------------------------------------------------------------
    # 6. Excess macropore water+ice correction during freezing
    # -------------------------------------------------------------------------

    excess_macropore_correction =
        if is_vertical_boundary
            min(
                0.0,
                destination_excess_macropore_water_ice
            )
        else
            0.0
        end

    limited_flux +=
        excess_macropore_correction


    # -------------------------------------------------------------------------
    # 7. Advective heat carried by macropore water
    # -------------------------------------------------------------------------

    heat_flux =
        advective_heat_by_water_flux(
            limited_flux,
            source_temperature,
            destination_temperature,
            cpw
        )


    # -------------------------------------------------------------------------
    # 8. Return diagnostics
    # -------------------------------------------------------------------------

    return (
        is_active = true,

        source_macropore_potential =
            source_macropore_potential,

        destination_macropore_potential =
            destination_macropore_potential,

        source_macropore_fraction =
            source_macropore_fraction,

        destination_macropore_fraction =
            destination_macropore_fraction,

        potential_flux =
            potential_flux,

        limited_flux =
            limited_flux,

        heat_flux =
            heat_flux,

        excess_macropore_correction =
            excess_macropore_correction,

        updated_macropore_block_flag =
            macropore_block_flag
    )
end
