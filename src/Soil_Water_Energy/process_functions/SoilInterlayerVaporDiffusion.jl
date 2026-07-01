# =============================================================================
# SoilInterlayerVaporDiffusion.jl
#
# Process-level function for vapor diffusion between adjacent soil micropore
# air domains.
#
# Design:
#   - pure function
#   - one source-destination soil interface at a time
#   - no mutation
#   - no @unpack
#   - no direct access to FLVL / HFLWL arrays
# =============================================================================


@doc raw"""
    soil_interlayer_vapor_diffusion_process(...)

Calculate vapor diffusion between adjacent soil micropore air domains.

```julia
# VAPOR PRESSURE AND DIFFUSIVITY IN EACH GRID CELL
```

The process is active only when both source and destination micropore air
volumes exceed a small numerical gate:

```math
V_{a,1} > \epsilon_a
\quad \text{and} \quad
V_{a,2} > \epsilon_a
```

where ``V_a`` is the substep micropore air volume used by the vapor routine.

The vapor concentrations are

```math
C_{v,1} =
\max\left(0, \frac{V_{v,1}}{V_{a,1}}\right)
```

and

```math
C_{v,2} =
\max\left(0, \frac{V_{v,2}}{V_{a,2}}\right)
```

where ``V_v`` is vapor storage and ``V_a`` is air volume.

The effective vapor diffusivities are

```math
D_{v,1} =
D_{g,1} V_{a,1}^{diff} q_p
\frac{V_{a,1}^{diff}}{\phi_1}
```

and

```math
D_{v,2} =
D_{g,2} V_{a,2}^{diff} q_p
\frac{V_{a,2}^{diff}}{\phi_2}
```

where ``D_g`` is the gas diffusivity coefficient, ``V_a^{diff}`` is the
air-volume term used for diffusivity, ``q_p`` is the legacy porosity/tortuosity
coefficient, and ``\phi`` is porosity.

The interfacial vapor conductance is

```math
D_{v,12} =
\frac{
    2 D_{v,1} D_{v,2}
}{
    D_{v,1} L_2 + D_{v,2} L_1
}
```

The potential diffusive vapor flux is

```math
F_v^* =
D_{v,12} A (C_{v,1} - C_{v,2}) \Delta t
```

The mixed equilibrium vapor concentration is

```math
C_{v,y} =
\frac{
    C_{v,1} V_{a,1} + C_{v,2} V_{a,2}
}{
    V_{a,1} + V_{a,2}
}
```

and the storage-limited source-side vapor exchange is

```math
F_{v,x} =
(C_{v,1} - C_{v,y}) V_{a,1} f_r
```

The final vapor flux is

```math
F_v =
\begin{cases}
\max(0, \min(F_v^*, F_{v,x})), & F_v^* \ge 0 \\
\min(0, \max(F_v^*, F_{v,x})), & F_v^* < 0
\end{cases}
```

The advective heat carried by vapor movement is

```math
H_v =
\begin{cases}
c_w T_1 F_v, & F_v > 0 \\
c_w T_2 F_v, & F_v < 0 \\
0, & F_v = 0
\end{cases}
```

implemented by

```julia
advective_heat_by_vapor_flux(F_v, T_1, T_2, cpw)
```

# Sign convention

- `vapor_flux > 0`: vapor moves from source layer to destination layer.
- `vapor_flux < 0`: vapor moves from destination layer to source layer.

# Returns

A named tuple with activation status, vapor concentrations, vapor diffusivities,
interfacial conductance, potential vapor flux, limited vapor flux, and advective
heat flux.
"""
function soil_interlayer_vapor_diffusion_process(;
    # -------------------------------------------------------------------------
    # Activation gate
    # -------------------------------------------------------------------------
    source_air_volume_gate,
    destination_air_volume_gate,
    minimum_air_volume,

    # -------------------------------------------------------------------------
    # Vapor storage and air volumes
    # -------------------------------------------------------------------------
    source_vapor_storage,
    destination_vapor_storage,

    source_air_volume_for_diffusion,
    destination_air_volume_for_diffusion,

    # -------------------------------------------------------------------------
    # Gas diffusivity and porosity
    # -------------------------------------------------------------------------
    source_gas_diffusivity,
    destination_gas_diffusivity,

    source_porosity,
    destination_porosity,

    porosity_tortuosity_factor,

    # -------------------------------------------------------------------------
    # Geometry and time
    # -------------------------------------------------------------------------
    source_flow_length,
    destination_flow_length,

    exchange_area,
    process_time_factor,
    storage_removal_fraction,

    # -------------------------------------------------------------------------
    # Heat
    # -------------------------------------------------------------------------
    source_temperature,
    destination_temperature,
    cpw,

    # -------------------------------------------------------------------------
    # Numerical controls
    # -------------------------------------------------------------------------
    tiny = tiny_num
)
    # -------------------------------------------------------------------------
    # 1. Activation gate
    #
    # Preserve WF8 bug-fix gate:
    #
    #     soilAir_MicporeM[M, N] > tiny_num * Area[3, NU]
    #
    # The already calculated threshold is passed as `minimum_air_volume`.
    # -------------------------------------------------------------------------

    is_active =
        source_air_volume_gate > minimum_air_volume &&
        destination_air_volume_gate > minimum_air_volume

    if !is_active
        return (
            is_active = false,

            source_vapor_concentration = 0.0,
            destination_vapor_concentration = 0.0,

            source_vapor_diffusivity = 0.0,
            destination_vapor_diffusivity = 0.0,

            vapor_interface_conductance = 0.0,

            potential_vapor_flux = 0.0,
            mixed_vapor_concentration = 0.0,
            source_limited_vapor_flux = 0.0,

            vapor_flux = 0.0,
            heat_flux = 0.0
        )
    end


    # -------------------------------------------------------------------------
    # 2. Vapor concentrations
    # -------------------------------------------------------------------------

    source_vapor_concentration =
        max(
            0.0,
            source_vapor_storage / source_air_volume_gate
        )

    destination_vapor_concentration =
        max(
            0.0,
            destination_vapor_storage / destination_air_volume_gate
        )


    # -------------------------------------------------------------------------
    # 3. Cell vapor diffusivities
    #
    # Preserve legacy WF8 form:
    #
    # CNV = soil_gasDiffuty[8, L] *
    #       soilAir_Vol_M[M, L] *
    #       POROQ *
    #       soilAir_Vol_M[M, L] /
    #       POROS[L]
    # -------------------------------------------------------------------------

    source_vapor_diffusivity =
        source_gas_diffusivity *
        source_air_volume_for_diffusion *
        porosity_tortuosity_factor *
        source_air_volume_for_diffusion /
        source_porosity

    destination_vapor_diffusivity =
        destination_gas_diffusivity *
        destination_air_volume_for_diffusion *
        porosity_tortuosity_factor *
        destination_air_volume_for_diffusion /
        destination_porosity


    # -------------------------------------------------------------------------
    # 4. Interfacial vapor conductance
    #
    # Preserve legacy WF8 form:
    #
    # ATCNVL = 2 CNV1 CNVL /
    #          (CNV1 * soil_cube[N, N6] + CNVL * soil_cube[N, N3])
    #
    # Use series_interface_conductance(...) for safety and consistency.
    # -------------------------------------------------------------------------

    vapor_interface_conductance =
        series_interface_conductance(
            source_vapor_diffusivity,
            destination_vapor_diffusivity,
            source_flow_length,
            destination_flow_length;
            tiny = tiny
        )


    # -------------------------------------------------------------------------
    # 5. Potential diffusive vapor flux
    # -------------------------------------------------------------------------

    potential_vapor_flux =
        diffusive_vapor_exchange(
            source_vapor_concentration,
            destination_vapor_concentration,
            vapor_interface_conductance,
            exchange_area,
            process_time_factor
        )


    # -------------------------------------------------------------------------
    # 6. Mixed concentration and source-side storage-limited vapor exchange
    # -------------------------------------------------------------------------

    total_air_volume =
        source_air_volume_gate +
        destination_air_volume_gate

    mixed_vapor_concentration =
        if total_air_volume > tiny
            (
                source_vapor_concentration * source_air_volume_gate +
                destination_vapor_concentration * destination_air_volume_gate
            ) /
            total_air_volume
        else
            0.0
        end

    source_limited_vapor_flux =
        (
            source_vapor_concentration -
            mixed_vapor_concentration
        ) *
        source_air_volume_gate *
        storage_removal_fraction


    # -------------------------------------------------------------------------
    # 7. Limit vapor flux
    #
    # Preserve WF8:
    #
    # if FLVC >= 0
    #     FLVQ = max(0, min(FLVC, FLVX))
    # else
    #     FLVQ = min(0, max(FLVC, FLVX))
    # end
    # -------------------------------------------------------------------------

    vapor_flux =
        if potential_vapor_flux >= 0.0
            max(
                0.0,
                min(
                    potential_vapor_flux,
                    source_limited_vapor_flux
                )
            )
        else
            min(
                0.0,
                max(
                    potential_vapor_flux,
                    source_limited_vapor_flux
                )
            )
        end


    # -------------------------------------------------------------------------
    # 8. Advective heat carried by vapor movement
    # -------------------------------------------------------------------------

    heat_flux =
        advective_heat_by_vapor_flux(
            vapor_flux,
            source_temperature,
            destination_temperature,
            cpw
        )


    # -------------------------------------------------------------------------
    # 9. Return diagnostics
    # -------------------------------------------------------------------------

    return (
        is_active = true,

        source_vapor_concentration =
            source_vapor_concentration,

        destination_vapor_concentration =
            destination_vapor_concentration,

        source_vapor_diffusivity =
            source_vapor_diffusivity,

        destination_vapor_diffusivity =
            destination_vapor_diffusivity,

        vapor_interface_conductance =
            vapor_interface_conductance,

        potential_vapor_flux =
            potential_vapor_flux,

        mixed_vapor_concentration =
            mixed_vapor_concentration,

        source_limited_vapor_flux =
            source_limited_vapor_flux,

        vapor_flux =
            vapor_flux,

        heat_flux =
            heat_flux
    )
end