# =============================================================================
# LitterSoilCapillaryExchange.jl
#
# Process-level function for capillary water exchange between surface litter
# and the upper soil micropore domain.
#
# Design:
#   - pure function
#   - one litter-soil exchange substep at a time
#   - no mutation
#   - no @unpack
#   - no direct access to waterVar_copy
# =============================================================================


@doc raw"""
    soil_litter_capillary_exchange_process(...)

Calculate one substep of capillary water exchange between surface litter and
the upper soil layer.

This process represents water movement across the litter--soil interface
driven by total water-potential differences.

The bounded volumetric water contents are

```math
\theta_l =
\min\left[
    \theta_{l,\max},
    \max\left(0, \frac{W_l}{V_l}\right)
\right]
````

and

```math
\theta_s =
\min\left[
    \phi_s,
    \max\left(\theta_{s,\min}, \frac{W_s}{V_s}\right)
\right]
```

where `W_l` and `W_s` are litter and soil water storages,
`V_l` is dry litter volume, `V_s` is soil reference volume,
`\theta_{l,\max}` is the litter saturation water content, and
`\phi_s` is soil porosity.

The total water potentials are

```math
\Psi_l = \psi_{m,l} + \psi_{g,l} + \psi_{o,l}
```

and

```math
\Psi_s = \psi_{m,s} + \psi_{g,s} + \psi_{o,s}
```

The litter--soil interface conductance is

```math
C_{ls} =
\frac{2 K_l K_s}{K_l L_s + K_s L_l}
```

where `K_l` and `K_s` are litter and soil hydraulic conductivities,
and `L_l` and `L_s` are effective flow lengths.

The potential capillary flux is

```math
F^* =
C_{ls} A f_w (\Psi_l - \Psi_s) \Delta t
```

where `A` is exchange area, `f_w` is the water-covered fraction, and
`\Delta t` is the capillary exchange time factor.

When the potential flux is near zero, a deadband is applied:

```math
F^*_0 =
\begin{cases}
0, & |F^*| < \epsilon \\
F^*, & |F^*| \ge \epsilon
\end{cases}
```

The air-entry correction is

```math
F_z =
\begin{cases}
F^*_0 +
\min\left[
    (\theta_l - \theta_{ae,l}) V_l,
    \max(0, (\theta_{ae,s} - \theta_s) V_s)
\right] f_r,
& F^*_0 > 0,\ \theta_l > \theta_{ae,l}
\\
F^*_0 +
\max\left[
    (\theta_{ae,s} - \theta_s) V_s,
    \min(0, (\theta_l - \theta_{ae,l}) V_l)
\right] f_r,
& F^*_0 < 0,\ \theta_s > \theta_{ae,s}
\\
F^*_0,
& \text{otherwise}
\end{cases}
```

The final flux is storage-limited:

```math
F =
\begin{cases}
\min(F_z, W_l f_r, P_s), & F_z > 0 \\
\max(F_z, -W_s f_r, -P_l), & F_z < 0
\end{cases}
```

where positive flux means litter-to-soil flow, `P_l` is litter pore
capacity, and `P_s` is soil pore capacity.

If the soil micropore domain has negative capacity, an additional correction
is applied:

```math
F \leftarrow F +
\min\left[
    0,
    \max(-W_s f_r, P_{s,raw})
\right]
\quad \text{if } P_{s,raw} < 0
```

The advective heat flux is

```math
H_F =
\begin{cases}
c_w F T_l, & F > 0 \\
c_w F T_s, & F < 0
\end{cases}
```

implemented by

```julia
advective_heat_by_water_flux(F, T_l, T_s, cpw)
```

# Sign convention

* `limited_flux > 0`: water moves from surface litter to soil.
* `limited_flux < 0`: water moves from soil to surface litter.

# Returns

A named tuple containing water contents, potentials, conductance, potential
flux, corrected flux, limited flux, heat flux, and updated substep storages.
"""
function soil_litter_capillary_exchange_process(;
# -------------------------------------------------------------------------
# Current substep storages
# -------------------------------------------------------------------------
litter_water,
soil_water,

litter_ice,
soil_ice,

litter_water_capacity,
soil_pore_volume,
soil_air_capacity,
soil_raw_air_capacity,

# -------------------------------------------------------------------------
# Litter and soil physical properties
# -------------------------------------------------------------------------
litter_dry_volume,
soil_reference_volume,
litter_saturated_porosity,
soil_porosity,

litter_field_capacity,
soil_field_capacity,
field_capacity_potential,

litter_wilting_point,
soil_wilting_point,
wilting_point_potential,

litter_shape_parameter,
soil_shape_parameter,

litter_air_entry_water_content,
soil_air_entry_water_content,

# -------------------------------------------------------------------------
# Potentials and conductivities
# -------------------------------------------------------------------------
litter_gravity_potential,
soil_gravity_potential,
litter_osmotic_potential,
soil_osmotic_potential,

litter_hydraulic_conductivity_table,
soil_hydraulic_conductivity_table,
# CLAUDE DEBUG BEGIN: soil_HydroConductivity is 3-D (direction, k_index, layer).
# Caller must pass the soil-layer index so the lookup is well-formed.
soil_layer_index,
# CLAUDE DEBUG END
hydraulic_conductivity_multiplier,

# -------------------------------------------------------------------------
# Geometry and time
# -------------------------------------------------------------------------
litter_flow_length,
soil_flow_length,
exchange_area,
water_cover_fraction,
capillary_time_factor,
storage_removal_fraction,

# -------------------------------------------------------------------------
# Heat
# -------------------------------------------------------------------------
litter_temperature,
soil_temperature,
cpw,

# -------------------------------------------------------------------------
# Numerical controls
# -------------------------------------------------------------------------
dry_volume_fallback_water_content,
soil_minimum_water_content,
cap_deadband = 1.0e-6,
flux_deadband = tiny_num,
tiny = tiny_num2

)
# -------------------------------------------------------------------------
# 1. Litter volumetric water content
#
# Preserve the current WF5 cap deadband:
# if litter water is within cap_deadband of holding capacity, clamp
# to holding capacity before calculating θ_l.
# -------------------------------------------------------------------------

if litter_dry_volume > tiny
    litter_water_for_theta =
        litter_water >= litter_water_capacity - cap_deadband ?
        litter_water_capacity :
        litter_water

    litter_water_content =
        bounded_volumetric_water_content(
            litter_water_for_theta,
            litter_dry_volume,
            0.0,
            litter_water_capacity / litter_dry_volume
        )
else
    litter_water_content =
        dry_volume_fallback_water_content
end


# -------------------------------------------------------------------------
# 2. Litter matric potential
# -------------------------------------------------------------------------

litter_matric_potential, litter_water_content =
    cal_matricWaterPotential_MCM(
        litter_water_content,
        litter_field_capacity,
        field_capacity_potential,
        litter_wilting_point,
        wilting_point_potential,
        litter_saturated_porosity,
        litter_shape_parameter
    )


# -------------------------------------------------------------------------
# 3. Soil volumetric water content and matric potential
# -------------------------------------------------------------------------

soil_water_content =
    bounded_volumetric_water_content(
        soil_water,
        soil_reference_volume,
        soil_minimum_water_content,
        soil_porosity
    )

soil_matric_potential, _ =
    cal_matricWaterPotential_MCM(
        soil_water_content,
        soil_field_capacity,
        field_capacity_potential,
        soil_wilting_point,
        wilting_point_potential,
        soil_porosity,
        soil_shape_parameter
    )


# -------------------------------------------------------------------------
# 4. Hydraulic conductivity lookup
# -------------------------------------------------------------------------

litter_conductivity_index =
    conductivity_table_index(
        litter_water_content,
        litter_saturated_porosity
    )

soil_conductivity_index =
    conductivity_table_index(
        soil_water_content,
        soil_porosity
    )

litter_hydraulic_conductivity =
    litter_hydraulic_conductivity_table[3, litter_conductivity_index]

# CLAUDE DEBUG BEGIN: index soil_HydroConductivity with layer dim.
soil_hydraulic_conductivity =
    soil_hydraulic_conductivity_table[3, soil_conductivity_index, soil_layer_index] *
    hydraulic_conductivity_multiplier
# CLAUDE DEBUG END


# -------------------------------------------------------------------------
# 5. Interface conductance
# -------------------------------------------------------------------------

interface_conductance_value =
    interface_conductance(
        litter_hydraulic_conductivity,
        soil_hydraulic_conductivity,
        litter_flow_length,
        soil_flow_length
    )


# -------------------------------------------------------------------------
# 6. Total water potentials
# -------------------------------------------------------------------------

litter_total_potential =
    total_water_potential(
        litter_matric_potential,
        litter_gravity_potential,
        litter_osmotic_potential
    )

soil_total_potential =
    total_water_potential(
        soil_matric_potential,
        soil_gravity_potential,
        soil_osmotic_potential
    )


# -------------------------------------------------------------------------
# 7. Potential capillary flux
# -------------------------------------------------------------------------

potential_flux =
    water_flux_from_potential(
        litter_total_potential,
        soil_total_potential,
        interface_conductance_value,
        exchange_area * water_cover_fraction,
        capillary_time_factor
    )


# -------------------------------------------------------------------------
# 8. Air-entry correction
# -------------------------------------------------------------------------

if abs(potential_flux) < flux_deadband
    corrected_potential_flux =
        0.0

elseif potential_flux >= 0.0
    if litter_water_content > litter_air_entry_water_content
        corrected_potential_flux =
            potential_flux +
            min(
                (litter_water_content - litter_air_entry_water_content) *
                litter_dry_volume,

                max(
                    0.0,
                    (soil_air_entry_water_content - soil_water_content) *
                    soil_reference_volume
                )
            ) *
            storage_removal_fraction
    else
        corrected_potential_flux =
            potential_flux
    end

else
    if soil_water_content > soil_air_entry_water_content
        corrected_potential_flux =
            potential_flux +
            max(
                (soil_air_entry_water_content - soil_water_content) *
                soil_reference_volume,

                min(
                    0.0,
                    (litter_water_content - litter_air_entry_water_content) *
                    litter_dry_volume
                )
            ) *
            storage_removal_fraction
    else
        corrected_potential_flux =
            potential_flux
    end
end


# -------------------------------------------------------------------------
# 9. Storage-limited exchange flux
# -------------------------------------------------------------------------

litter_available_water =
    max(0.0, litter_water)

soil_available_water =
    max(0.0, soil_water)

litter_pore_capacity =
    max(
        0.0,
        litter_water_capacity - litter_water
    )

limited_flux =
    limit_bidirectional_flux(
        corrected_potential_flux,
        litter_available_water * storage_removal_fraction,
        soil_air_capacity,
        soil_available_water * storage_removal_fraction,
        litter_pore_capacity
    )


# -------------------------------------------------------------------------
# 10. Correct for negative raw soil air capacity
# -------------------------------------------------------------------------

if soil_raw_air_capacity < 0.0
    limited_flux +=
        min(
            0.0,
            max(
                -soil_available_water * storage_removal_fraction,
                soil_raw_air_capacity
            )
        )
end


# -------------------------------------------------------------------------
# 11. Advective heat carried by capillary water exchange
# -------------------------------------------------------------------------

heat_flux =
    advective_heat_by_water_flux(
        limited_flux,
        litter_temperature,
        soil_temperature,
        cpw
    )


# -------------------------------------------------------------------------
# 12. Updated substep storages
# -------------------------------------------------------------------------

updated_litter_water =
    litter_water - limited_flux

updated_soil_water =
    soil_water + limited_flux

updated_litter_pore_capacity =
    max(
        0.0,
        litter_water_capacity - updated_litter_water
    )

updated_soil_raw_air_capacity =
    soil_pore_volume -
    updated_soil_water -
    soil_ice

updated_soil_air_capacity =
    max(
        0.0,
        updated_soil_raw_air_capacity
    )


# -------------------------------------------------------------------------
# 13. Return process diagnostics
# -------------------------------------------------------------------------

return (
    litter_water_content =
        litter_water_content,

    soil_water_content =
        soil_water_content,

    litter_matric_potential =
        litter_matric_potential,

    soil_matric_potential =
        soil_matric_potential,

    litter_total_potential =
        litter_total_potential,

    soil_total_potential =
        soil_total_potential,

    litter_conductivity_index =
        litter_conductivity_index,

    soil_conductivity_index =
        soil_conductivity_index,

    litter_hydraulic_conductivity =
        litter_hydraulic_conductivity,

    soil_hydraulic_conductivity =
        soil_hydraulic_conductivity,

    interface_conductance =
        interface_conductance_value,

    potential_flux =
        potential_flux,

    corrected_potential_flux =
        corrected_potential_flux,

    limited_flux =
        limited_flux,

    heat_flux =
        heat_flux,

    updated_litter_water =
        updated_litter_water,

    updated_soil_water =
        updated_soil_water,

    updated_litter_pore_capacity =
        updated_litter_pore_capacity,

    updated_soil_raw_air_capacity =
        updated_soil_raw_air_capacity,

    updated_soil_air_capacity =
        updated_soil_air_capacity
)


end
