# =============================================================================
# SoilInterlayerMicroporeFlow.jl
#
# Process-level function for liquid-water exchange between adjacent soil
# micropore domains.
#
# Design:
#   - pure function
#   - one source-destination soil interface at a time
#   - no mutation
#   - no @unpack
#   - no direct access to FLWL / FLWLX / HFLWL / soilW_cp2 arrays
# =============================================================================


@doc raw"""
    soil_interlayer_micropore_flow_process(...)

Calculate liquid-water flow between two adjacent soil micropore domains.
It includes:
1. bounded source and destination water contents,
2. current matric potentials,
3. Darcy / Green--Ampt / Richards branch selection,
4. hydraulic-conductivity table lookup,
5. total water-potential gradient,
6. interface conductance,
7. potential liquid-water flux,
8. air-entry redistribution correction,
9. storage-limited water flux,
10. excess water+ice correction during freezing,
11. advective heat carried by liquid-water flow,
12. updated temporary available water storages.

The current bounded water contents are

```math
\theta_1 =
\min\left[
    \phi_1,
    \max\left(\theta_{min,1}, \frac{W_1}{V_1}\right)
\right]
````

and

```math
\theta_2 =
\min\left[
    \phi_2,
    \max\left(\theta_{min,2}, \frac{W_2}{V_2}\right)
\right]
```

where `W` is liquid-water storage, `V` is soil reference volume,
`\phi` is porosity, and `\theta_{min}` is the minimum water content
used for hydraulic lookup.

The total water potentials used for liquid flow are

```math
\Psi_1 =
\psi_{m,1} + \psi_{g,1} + \psi_{o,1}
```

and

```math
\Psi_2 =
\psi_{m,2} + \psi_{g,2} + \psi_{o,2}
```

The interface conductance is

```math
C_{12} =
\frac{2 K_1 K_2}{K_1 L_2 + K_2 L_1}
```

The potential flux is

```math
F^* =
C_{12} A (\Psi_1 - \Psi_2) \Delta t
```

A near-zero deadband can optionally preserve legacy symmetric behavior near
saturated interfaces:

```math
F = F^* \quad \text{if } |F^*| < \epsilon
```

Otherwise, the potential flux is corrected by air-entry redistribution and
then limited by available source water and destination air capacity.

# Sign convention

* `limited_flux > 0`: water moves from source layer to destination layer.
* `limited_flux < 0`: water moves from destination layer to source layer.

# Returns

A named tuple containing water contents, matric potentials, vapor potentials,
conductivities, fluxes, advective heat, and updated temporary available water.
"""
function soil_interlayer_micropore_flow_process(;
# -------------------------------------------------------------------------
# Source and destination identity
# -------------------------------------------------------------------------
source_layer_index,
destination_layer_index,
boundary_axis_index,
is_vertical_boundary,
is_surface_soil_layer,


# -------------------------------------------------------------------------
# Current water states used for water-content and matric-potential calculation
# -------------------------------------------------------------------------
source_liquid_water,
destination_liquid_water,

source_wet_front_water,
destination_wet_front_water,

source_available_water,
destination_available_water,

source_air_capacity,
destination_air_capacity,

destination_excess_water_ice,

# -------------------------------------------------------------------------
# Soil volumes and pond/ice-state terms
# -------------------------------------------------------------------------
source_reference_volume,
destination_reference_volume,

source_existing_volume,
destination_existing_volume,

source_liquid_water_volume,
destination_liquid_water_volume,

source_ice_volume,
destination_ice_volume,

# -------------------------------------------------------------------------
# Soil physical properties
# -------------------------------------------------------------------------
source_soil_mass,
destination_soil_mass,

source_porosity,
destination_porosity,

source_minimum_water_content,
destination_minimum_water_content,

source_field_capacity,
destination_field_capacity,

source_wilting_point,
destination_wilting_point,

source_shape_parameter,
destination_shape_parameter,

source_air_entry_water_content,
destination_air_entry_water_content,

source_air_entry_potential,
destination_air_entry_potential,

source_saturation_potential,
destination_saturation_potential,

# -------------------------------------------------------------------------
# Potentials
# -------------------------------------------------------------------------
field_capacity_potential,
wilting_point_potential,

source_gravity_potential,
destination_gravity_potential,

source_osmotic_potential,
destination_osmotic_potential,

# -------------------------------------------------------------------------
# Hydraulic conductivity
# -------------------------------------------------------------------------
hydraulic_conductivity_table,
source_hydraulic_conductivity_multiplier,

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
# Constants / numerical controls
# -------------------------------------------------------------------------
ice_field_capacity_coefficient = FCI,
ice_wilting_point_coefficient = WPI,
potential_flux_deadband = 1.0e-9,
tiny = tiny_num,
tiny_volume = tiny_num2

)
# -------------------------------------------------------------------------
# 1. Current bounded water contents
# -------------------------------------------------------------------------


source_water_content_current =
    bounded_volumetric_water_content(
        source_liquid_water,
        source_reference_volume,
        source_minimum_water_content,
        source_porosity
    )

destination_water_content_current =
    bounded_volumetric_water_content(
        destination_liquid_water,
        destination_reference_volume,
        destination_minimum_water_content,
        destination_porosity
    )


# -------------------------------------------------------------------------
# 2. Current source matric potential
#
#   - soil cell: matric potential from θ
#   - subsurface pond/ice cell: matric potential from liquid-water volume
#   - otherwise: saturation potential
# -------------------------------------------------------------------------

source_current_matric_potential =
    if source_soil_mass > tiny
        first(
            cal_matricWaterPotential_MCM(
                source_water_content_current,
                source_field_capacity,
                field_capacity_potential,
                source_wilting_point,
                wilting_point_potential,
                source_porosity,
                source_shape_parameter
            )
        )
    elseif source_existing_volume > tiny_volume &&
           source_ice_volume > tiny_volume

        source_ice_field_capacity =
            ice_field_capacity_coefficient * source_ice_volume

        source_ice_wilting_point =
            ice_wilting_point_coefficient * source_ice_volume

        first(
            cal_matricWaterPotential_MCM(
                source_liquid_water_volume,
                source_ice_field_capacity,
                field_capacity_potential,
                source_ice_wilting_point,
                wilting_point_potential,
                source_porosity,
                1.0
            )
        )
    else
        source_saturation_potential
    end


# -------------------------------------------------------------------------
# 3. Current destination matric potential
# -------------------------------------------------------------------------

destination_current_matric_potential =
    if destination_soil_mass > tiny
        first(
            cal_matricWaterPotential_MCM(
                destination_water_content_current,
                destination_field_capacity,
                field_capacity_potential,
                destination_wilting_point,
                wilting_point_potential,
                destination_porosity,
                destination_shape_parameter
            )
        )
    elseif destination_existing_volume > tiny_volume &&
           destination_ice_volume > tiny_volume

        destination_ice_field_capacity =
            ice_field_capacity_coefficient * destination_ice_volume

        destination_ice_wilting_point =
            ice_wilting_point_coefficient * destination_ice_volume

        first(
            cal_matricWaterPotential_MCM(
                destination_liquid_water_volume,
                destination_ice_field_capacity,
                field_capacity_potential,
                destination_ice_wilting_point,
                wilting_point_potential,
                destination_porosity,
                1.0
            )
        )
    else
        destination_saturation_potential
    end


# -------------------------------------------------------------------------
# 4. Darcy / Green-Ampt / Richards branch selection
# -------------------------------------------------------------------------

source_is_above_air_entry =
    source_current_matric_potential > source_air_entry_potential

destination_is_above_air_entry =
    destination_current_matric_potential > destination_air_entry_potential

if source_is_above_air_entry && destination_is_above_air_entry
    # Darcy flow: both cells saturated / above air entry.
    source_flow_water_content =
        source_water_content_current

    destination_flow_water_content =
        destination_water_content_current

    source_conductivity_index =
        conductivity_table_index(
            source_flow_water_content,
            source_porosity
        )

    destination_conductivity_index =
        conductivity_table_index(
            destination_flow_water_content,
            destination_porosity
        )

    source_flow_matric_potential =
        source_current_matric_potential

    destination_flow_matric_potential =
        destination_current_matric_potential

    flow_regime =
        :darcy_both_above_air_entry

elseif source_is_above_air_entry
    # Green-Ampt: source above air entry, destination uses wet front.
    source_flow_water_content =
        source_water_content_current

    destination_flow_water_content =
        bounded_volumetric_water_content(
            destination_wet_front_water,
            destination_reference_volume,
            destination_minimum_water_content,
            destination_porosity
        )

    source_conductivity_index =
        conductivity_table_index(
            source_flow_water_content,
            source_porosity
        )

    destination_conductivity_index =
        conductivity_table_index_airentry(
            destination_flow_water_content,
            destination_porosity,
            destination_air_entry_water_content
        )

    source_flow_matric_potential =
        source_current_matric_potential

    if destination_soil_mass > tiny
        destination_flow_matric_potential,
        destination_flow_water_content =
            cal_matricWaterPotential_MCM(
                destination_flow_water_content,
                destination_field_capacity,
                field_capacity_potential,
                destination_wilting_point,
                wilting_point_potential,
                destination_porosity,
                destination_shape_parameter
            )
    else
        destination_flow_water_content =
            destination_porosity

        destination_flow_matric_potential =
            destination_saturation_potential
    end

    flow_regime =
        :green_ampt_source_above_air_entry

elseif destination_is_above_air_entry
    # Green-Ampt: destination above air entry, source uses wet front.
    source_flow_water_content =
        bounded_volumetric_water_content(
            source_wet_front_water,
            source_reference_volume,
            source_minimum_water_content,
            source_porosity
        )

    destination_flow_water_content =
        destination_water_content_current

    source_conductivity_index =
        conductivity_table_index_airentry(
            source_flow_water_content,
            source_porosity,
            source_air_entry_water_content
        )

    destination_conductivity_index =
        conductivity_table_index(
            destination_flow_water_content,
            destination_porosity
        )

    if source_soil_mass > tiny
        source_flow_matric_potential,
        source_flow_water_content =
            cal_matricWaterPotential_MCM(
                source_flow_water_content,
                source_field_capacity,
                field_capacity_potential,
                source_wilting_point,
                wilting_point_potential,
                source_porosity,
                source_shape_parameter
            )
    else
        source_flow_water_content =
            source_porosity

        source_flow_matric_potential =
            source_saturation_potential
    end

    destination_flow_matric_potential =
        destination_current_matric_potential

    flow_regime =
        :green_ampt_destination_above_air_entry

else
    # Richards flow: neither cell above air entry.
    source_flow_water_content =
        source_water_content_current

    destination_flow_water_content =
        destination_water_content_current

    source_conductivity_index =
        conductivity_table_index(
            source_flow_water_content,
            source_porosity
        )

    destination_conductivity_index =
        conductivity_table_index(
            destination_flow_water_content,
            destination_porosity
        )

    source_flow_matric_potential =
        source_current_matric_potential

    destination_flow_matric_potential =
        destination_current_matric_potential

    flow_regime =
        :richards_neither_above_air_entry
end


# -------------------------------------------------------------------------
# 5. Hydraulic conductivity and interface conductance
# -------------------------------------------------------------------------

source_hydraulic_conductivity =
    hydraulic_conductivity_table[
        boundary_axis_index,
        source_conductivity_index,
        source_layer_index
    ] *
    source_hydraulic_conductivity_multiplier

destination_hydraulic_conductivity =
    hydraulic_conductivity_table[
        boundary_axis_index,
        destination_conductivity_index,
        destination_layer_index
    ]

interface_conductance_value =
    interface_conductance(
        source_hydraulic_conductivity,
        destination_hydraulic_conductivity,
        source_flow_length,
        destination_flow_length
    )


# -------------------------------------------------------------------------
# 6. Total and vapor water potentials
# -------------------------------------------------------------------------

source_total_potential =
    total_water_potential(
        source_flow_matric_potential,
        source_gravity_potential,
        source_osmotic_potential
    )

destination_total_potential =
    total_water_potential(
        destination_flow_matric_potential,
        destination_gravity_potential,
        destination_osmotic_potential
    )

source_vapor_potential =
    vapor_water_potential(
        source_flow_matric_potential,
        source_osmotic_potential
    )

destination_vapor_potential =
    vapor_water_potential(
        destination_flow_matric_potential,
        destination_osmotic_potential
    )


# -------------------------------------------------------------------------
# 7. Potential liquid-water flux
# -------------------------------------------------------------------------

potential_flux =
    water_flux_from_potential(
        source_total_potential,
        destination_total_potential,
        interface_conductance_value,
        exchange_area,
        process_time_factor
    )


# -------------------------------------------------------------------------
# 8. Air-entry correction and storage-limited fluxes
#
# In the deadband:
#     FLQZ = FLQX
#     FLQL = FLQX
#     FLQ2 = FLQX
# -------------------------------------------------------------------------

if abs(potential_flux) < potential_flux_deadband
    corrected_potential_flux =
        potential_flux

    limited_flux =
        potential_flux

    mobile_limited_flux =
        potential_flux

else
    if potential_flux >= 0.0
        if source_flow_water_content > source_air_entry_water_content
            corrected_potential_flux =
                potential_flux +
                min(
                    (source_flow_water_content - source_air_entry_water_content) *
                    source_reference_volume,

                    max(
                        0.0,
                        (destination_air_entry_water_content -
                         destination_flow_water_content) *
                        destination_reference_volume
                    )
                ) *
                storage_removal_fraction
        else
            corrected_potential_flux =
                potential_flux
        end
    else
        if destination_flow_water_content > destination_air_entry_water_content
            corrected_potential_flux =
                potential_flux +
                max(
                    (destination_air_entry_water_content -
                     destination_flow_water_content) *
                    destination_reference_volume,

                    min(
                        0.0,
                        (source_flow_water_content -
                         source_air_entry_water_content) *
                        source_reference_volume
                    )
                ) *
                storage_removal_fraction
        else
            corrected_potential_flux =
                potential_flux
        end
    end

    limited_flux =
        limit_bidirectional_flux(
            corrected_potential_flux,
            source_available_water * storage_removal_fraction,
            destination_air_capacity * storage_removal_fraction,
            destination_available_water * storage_removal_fraction,
            source_air_capacity * storage_removal_fraction
        )

    mobile_limited_flux =
        limit_bidirectional_flux(
            potential_flux,
            source_available_water * storage_removal_fraction,
            destination_air_capacity * storage_removal_fraction,
            destination_available_water * storage_removal_fraction,
            source_air_capacity * storage_removal_fraction
        )
end


# -------------------------------------------------------------------------
# 9. Excess water+ice correction during freezing
#
# Preserve WF8 behavior:
# only apply this correction for vertical interfaces.
# -------------------------------------------------------------------------

excess_water_ice_correction =
    if is_vertical_boundary && destination_excess_water_ice < 0.0
        min(
            0.0,
            max(
                -destination_available_water * storage_removal_fraction,
                destination_excess_water_ice
            )
        )
    else
        0.0
    end

limited_flux +=
    excess_water_ice_correction

mobile_limited_flux +=
    excess_water_ice_correction


# -------------------------------------------------------------------------
# 10. Advective heat and updated temporary available water
# -------------------------------------------------------------------------

advective_heat_flux =
    advective_heat_by_water_flux(
        limited_flux,
        source_temperature,
        destination_temperature,
        cpw
    )

updated_source_available_water =
    source_available_water - limited_flux

updated_destination_available_water =
    destination_available_water + limited_flux


# -------------------------------------------------------------------------
# 11. Return process diagnostics
# -------------------------------------------------------------------------

return (
    flow_regime =
        flow_regime,

    source_water_content_current =
        source_water_content_current,

    destination_water_content_current =
        destination_water_content_current,

    source_flow_water_content =
        source_flow_water_content,

    destination_flow_water_content =
        destination_flow_water_content,

    source_current_matric_potential =
        source_current_matric_potential,

    destination_current_matric_potential =
        destination_current_matric_potential,

    source_flow_matric_potential =
        source_flow_matric_potential,

    destination_flow_matric_potential =
        destination_flow_matric_potential,

    source_total_potential =
        source_total_potential,

    destination_total_potential =
        destination_total_potential,

    source_vapor_potential =
        source_vapor_potential,

    destination_vapor_potential =
        destination_vapor_potential,

    source_conductivity_index =
        source_conductivity_index,

    destination_conductivity_index =
        destination_conductivity_index,

    source_hydraulic_conductivity =
        source_hydraulic_conductivity,

    destination_hydraulic_conductivity =
        destination_hydraulic_conductivity,

    interface_conductance =
        interface_conductance_value,

    potential_flux =
        potential_flux,

    corrected_potential_flux =
        corrected_potential_flux,

    excess_water_ice_correction =
        excess_water_ice_correction,

    limited_flux =
        limited_flux,

    mobile_limited_flux =
        mobile_limited_flux,

    advective_heat_flux =
        advective_heat_flux,

    updated_source_available_water =
        updated_source_available_water,

    updated_destination_available_water =
        updated_destination_available_water
)

end