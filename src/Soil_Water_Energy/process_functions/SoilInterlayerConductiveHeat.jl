# =============================================================================
# SoilInterlayerConductiveHeat.jl
#
# Process-level function for conductive heat exchange between adjacent soil
# layers or adjacent soil grid cells.
#
# Design:
#   - pure function
#   - one source-destination soil interface at a time
#   - no mutation
#   - no @unpack
#   - no direct access to HFLWL arrays
# =============================================================================


@doc raw"""
    soil_interlayer_conductive_heat_process(...)

Calculate conductive heat exchange between adjacent soil cells.

```julia
# THERMAL CONDUCTIVITY IN EACH GRID CELL
# HEAT FLOW FROM THERMAL CONDUCTIVITY AND TEMPERATURE GRADIENT
```

The source and destination thermal conductivities are calculated from mineral,
rock, water, ice, and air contributions. The dynamic water and air conductivity
enhancement terms are controlled by the local temperature-gradient scale:

```math
\Delta T_x =
10^{-6} |T_1 - T_2|
```

For each cell, the water and air enhancement factors are

```math
R_w =
\min\left(10^4, r_w \Delta T_x \max(0, \theta_w - \theta_{w,r})^3\right)
```

and

```math
R_a =
\min\left(10^4, r_a \Delta T_x \max(0, \theta_a - \theta_{a,r})^3\right)
```

with corresponding Nusselt-like factors

```math
N_w =
\max\left(1, 0.68 + \frac{0.67 R_w^{1/4}}{d_w}\right)
```

and

```math
N_a =
\max\left(1, 0.68 + \frac{0.67 R_a^{1/4}}{d_a}\right)
```

The effective cell thermal conductivity follows the legacy WF8 mixture form:

```math
k =
\frac{
    k_m + \theta_w k_w + 0.611 \theta_i k_i
    + w_\theta \theta_a k_a
}{
    k_r + \theta_w + 0.611 \theta_i + w_\theta \theta_a
}
```

where

```math
w_\theta = 1.467 - 0.467 \phi_a
```

The interface conductance is calculated as a series conductance:

```math
G_{12} =
\frac{2 k_1 k_2}{k_1 L_2 + k_2 L_1}
```

Before conductive exchange is calculated, the temperatures are adjusted for
vapor heat transport. For the source cell,

```math
T_1' =
T_1 - \frac{H_v}{C_1}
```

or, for the top soil layer without snow cover,

```math
T_1' =
T_1 - \frac{H_v - H_g}{C_1}
```

where ``H_v`` is vapor advective heat and ``H_g`` is ground heat exchange.
The destination cell is adjusted as

```math
T_2' =
T_2 + \frac{H_v}{C_2}
```

The equilibrium temperature is

```math
T_y =
\frac{C_1 T_1' + C_2 T_2'}{C_1 + C_2}
```

The source-side sensible heat availability is

```math
H_x =
(T_1' - T_y) C_1 f_r
```

The potential conductive heat exchange is

```math
H_c =
G_{12} A (T_1' - T_2') \Delta t
```

The final limited conductive heat exchange is

```math
H =
\begin{cases}
\max(0, \min(H_x, H_c)), & H_c \ge 0 \\
\min(0, \max(H_x, H_c)), & H_c < 0
\end{cases}
```

# Sign convention

- `conductive_heat_flux > 0`: heat moves from source to destination.
- `conductive_heat_flux < 0`: heat moves from destination to source.

# Returns

A named tuple with source/destination thermal conductivities, adjusted
temperatures, equilibrium temperature, potential conductive heat, sensible
heat limit, and final limited conductive heat.
"""
function soil_interlayer_conductive_heat_process(;
    # -------------------------------------------------------------------------
    # Source and destination state
    # -------------------------------------------------------------------------
    source_bulk_density,
    destination_bulk_density,

    source_water_volume,
    destination_water_volume,

    source_ice_volume,
    destination_ice_volume,

    source_air_volume,
    destination_air_volume,

    source_air_porosity,
    destination_air_porosity,

    # -------------------------------------------------------------------------
    # Static soil/layer thermal properties
    # -------------------------------------------------------------------------
    source_mineral_thermal_conductivity,
    destination_mineral_thermal_conductivity,

    source_rock_thermal_conductivity,
    destination_rock_thermal_conductivity,

    # -------------------------------------------------------------------------
    # Temperatures and heat capacities
    # -------------------------------------------------------------------------
    source_temperature,
    destination_temperature,

    source_heat_capacity,
    destination_heat_capacity,

    minimum_soil_heat_capacity,

    # -------------------------------------------------------------------------
    # Coupled vapor/surface heat terms
    # -------------------------------------------------------------------------
    vapor_heat_flux,
    ground_heat_flux,

    is_surface_source_layer,
    snow_heat_capacity,
    minimum_snow_heat_capacity,

    # -------------------------------------------------------------------------
    # Geometry and time
    # -------------------------------------------------------------------------
    source_flow_length,
    destination_flow_length,

    exchange_area,
    process_time_factor,
    storage_removal_fraction,

    # -------------------------------------------------------------------------
    # Legacy thermal enhancement factors
    # -------------------------------------------------------------------------
    water_rayleigh_multiplier,
    air_rayleigh_multiplier,

    water_diffusion_denominator,
    air_diffusion_denominator,

    residual_water_volume,
    residual_air_volume,

    # -------------------------------------------------------------------------
    # Numerical constants
    # -------------------------------------------------------------------------
    tiny = tiny_num,
    temperature_gradient_scale = 1.0e-6,
    max_rayleigh_number = 1.0e4,
    water_base_conductivity = 2.067e-3,
    air_base_conductivity = 9.050e-5,
    ice_thermal_conductivity = 7.844e-3,
    ice_geometry_factor = 0.611
)
    # -------------------------------------------------------------------------
    # 1. Temperature-gradient scale
    # -------------------------------------------------------------------------

    temperature_gradient_factor =
        abs(source_temperature - destination_temperature) *
        temperature_gradient_scale


    # -------------------------------------------------------------------------
    # 2. Source thermal conductivity
    # -------------------------------------------------------------------------

    source_is_thermally_active =
        source_bulk_density > tiny ||
        source_water_volume + source_ice_volume > tiny

    if source_is_thermally_active
        source_water_excess =
            max(
                0.0,
                source_water_volume - residual_water_volume
            )^3

        source_air_excess =
            max(
                0.0,
                source_air_volume - residual_air_volume
            )^3

        source_water_rayleigh =
            min(
                max_rayleigh_number,
                water_rayleigh_multiplier *
                temperature_gradient_factor *
                source_water_excess
            )

        source_air_rayleigh =
            min(
                max_rayleigh_number,
                air_rayleigh_multiplier *
                temperature_gradient_factor *
                source_air_excess
            )

        source_water_nusselt =
            max(
                1.0,
                0.68 +
                0.67 *
                source_water_rayleigh^0.25 /
                water_diffusion_denominator
            )

        source_air_nusselt =
            max(
                1.0,
                0.68 +
                0.67 *
                source_air_rayleigh^0.25 /
                air_diffusion_denominator
            )

        source_water_thermal_conductivity =
            water_base_conductivity *
            source_water_nusselt

        source_air_thermal_conductivity =
            air_base_conductivity *
            source_air_nusselt

        source_air_weight =
            1.467 -
            0.467 *
            source_air_porosity

        source_conductivity_numerator =
            source_mineral_thermal_conductivity +
            source_water_volume *
            source_water_thermal_conductivity +
            ice_geometry_factor *
            source_ice_volume *
            ice_thermal_conductivity +
            source_air_weight *
            source_air_volume *
            source_air_thermal_conductivity

        source_conductivity_denominator =
            source_rock_thermal_conductivity +
            source_water_volume +
            ice_geometry_factor *
            source_ice_volume +
            source_air_weight *
            source_air_volume

        source_thermal_conductivity =
            if source_conductivity_denominator > tiny
                source_conductivity_numerator /
                source_conductivity_denominator
            else
                0.0
            end
    else
        source_water_excess = 0.0
        source_air_excess = 0.0
        source_water_rayleigh = 0.0
        source_air_rayleigh = 0.0
        source_water_nusselt = 1.0
        source_air_nusselt = 1.0
        source_water_thermal_conductivity = 0.0
        source_air_thermal_conductivity = 0.0
        source_air_weight = 0.0
        source_thermal_conductivity = 0.0
    end


    # -------------------------------------------------------------------------
    # 3. Destination thermal conductivity
    # -------------------------------------------------------------------------

    destination_is_thermally_active =
        destination_bulk_density > tiny ||
        destination_water_volume + destination_ice_volume > tiny

    if destination_is_thermally_active
        destination_water_excess =
            max(
                0.0,
                destination_water_volume - residual_water_volume
            )^3

        destination_air_excess =
            max(
                0.0,
                destination_air_volume - residual_air_volume
            )^3

        destination_water_rayleigh =
            min(
                max_rayleigh_number,
                water_rayleigh_multiplier *
                temperature_gradient_factor *
                destination_water_excess
            )

        destination_air_rayleigh =
            min(
                max_rayleigh_number,
                air_rayleigh_multiplier *
                temperature_gradient_factor *
                destination_air_excess
            )

        destination_water_nusselt =
            max(
                1.0,
                0.68 +
                0.67 *
                destination_water_rayleigh^0.25 /
                water_diffusion_denominator
            )

        destination_air_nusselt =
            max(
                1.0,
                0.68 +
                0.67 *
                destination_air_rayleigh^0.25 /
                air_diffusion_denominator
            )

        destination_water_thermal_conductivity =
            water_base_conductivity *
            destination_water_nusselt

        destination_air_thermal_conductivity =
            air_base_conductivity *
            destination_air_nusselt

        destination_air_weight =
            1.467 -
            0.467 *
            destination_air_porosity

        destination_conductivity_numerator =
            destination_mineral_thermal_conductivity +
            destination_water_volume *
            destination_water_thermal_conductivity +
            ice_geometry_factor *
            destination_ice_volume *
            ice_thermal_conductivity +
            destination_air_weight *
            destination_air_volume *
            destination_air_thermal_conductivity

        destination_conductivity_denominator =
            destination_rock_thermal_conductivity +
            destination_water_volume +
            ice_geometry_factor *
            destination_ice_volume +
            destination_air_weight *
            destination_air_volume

        destination_thermal_conductivity =
            if destination_conductivity_denominator > tiny
                destination_conductivity_numerator /
                destination_conductivity_denominator
            else
                0.0
            end
    else
        destination_water_excess = 0.0
        destination_air_excess = 0.0
        destination_water_rayleigh = 0.0
        destination_air_rayleigh = 0.0
        destination_water_nusselt = 1.0
        destination_air_nusselt = 1.0
        destination_water_thermal_conductivity = 0.0
        destination_air_thermal_conductivity = 0.0
        destination_air_weight = 0.0
        destination_thermal_conductivity = 0.0
    end


    # -------------------------------------------------------------------------
    # 4. Interlayer thermal conductance
    # -------------------------------------------------------------------------

    thermal_interface_conductance =
        series_interface_conductance(
            source_thermal_conductivity,
            destination_thermal_conductivity,
            source_flow_length,
            destination_flow_length;
            tiny = tiny
        )


    # -------------------------------------------------------------------------
    # 5. Vapor-heat-adjusted source temperature
    #
    # Preserve WF8:
    #   top soil + no snow uses vapor_heat_flux - ground_heat_flux
    #   all other source cells use vapor_heat_flux only
    # -------------------------------------------------------------------------

    if source_heat_capacity > minimum_soil_heat_capacity
        if is_surface_source_layer &&
           snow_heat_capacity <= minimum_snow_heat_capacity

            adjusted_source_temperature =
                source_temperature -
                (vapor_heat_flux - ground_heat_flux) /
                source_heat_capacity
        else
            adjusted_source_temperature =
                source_temperature -
                vapor_heat_flux /
                source_heat_capacity
        end
    else
        adjusted_source_temperature =
            source_temperature
    end


    # -------------------------------------------------------------------------
    # 6. Vapor-heat-adjusted destination temperature
    # -------------------------------------------------------------------------

    if destination_heat_capacity > minimum_soil_heat_capacity
        adjusted_destination_temperature =
            destination_temperature +
            vapor_heat_flux /
            destination_heat_capacity
    else
        adjusted_destination_temperature =
            destination_temperature
    end


    # -------------------------------------------------------------------------
    # 7. Equilibrium temperature and source-side sensible heat limit
    # -------------------------------------------------------------------------

    total_heat_capacity =
        source_heat_capacity +
        destination_heat_capacity

    equilibrium_temperature_value =
        if total_heat_capacity > tiny
            (
                source_heat_capacity *
                adjusted_source_temperature +
                destination_heat_capacity *
                adjusted_destination_temperature
            ) /
            total_heat_capacity
        else
            0.5 *
            (
                adjusted_source_temperature +
                adjusted_destination_temperature
            )
        end

    source_sensible_heat_limit =
        (
            adjusted_source_temperature -
            equilibrium_temperature_value
        ) *
        source_heat_capacity *
        storage_removal_fraction


    # -------------------------------------------------------------------------
    # 8. Potential conductive heat exchange
    # -------------------------------------------------------------------------

    potential_conductive_heat_flux =
        conductive_heat_exchange(
            adjusted_source_temperature,
            adjusted_destination_temperature,
            thermal_interface_conductance,
            exchange_area,
            process_time_factor
        )


    # -------------------------------------------------------------------------
    # 9. Storage-limited conductive heat exchange
    # -------------------------------------------------------------------------

    conductive_heat_flux =
        if potential_conductive_heat_flux >= 0.0
            max(
                0.0,
                min(
                    source_sensible_heat_limit,
                    potential_conductive_heat_flux
                )
            )
        else
            min(
                0.0,
                max(
                    source_sensible_heat_limit,
                    potential_conductive_heat_flux
                )
            )
        end


    # -------------------------------------------------------------------------
    # 10. Return diagnostics
    # -------------------------------------------------------------------------

    return (
        source_is_thermally_active =
            source_is_thermally_active,

        destination_is_thermally_active =
            destination_is_thermally_active,

        temperature_gradient_factor =
            temperature_gradient_factor,

        source_water_rayleigh =
            source_water_rayleigh,

        destination_water_rayleigh =
            destination_water_rayleigh,

        source_air_rayleigh =
            source_air_rayleigh,

        destination_air_rayleigh =
            destination_air_rayleigh,

        source_water_nusselt =
            source_water_nusselt,

        destination_water_nusselt =
            destination_water_nusselt,

        source_air_nusselt =
            source_air_nusselt,

        destination_air_nusselt =
            destination_air_nusselt,

        source_water_thermal_conductivity =
            source_water_thermal_conductivity,

        destination_water_thermal_conductivity =
            destination_water_thermal_conductivity,

        source_air_thermal_conductivity =
            source_air_thermal_conductivity,

        destination_air_thermal_conductivity =
            destination_air_thermal_conductivity,

        source_thermal_conductivity =
            source_thermal_conductivity,

        destination_thermal_conductivity =
            destination_thermal_conductivity,

        thermal_interface_conductance =
            thermal_interface_conductance,

        adjusted_source_temperature =
            adjusted_source_temperature,

        adjusted_destination_temperature =
            adjusted_destination_temperature,

        equilibrium_temperature =
            equilibrium_temperature_value,

        source_sensible_heat_limit =
            source_sensible_heat_limit,

        potential_conductive_heat_flux =
            potential_conductive_heat_flux,

        conductive_heat_flux =
            conductive_heat_flux
    )
end