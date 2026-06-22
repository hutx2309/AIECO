# ============================================================
# Hydraulics.jl
#
# Pure hydraulic helper functions for water movement.
#
# These functions should:
#   - not mutate model state
#   - not access water variables directly
#   - not unpack structs
#   - return values only
#
# ============================================================


"""
    bounded_volumetric_water_content(water_storage, reference_volume, θ_min, θ_max)

Convert water storage to bounded volumetric water content.

Equation:
    θ = water_storage / reference_volume

Inputs:
    water_storage    : water volume or water-equivalent storage
    reference_volume : bulk/control volume used to define volumetric water content
    θ_min            : lower bound
    θ_max            : upper bound, usually porosity

Output:
    θ : bounded volumetric water content [m³ water m⁻³ bulk volume]

Notes:
    `reference_volume` is usually soil bulk volume or litter reference volume. 
    If the denominator were pore volume, the result would be
    saturation degree, not volumetric water content.
"""
function bounded_volumetric_water_content(
    water_storage,
    reference_volume,
    θ_min,
    θ_max
)
    if reference_volume > tiny_num2
        return max(θ_min, min(θ_max, water_storage / reference_volume))
    else
        return θ_max
    end
end


"""
    conductivity_table_index(θ, porosity)

Convert volumetric water content θ to a hydraulic-conductivity lookup-table index.

The table uses 1:100 classes based on air-filled fraction:

    index = floor(100 * max(0, porosity - θ) / porosity) + 1

Inputs:
    θ        : volumetric water content [m³ water m⁻³ bulk volume]
    porosity : porosity [m³ pore m⁻³ bulk volume]

Output:
    integer index in 1:100
"""
function conductivity_table_index(θ, porosity)
    if porosity > tiny_num
        air_fraction = max(0.0, porosity - θ) / porosity
        return clamp(floor(Int, 100.0 * air_fraction) + 1, 1, 100)
    else
        return 100
    end
end


"""
    conductivity_table_index_airentry(θ, porosity, θ_airentry)

Variant used in Green-Ampt-type cases where water content is capped
by the air-entry water content before selecting conductivity class.

This reproduces patterns like:

    floor(100 * (porosity - min(θ_airentry, θ)) / porosity) + 1
"""
function conductivity_table_index_airentry(θ, porosity, θ_airentry)
    θ_eff = min(θ_airentry, θ)
    return conductivity_table_index(θ_eff, porosity)
end


"""
    total_water_potential(ψ_matric, ψ_gravity, ψ_osmotic)

Total water potential used to drive liquid water flow.

Equation:
    Ψ_total = ψ_matric + ψ_gravity + ψ_osmotic
"""
total_water_potential(ψ_matric, ψ_gravity, ψ_osmotic) =
    ψ_matric + ψ_gravity + ψ_osmotic


"""
    vapor_water_potential(ψ_matric, ψ_osmotic)

Water potential used for vapor-equilibrium calculations.
Gravity is omitted because local vapor pressure and freezing point depression depend
on local water chemical potential, not vertical gravitational potential.

Equation:
    Ψ_vapor = ψ_matric + ψ_osmotic
"""
vapor_water_potential(ψ_matric, ψ_osmotic) =
    ψ_matric + ψ_osmotic


"""
    interface_conductance(K_source, K_dest, length_source, length_dest; tiny=tiny_num)

Effective conductance between two adjacent layers or compartments.

This is the harmonic-style conductance for two conductors in series:

    C = 2 K_source K_dest /
        (K_source * length_dest + K_dest * length_source)

Inputs:
    K_source      : conductivity of source compartment
    K_dest        : conductivity of destination compartment
    length_source : flow-path length / thickness of source compartment
    length_dest   : flow-path length / thickness of destination compartment

Output:
    effective interface conductance
"""
function interface_conductance(
    K_source,
    K_dest,
    length_source,
    length_dest;
    tiny = tiny_num
)
    if (K_source > tiny) && (K_dest > tiny)
        denom = K_source * length_dest + K_dest * length_source
        return denom > tiny ? 2.0 * K_source * K_dest / denom : 0.0
    else
        return 0.0
    end
end


"""
    water_flux_from_potential(Ψ_source, Ψ_dest, conductance, active_area, dt)

Raw potential-driven water flux between two compartments.

Positive flux means source -> destination.

Equation:
    q = C * (Ψ_source - Ψ_dest) * A * f

where:
    C           = interface conductance
    Ψ_source    = total water potential of source compartment
    Ψ_dest      = total water potential of destination compartment
    active_area = exchange area
    dt          = time step 
"""
function water_flux_from_potential(
    Ψ_source,
    Ψ_dest,
    conductance,
    active_area,
    dt
)
    return conductance * (Ψ_source - Ψ_dest) * active_area * dt
end

"""
    advective_heat_by_water_flux(q_water, T_source, T_dest, cpw)

Heat carried by a water flux.

Positive q_water means:
    source -> destination

Negative q_water means:
    destination -> source

The heat flux uses the temperature of the compartment where the moving
water originates.
"""
function advective_heat_by_water_flux(q_water, T_source, T_dest, cpw)
    if q_water > 0.0
        return cpw * T_source * q_water
    elseif q_water < 0.0
        return cpw * T_dest * q_water
    else
        return 0.0
    end
end


"""
    advective_heat_by_vapor_flux(q_vapor, T_source, T_dest, cpw)

Heat carried by a vapor flux.

Positive q_vapor means:
    source -> destination

Negative q_vapor means:
    destination -> source

The heat flux uses the temperature of the compartment where the moving
vapor originates.

Note:
    This is sensible/advective heat carried by vapor mass, not latent heat.
    Latent heat from evaporation/condensation should be handled separately.
"""
function advective_heat_by_vapor_flux(q_vapor, T_source, T_dest, cpw)
    if q_vapor > 0.0
        return cpw * T_source * q_vapor
    elseif q_vapor < 0.0
        return cpw * T_dest * q_vapor
    else
        return 0.0
    end
end