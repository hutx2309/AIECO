# =============================================================================
# Pure transport helper equations for water-energy exchange.

# Scope:
# - interface conductance between two adjacent domains
# - conductive heat exchange
# - diffusive vapor exchange
# - generic gradient-driven exchange

# Design rules:
# - no mutation
# - no @unpack
# - no direct access to waterVar_copy
# - no process-level state update
# - functions return algebraic transport quantities only
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Interface conductance
# -----------------------------------------------------------------------------

@doc raw"""
series_interface_conductance(
conductance_source,
conductance_destination,
length_source,
length_destination;
tiny = tiny_num
)

Effective interface conductance between two adjacent domains connected in series.

Legacy form:

```math
G_{interface} =
\frac{2G_sG_d}{G_sL_d + G_dL_s}
```

where `G_s` and `G_d` are source and destination conductances, and `L_s` and
`L_d` are source and destination path lengths.

This helper is suitable for thermal conduction, vapor diffusion, and hydraulic
conductance when the same harmonic/series averaging structure is used.

For hydraulic-only code, you may keep a more domain-specific wrapper in
Hydraulics.jl. This generic version is useful when the same equation appears
for heat and vapor transport.
"""
function series_interface_conductance(
    conductance_source,
    conductance_destination,
    length_source,
    length_destination;
    tiny = tiny_num
    )
    if conductance_source > tiny && conductance_destination > tiny
        denominator = conductance_source * length_destination +
        conductance_destination * length_source

        if denominator > tiny
            return 2.0 *
                conductance_source *
                conductance_destination /
                denominator
        else
            return 0.0
        end
    else
        return 0.0
    end
end

# -----------------------------------------------------------------------------
# 2. Generic gradient-driven exchange
# -----------------------------------------------------------------------------

@doc raw"""
gradient_exchange(
source_value,
destination_value,
interface_conductance,
area,
time_factor;
multiplier = 1.0
)

Generic exchange amount driven by a difference between two domains.

Physical form:

```math
F = G_{interface}(X_s - X_d)A\Delta t\,m
```

where `F` is the exchanged amount, `G_{interface}` is the interface conductance,
`X_s - X_d` is the driving gradient, `A` is exchange area, `\Delta t` is the
substep time factor, and `m` is an optional multiplier.

Sign convention:
positive F means movement from source to destination.
negative F means movement from destination to source.

"""
function gradient_exchange(
    source_value,
    destination_value,
    interface_conductance,
    area,
    time_factor;
    multiplier = 1.0
    )
    return interface_conductance *
    (source_value - destination_value) *
    area *
    time_factor *
    multiplier
end

# -----------------------------------------------------------------------------
# 3. Conductive heat exchange
# -----------------------------------------------------------------------------

@doc raw"""
conductive_heat_exchange(
source_temperature,
destination_temperature,
thermal_interface_conductance,
area,
time_factor;
multiplier = 1.0
)

Conductive heat exchange between two adjacent domains.

Physical form:

```math
H = G_T(T_s - T_d)A\Delta t\,m
```

where `H` is heat exchanged over the model substep, `G_T` is thermal interface
conductance, `T_s - T_d` is the temperature difference, `A` is exchange area,
`\Delta t` is the substep time factor, and `m` is an optional multiplier.

Sign convention:
positive H means heat moves from source to destination.
negative H means heat moves from destination to source.
"""
function conductive_heat_exchange(
    source_temperature,
    destination_temperature,
    thermal_interface_conductance,
    area,
    time_factor;
    multiplier = 1.0
    )
    return thermal_interface_conductance *
           (source_temperature - destination_temperature) *
           area *
           time_factor *
           multiplier
end

# -----------------------------------------------------------------------------
# 4. Diffusive vapor exchange
# -----------------------------------------------------------------------------

@doc raw"""
diffusive_vapor_exchange(
source_vapor_concentration,
destination_vapor_concentration,
vapor_interface_conductance,
area,
time_factor;
multiplier = 1.0
)

Diffusive vapor exchange between two adjacent domains.

Physical form:

```math
F_v = G_v(C_s - C_d)A\Delta t\,m
```

where `F_v` is vapor exchanged over the model substep, `G_v` is vapor interface
conductance or a diffusivity-like coefficient, `C_s - C_d` is the vapor
concentration difference, `A` is exchange area, `\Delta t` is the substep time
factor, and `m` is an optional multiplier.

Sign convention:
positive F_v means vapor moves from source to destination.
negative F_v means vapor moves from destination to source.

This function calculates vapor movement only. It does not calculate latent heat.
If vapor exchange causes condensation or evaporation, latent heat should be
calculated separately using `latent_heat_from_phase_change(...)`.
"""
function diffusive_vapor_exchange(
    source_vapor_concentration,
    destination_vapor_concentration,
    vapor_interface_conductance,
    area,
    time_factor;
    multiplier = 1.0
    )
    return vapor_interface_conductance *
           (source_vapor_concentration - destination_vapor_concentration) *
           area *
           time_factor *
           multiplier
end
