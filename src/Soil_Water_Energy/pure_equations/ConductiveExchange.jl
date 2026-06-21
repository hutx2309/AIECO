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

"""
series_interface_conductance(
conductance_source,
conductance_destination,
length_source,
length_destination;
tiny = tiny_num
)

Effective interface conductance between two adjacent domains connected in series.

Legacy form:

```
G_interface =
    2 * G_source * G_destination /
    (G_source * length_destination + G_destination * length_source)
```

where

```
G_source       = conductance or conductivity-like coefficient in source domain
G_destination  = conductance or conductivity-like coefficient in destination domain
length_source  = source-domain path length
length_destination = destination-domain path length
```

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

"""
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

```
F = G_interface * (X_source - X_destination) * A * Δt * m
```

where

```
F        = exchanged amount over the model substep
G        = interface conductance
X_source = source-side potential, temperature, or concentration
X_destination = destination-side potential, temperature, or concentration
A        = exchange area
Δt       = substep time factor
m        = optional multiplier
```

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

"""
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

```
H = G_T * (T_source - T_destination) * A * Δt * m
```

where

```
H      = heat exchanged over the model substep
G_T    = thermal interface conductance
T      = temperature in K
A      = exchange area
Δt     = substep time factor
m      = optional multiplier
```

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

"""
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

```
F_v = G_v * (C_source - C_destination) * A * Δt * m
```

where

```
F_v = vapor exchanged over the model substep
G_v = vapor interface conductance or diffusivity-like coefficient
C   = vapor concentration
A   = exchange area
Δt  = substep time factor
m   = optional multiplier
```

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