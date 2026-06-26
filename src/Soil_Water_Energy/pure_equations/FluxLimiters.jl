# =============================================================================
# FluxLimiters.jl

# Pure limiter helper equations for water, vapor, and energy exchange.

#
# Scope:
# - limit positive flux by source availability and destination capacity
# - limit negative flux by reverse-source availability and reverse-destination capacity
# - provide readable wrappers for bidirectional exchange limits
#
# Design rules:
# - no mutation
# - no @unpack
# - no direct access to waterVar_copy
# - no process-level state update
# - functions only constrain candidate fluxes
#

# Sign convention:
# q > 0 means movement from source to destination
# q < 0 means movement from destination back to source
# =============================================================================

# -----------------------------------------------------------------------------
# 1. One-direction positive flux limiter
# -----------------------------------------------------------------------------

@doc raw"""
limit_positive_flux(
potential_flux,
source_available,
destination_capacity
)

Limit a positive flux from source to destination.

Physical meaning: for `q > 0`, the source loses material or energy and the
destination gains it.

The actual flux cannot exceed:
1. the potential positive flux
2. source available storage
3. destination available capacity

Legacy form:

```math
q_{limited} = \max\left(0, \min(q_{potential}, S_{avail}, D_{cap})\right)
```

If potential_flux is negative, this function returns 0.
"""
function limit_positive_flux(
    potential_flux,
    source_available,
    destination_capacity
    )
    return max(
        0.0,
        min(
            potential_flux,
            max(0.0, source_available),
            max(0.0, destination_capacity)
        )
    )
end

# -----------------------------------------------------------------------------
# 2. One-direction negative flux limiter
# -----------------------------------------------------------------------------

@doc raw"""
limit_negative_flux(
potential_flux,
destination_available,
source_capacity
)

Limit a negative flux from destination back to source.

Physical meaning: for `q < 0`, the destination loses material or energy and
the source gains it.

The actual negative flux cannot be more negative than:
1. the potential negative flux
2. available storage in the destination
3. available capacity in the source

Legacy form:

```math
q_{limited} = \min\left(0, \max(q_{potential}, -D_{avail}, -S_{cap})\right)
```

If potential_flux is positive, this function returns 0.
"""
function limit_negative_flux(
    potential_flux,
    destination_available,
    source_capacity
    )
    return min(
        0.0,
        max(
        potential_flux,
        -max(0.0, destination_available),
        -max(0.0, source_capacity)
        )
    )
end

# -----------------------------------------------------------------------------
# 3. Bidirectional storage-capacity limiter
# -----------------------------------------------------------------------------

@doc raw"""
limit_bidirectional_flux(
potential_flux,
source_available,
destination_capacity,
destination_available,
source_capacity
)

Limit a bidirectional flux using source availability and destination capacity.

Sign convention: `potential_flux > 0` moves source to destination;
`potential_flux < 0` moves destination to source.

For positive flux:

```math
q_{limited} = \operatorname{limit\_positive}(q_p, S_{avail}, D_{cap})
```

For negative flux:

```math
q_{limited} = \operatorname{limit\_negative}(q_p, D_{avail}, S_{cap})
```

This helper is useful when the same potential flux can move in either
direction, such as capillary exchange or vertical soil water redistribution.
"""
function limit_bidirectional_flux(
    potential_flux,
    source_available,
    destination_capacity,
    destination_available,
    source_capacity
    )
    if potential_flux > 0.0
    return limit_positive_flux(
    potential_flux,
    source_available,
    destination_capacity
    )
    elseif potential_flux < 0.0
    return limit_negative_flux(
    potential_flux,
    destination_available,
    source_capacity
    )
    else
    return 0.0
    end
end

# -----------------------------------------------------------------------------
# 5. Evaporation-style negative source limiter
# -----------------------------------------------------------------------------

@doc raw"""
limit_negative_flux_by_source(
potential_flux,
source_available
)

Limit a negative flux by available source storage.

This is useful for evaporation-like terms where negative flux removes water
from a storage pool.

Legacy form:

```math
q_{limited} = \max(q_{potential}, -S_{avail})
```

If potential_flux is positive, this function returns the positive value
unchanged. If you want condensation to be allowed but evaporation limited,
this helper can be used directly.
"""
function limit_negative_flux_by_source(
    potential_flux,
    source_available
    )
    return max(
    potential_flux,
    -max(0.0, source_available)
    )
end

# -----------------------------------------------------------------------------
# 7. Fractional removal limiter
# -----------------------------------------------------------------------------

@doc raw"""
removable_storage(
storage,
removal_fraction
)

Calculate the amount of a storage pool removable during a model substep.

Legacy pattern:

```math
S_{removable} = \max(0, S f_r)
```

where XNPXX, XNPAX, etc. represent substep removal fractions.

This helper is intentionally simple, but it makes flux limiters easier to read.
"""
function removable_storage(
    storage,
    removal_fraction
    )
    return max(0.0, storage * removal_fraction)
end

@doc raw"""
limit_negative_flux_by_fractional_storage(
potential_flux,
storage,
removal_fraction
)

Limit a negative flux by a fraction of available storage.

Legacy form:

```math
q_{limited} = \max\left(q_{potential}, -\max(0, S f_r)\right)
```

This is useful for evaporation and vapor-condensation terms.
"""
function limit_negative_flux_by_fractional_storage(
    potential_flux,
    storage,
    removal_fraction
    )
    return limit_negative_flux_by_source(
        potential_flux,
        removable_storage(storage, removal_fraction)
    )
end
