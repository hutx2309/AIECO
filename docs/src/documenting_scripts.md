# Documenting scripts

Every Julia source file under `src/` must have a matching Markdown page under
`docs/src/api/`. The directory and base name must match exactly. For example:

```text
src/Soil_Water_Energy/pure_equations/Radiation.jl
docs/src/api/Soil_Water_Energy/pure_equations/Radiation.md
```

The documentation build checks this mapping and stops with a list of missing
pages when a new script is not documented.

## Workflow for a completed script

1. Copy `docs/templates/script.md` to the matching path under `docs/src/api/`.
2. Record the script purpose, equations, units, sign conventions, assumptions,
   bounds, and validation evidence.
3. Add Julia docstrings immediately above public functions, types, and constants.
4. Include the script in the package module. If it is not integrated yet, load
   it in an isolated module in `docs/make.jl`.
5. Add its page to the `pages` navigation in `docs/make.jl`.
6. Build the site and resolve warnings.

## Build locally

From the repository root:

```julia
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Open `docs/build/index.html` after the build. The generated `build/` directory
is ignored by Git.

## Docstring example

````julia
"""
    function_name(input; keyword = default)

One-sentence description.

# Arguments
- `input`: Meaning and units.
- `keyword`: Meaning, units, and allowed range.

# Returns
Meaning, units, shape, and sign convention of the return value.

# Notes
State the governing equation, assumptions, bounds, and edge-case behavior.
"""
function function_name(input; keyword = default)
    # implementation
end
````

