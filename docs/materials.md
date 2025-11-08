# [materials]

Defines the initial state of a material that can be used to fill a mesh cell. The [materials] block is repeatable to define multiple materials.

| Field | Type | Description |
|--------|------|-------------|
| `ion_temperature` | float | Initial ion temperature (keV) |
| `electron_temperature` | float | Initial electron temperature (keV) |
| `particles` | list[str] | List of ion species (e.g. `d t a`) |
| `densities` | list[float] | Densities of the various ion species (g/cc in the order of `particles`) |
| `mat_id` | int | Material ID (unique) |

---

### Example

```ini

[materials]
ion_temperature 25
electron_temperature 25
particles d t a
densities 5.0 5.0 0.0
mat_id 1

```
