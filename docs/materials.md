```markdown

# [materials]

Defines the physical and compositional properties of the simulation material.

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