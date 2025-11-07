# [materials]

Defines the physical and compositional properties of the simulation material.

| Field | Type | Description |
|--------|------|-------------|
| `ion_temperature` | float (keV) | Initial ion temperature |
| `electron_temperature` | float (keV) | Initial electron temperature |
| `particles` | list[str] | List of ion species (e.g. `d t a`) |
| `densities` | list[float] | Densities (atoms/cc or g/cc, consistent with code) |
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