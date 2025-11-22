# [source]

Defines the initial particle source conditions for the simulation. Currently, only a single monoenergetic, pulse, point, beam source of a single particle type is allowed.

| Field | Type | Description |
|--------|------|-------------|
| `particle` | str | Particle species to emit (e.g., `a`, `d`, `t`) |
| `space` | string+ | The first entry after space is the category which can be followed by numeric specifics  |
| `time` | string+ | The first entry after time is the category which can be followed by numeric specifics |
| `energy` | string+ | The first entry after energy is the category which can be followed by numeric specifics |
| `strength` | float | Number of particles emitted (or energy-weighted source strength) |
| `direction` | string+ | The first entry after direction is the category which can be followed by numeric specifics |

 <br/><br/>

| Space category | Type | Description |
|--------|------|-------------|
| `point` | 3 floats | Cartesian coordinates of the point source in cm |

 <br/><br/>

| Time category | Type | Description |
|--------|------|-------------|
| `point` | float | Emission time in shakes |

 <br/><br/>

| Energy category | Type | Description |
|--------|------|-------------|
| `point` | float | Initial particle energy (keV depending on code units) |

 <br/><br/>

| Direction category | Type | Description |
|--------|------|-------------|
| `isotropic` | N/A | No subsequent tokens |
| `beam` | 3 floats | An x,y,z unit vector for the direction (normalized internally) |


---

### Example

```ini

[source]
particle a
space point 0.0 0.0 0.0
time point 0.0
energy point 3500.0
strength 1e18
direction beam 1.0 0.0 0.0

```