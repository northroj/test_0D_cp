# [source]

Defines the initial particle source conditions for the simulation.

| Field | Type | Description |
|--------|------|-------------|
| `particle` | str | Particle species to emit (e.g., `a`, `d`, `t`) |
| `point` | 3 floats | Cartesian coordinates of the source in cm |
| `time` | float | Emission time in shakes |
| `energy` | float | Initial particle energy (keV or MeV, depending on code units) |
| `strength` | float | Number of particles emitted (or energy-weighted source strength) |

---

### Example
```ini
[source]
particle a
point 0.0 0.0 0.0
time 0.0
energy 3500.0
strength 1e18