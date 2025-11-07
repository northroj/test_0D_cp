# [settings]

Controls global simulation settings and run-time parameters.

| Field | Type | Description |
|--------|------|-------------|
| `num_particles` | int | Number of Monte Carlo particles to simulate |
| `num_t_steps` | int | Number of discrete time steps |
| `t_step_size` | float | Time step size (shakes) |

---

### Example
```ini
[settings]
num_particles 1000
num_t_steps 100
t_step_size 0.0001