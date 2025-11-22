# [tallies]

Specifies diagnostic tallies to record quantities such as energy deposition, particle counts, or flux spectra. For the binning parameters, providing 2+ space delimited ascending values should also work. Event triggering is not yet functional.

Each `[tallies]` block defines one independent tally.

| Field | Type | Description |
|--------|------|-------------|
| `name` | str | Unique name for the tally |
| `category` | str | What quantity to tally (e.g., `csd_energy_loss`, `average_particle_energy`, etc.) |
| `event` | str | At what point to trigger the tallying routines. NOT YET IMPLEMENTED |
| `particles` | list[str] | Species to include (e.g., `d t a`) each treated separately |
| `energy_bins` | float float int mode | Defines energy binning (start end N lin/log/custom) |
| `time_bins` | float float int mode or `timesteps` or N floats | Defines time binning or uses global simulation bins |
| `x_bins` | float float int mode or `cells` or N floats | Defines x binning or uses global simulation bins |
| `y_bins` | float float int mode or `cells` or N floats | Defines y binning or uses global simulation bins |
| `z_bins` | float float int mode or `cells` or N floats | Defines z binning or uses the geometry as bins |

<br>

---

Available tally triggering events

| `event` | Description |
|--------|-------------|
| `pathlength` | Triggered after every time a particle is moved |
| `timestep` | Triggered at the end of the timestep (after the census bank is dumped into the active bank) |

<br>

---
 
Available tally categories

| `category` | Description | Recommended Event |
|--------|-----|--------|
| `csd_energy_loss` | Tallies energy lost to continuous slowing down mechanisms | `pathlength` |
| `average_particle_energy` | Tallies the weighted average energy across particles in the active bank (currently only the particles filter works) | `timestep` |


<br>

---

### Example

```ini

[tallies]
name test_tally_1
category csd_energy_loss
particles d t a
energy_bins 0 10000 99 lin
time_bins 0 0.01 99 lin

```

```ini

[tallies]
name test_tally_2
category average_particle_energy
particles a
time_bins timesteps

```