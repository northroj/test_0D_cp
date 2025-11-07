```markdown

# [tallies]

Specifies diagnostic tallies to record quantities such as energy deposition, particle counts, or flux spectra.

Each `[tallies]` block defines one independent tally.

| Field | Type | Description |
|--------|------|-------------|
| `name` | str | Unique name for the tally |
| `category` | str | Type of measurement (e.g., `csd_energy_loss`, `average_particle_energy`, etc.) |
| `particles` | list[str] | Species to include (e.g., `d t a`) |
| `energy_bins` | float float int mode | Defines energy binning (start end N lin/log/custom) |
| `time_bins` | float float int mode or `timesteps` | Defines time binning or uses global simulation bins |

---

### Example
```ini
[tallies]
name test_tally_1
category csd_energy_loss
particles d t a
energy_bins 0 10000 99 lin
time_bins 0 0.01 99 lin

[tallies]
name test_tally_2
category average_particle_energy
particles a
time_bins timesteps