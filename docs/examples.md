```markdown

# Example: Full Input File

Below is a complete example input file combining all sections.

```ini
[materials]
ion_temperature 25
electron_temperature 25
particles d t a
densities 5.0 5.0 0.0
mat_id 1

[source]
particle a
point 0.0 0.0 0.0
time 0.0
energy 3500.0
strength 1e18

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

[settings]
num_particles 1000
num_t_steps 100
t_step_size 0.0001