# Example: Full Input File

Below is a complete example input file combining all sections.

```ini

[materials]
ion_temperature 25
electron_temperature 25
particles d t a
densities 5.0 5.0 0.0
mat_id 1

[geometry]
x -0.05 0.2 99 lin
y -1e20 1e20
z -1e20 1e20
fill 1

[source]
particle a
space point 0.0 0.0 0.0
time point 0.0
energy point 3500.0
strength 1e18
direction beam 1.0 0.0 0.0

[tallies]
name test_tally_1
category csd_energy_loss
particles d t a
energy_bins 0 10000 99 lin
time_bins 0 0.01 99 lin
x_bins cells

[tallies]
name test_tally_2
category average_particle_energy
particles a
time_bins timesteps

[settings]
num_particles 1000
timestep_bins 0.0 0.01 99 lin
dimensions 3
csd_step 10
csd_model spitzer


```