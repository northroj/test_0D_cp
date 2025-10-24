#pragma once

void simulate();

void simulate_timestep(int t_it);


void transport_particle(class Particle& p, double time_census);

void kp_alpha_csd(class Particle& p);

void spitzer_csd(class Particle& p, double& dedt_electron, double& dedx_electron, double& dedt_ion, double& dedx_ion);


void source_particles(int t_it);

void population_control();