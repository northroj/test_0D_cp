#pragma once

void simulate();

void simulate_timestep(int t_it);


void transport_particle(class Particle& p);


void source_particles(int t_it);

void population_control();