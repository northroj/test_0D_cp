#pragma once
#include <string>
#include <vector>
#include <unordered_map>

void simulate();

void simulate_timestep(int t_it);


void transport_particle(class Particle& p, double time_census);

void kp_alpha_csd(class Particle& p);

void spitzer_csd(class Particle& p, double& dedt_electron, double& dedx_electron, double& dedt_ion, double& dedx_ion);

std::unordered_map<std::string, double> average_energy_by_species(const std::vector<Particle>& bank);

void source_particles(double time_start, double time_census);

void population_control();

void time_step_setup(std::string mode);