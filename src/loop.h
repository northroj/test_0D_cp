#pragma once
#include <string>
#include <vector>
#include <unordered_map>

void simulate();

void simulate_timestep(int t_it);


void transport_particle(class Particle& p, double time_census);

// Face order convention: 0:-x, 1:+x, 2:-y, 3:+y, 4:-z, 5:+z
struct BoundaryCross {
    double distance;  // > 0 distance to the hit plane (cm)
    int    face;      // 0..5 (see above), or -1 if no hit
    int    next_zone; // neighbor cell id, or boundary code (-1 vacuum, -2 reflective)
};

// Returns true if a boundary was found; false if particle not in a valid cell or no forward hit.
bool distance_to_boundary(Particle& p, BoundaryCross& out);

void nudge_into_cell(double &coord, double lo, double hi, double raw_eps);

void spitzer_csd(class Particle& p, double& dedt_electron, double& dedx_electron, double& dedt_ion, double& dedx_ion);

std::unordered_map<std::string, double> average_energy_by_species(const std::vector<Particle>& bank);

void source_particles(double time_start, double time_census);

void population_control();

void time_step_setup(std::string mode);

// Returns flat cell id in garage.mesh_cells that contains (x,y,z), or -1 if out of bounds.
int locate_cell_id(double x, double y, double z);

// Variant that also returns (ix, iy, iz). Returns false if out of bounds.
bool locate_cell_indices(double x, double y, double z, int &ix, int &iy, int &iz);