#pragma once
#include <memory>
#include <vector>
#include "utilities.h" // for Material/Surface/Particle
#include <string>
#include <array>

struct Garage {
    // Categories
    std::vector<Material> materials;
    std::vector<Particle> active_bank;
    std::vector<Particle> census_bank;
    std::vector<Particle> secondary_bank;

    // ---- settings & bookkeeping ----
    int    num_particles     = 0;
    int    num_t_steps       = 0;
    double t_step_size       = 0.0;
    int    current_time_step = 0;   // initialize to 0

    // Tallies
    int    standard_tallies  = 0;

    // Source fields
    std::string         source_particle;        // e.g., "a"
    std::array<double,3> source_point{0.0,0.0,0.0}; // x,y,z
    double              source_time = 0.0;
    double              source_energy = 0.0;
};

extern Garage garage; // declaration only