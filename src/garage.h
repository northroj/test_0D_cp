#pragma once
#include <memory>
#include <vector>
#include "utilities.h"
#include <string>
#include <array>
#include <cstddef>

struct Garage {

    // RNG system
    R123Rng rng;
    // optional seed to initialize RNG
    uint64_t rng_seed = 123456789u;

    std::vector<Material> materials;

    // Particle banks
    std::vector<Particle> active_bank;
    std::vector<Particle> census_bank;
    std::vector<Particle> secondary_bank;

    // ---- geometry (cartesian) ----
    std::vector<double> x_bin_bounds; // size nx+1
    std::vector<double> y_bin_bounds; // size ny+1
    std::vector<double> z_bin_bounds; // size nz+1

    // derived after geometry finalize
    int nx = 0, ny = 0, nz = 0;
    std::array<size_t,3> mesh_strides{0,0,0}; // row-major: ix is fastest
    std::vector<MeshCell> mesh_cells;         // size = nx*ny*nz
    std::vector<int>      mesh_mat_ids;       // per-cell mat_id, same size


    // ---- settings & bookkeeping ----
    int    num_particles     = 0;
    int    num_t_steps       = 0;
    double t_step_size       = 0.0;
    int    current_time_step = 0;   // initialize to 0
    std::vector<double> time_step_bins;
    int lost_particles = 0;

    // Tallies
    std::vector<Tally>    tallies;
    std::vector<Tally> standard_tallies;
    int total_particles_created = 0;

    // Source fields
    std::string         source_particle;        // e.g., "a"
    std::array<double,3> source_point{0.0,0.0,0.0}; // x,y,z
    double              source_time = 0.0;
    double              source_energy = 0.0;
    std::array<double,3> source_direction{1.0,0.0,0.0};
    double source_strength = 1; // particles/cc (per second if applicable)
};

extern Garage garage; // declaration only