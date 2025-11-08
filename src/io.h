// io.h
#pragma once
#include <filesystem>
#include <vector>
#include <string>

#include "utilities.h"  // NDArray, Tally, etc.
#include "garage.h"     // Garage, globals used by parsing/writing

// Parse the input file into the global garage (defined in main.cpp).
bool parse_input_file(const std::filesystem::path& path);

// Reserve particle-bank capacity based on num_particles in garage.
void plan_particle_capacity_after_parse();

// finalize geometry into structured mesh (cells, neighbors, bounds, materials)
bool build_cartesian_mesh_after_parse();

// Write an HDF5 output file with input, tallies, and timing info.
bool write_output(const std::filesystem::path& out, double sim_s, double total_s);


size_t mesh_flatten(const Garage& g, int ix, int iy, int iz);
void   mesh_unflatten(const Garage& g, size_t id, int &ix, int &iy, int &iz);
bool   mesh_in_bounds(const Garage& g, int ix, int iy, int iz);
int    mesh_neighbor_id(const Garage& g, int cell_id, int a, int s);