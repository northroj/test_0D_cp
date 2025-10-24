#pragma once
#include <filesystem>
#include <string>
#include <vector>
#include <cstdint>

// ---------- classes ----------
class Material {
public:
    double ion_temperature      = 0.0;
    double electron_temperature = 0.0;
    double ion_molar_mass = 0.0;

    // Species store: parallel vectors
    std::vector<std::string> species;          // e.g., {"d","t","a"}
    std::vector<double>      densities; // e.g., {1e-20, 1e-20, 0.0}

    void add_species(const std::string& name, double nd) {
        species.push_back(name);
        densities.push_back(nd);
    }

    Material() = default;
};

class Particle {
public:
    std::string species;
    
    double weight = 1.0; // ?

    double x = 0.0; // cm
    double y = 0.0;
    double z = 0.0;
    int zone = 0;

    double t = 0.0;  // shk

    double energy = 0.0; // MeV
    double speed = 0.0; // cm/shk

    Particle() = default;

    // Easy initializer constructor
    Particle(const std::string &s,
             double x_, double y_, double z_,
             double t_,
             double energy_,
             double weight_ = 1.0);
};

// ---------- CLI ----------
struct Cli {
    std::filesystem::path input;
    std::filesystem::path output; // may be empty until defaulted
};

// Prints usage text
void print_usage(const char* exe);

// Replace any existing extension on a path with ".h5"
std::filesystem::path with_h5_extension(std::filesystem::path p);

// Parse command line flags into Cli
bool parse_args(int argc, char** argv, Cli& cli);


int32_t species_2_zaid(std::string ion_species);

double species_2_mass(std::string ion_species);

double species_2_molar_mass(std::string ion_species);

double compute_mixture_molar_mass(const std::vector<std::string>& species,
                                  const std::vector<double>& densities);

