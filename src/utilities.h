#pragma once
#include <filesystem>
#include <string>
#include <vector>


// ---------- classes ----------
class Material {
public:
    double ion_temperature      = 0.0;
    double electron_temperature = 0.0;

    // Species store: parallel vectors
    std::vector<std::string> species;          // e.g., {"d","t","a"}
    std::vector<double>      number_densities; // e.g., {1e-20, 1e-20, 0.0}

    void add_species(const std::string& name, double nd) {
        species.push_back(name);
        number_densities.push_back(nd);
    }

    Material() = default;
};

class Particle {
public:
    std::string species;

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    double t = 0.0;

    double energy = 0.0;

    Particle() = default;

    // Easy initializer constructor
    Particle(std::string s,
             double x_, double y_, double z_,
             double t_,
             double energy_ = 0.0)
        : species(std::move(s)), x(x_), y(y_), z(z_), t(t_), energy(energy_) {}
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
