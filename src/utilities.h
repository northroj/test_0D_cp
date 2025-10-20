#pragma once
#include <filesystem>
#include <string>
#include <vector>

// ---------- Placeholder classes (you'll flesh these out later) ----------
class Material {
public:
    // Thermodynamic-like scalars
    double ion_temperature      = 0.0;
    double electron_temperature = 0.0;

    // Species store: parallel vectors (simple, cache-friendly)
    std::vector<std::string> species;          // e.g., {"d","t","a"}
    std::vector<double>      number_densities; // e.g., {1e-20, 1e-20, 0.0}

    // Optional helpers you can use later
    void add_species(const std::string& name, double nd) {
        species.push_back(name);
        number_densities.push_back(nd);
    }

    // You can add constructors/methods later as needed
    Material() = default;
};

class Surface {
public:
    Surface() = default;
    // TODO: add data/methods
};

class Particle {
public:
    Particle() = default;
    // TODO: add data/methods
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
