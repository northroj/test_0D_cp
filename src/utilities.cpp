#include "utilities.h"

#include <iostream>
#include <cstdlib>
#include <string>
#include <cmath>

// HDF5 C API is a C library; wrap in extern "C" for C++ builds
extern "C" {
#include <hdf5.h>
}


// Constructor definition
Particle::Particle(const std::string &s,
                   double x_, double y_, double z_,
                   double t_,
                   double energy_,
                   double weight_)
    : species(s),
      weight(weight_),
      x(x_), y(y_), z(z_), t(t_),
      energy(energy_) {

    double projectile_mass = species_2_mass(species);
    speed = std::sqrt(2.0 * energy_ / (projectile_mass * 1.0e-3)) * 1.0e-8 * 1.0e2; // cm/shk
}



// ------------------------ helpers ------------------------

void print_usage(const char* exe) {
    std::cerr
        << "Usage:\n"
        << "  " << exe << " -i <input_file> [-o <output_path>]\n\n"
        << "Examples:\n"
        << "  " << exe << " -i path/to/input/input_file.txt -o path/to/output/output_file.h5\n"
        << "  " << exe << " -i input_file.txt\n"
        << "  " << exe << " -i input_file.txt -o output_file   (auto .h5)\n";
}

std::filesystem::path with_h5_extension(std::filesystem::path p) {
    p.replace_extension(".h5");
    return p;
}

bool parse_args(int argc, char** argv, Cli& cli) {
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-i" || a == "--input") {
            if (i + 1 >= argc) {
                std::cerr << "ERROR: Missing value after " << a << "\n";
                return false;
            }
            cli.input = std::filesystem::path(argv[++i]);
        } else if (a == "-o" || a == "--output") {
            if (i + 1 >= argc) {
                std::cerr << "ERROR: Missing value after " << a << "\n";
                return false;
            }
            cli.output = std::filesystem::path(argv[++i]);
        } else if (a == "-h" || a == "--help") {
            print_usage(argv[0]);
            std::exit(EXIT_SUCCESS);
        } else {
            std::cerr << "WARN: Unrecognized argument ignored: " << a << "\n";
        }
    }

    if (cli.input.empty()) {
        std::cerr << "ERROR: No input file provided. Use -i <input_file>.\n";
        return false;
    }

    // If -o not provided: default to input path with .h5 extension.
    // If -o provided: force .h5 extension regardless of user suffix.
    if (cli.output.empty()) {
        cli.output = with_h5_extension(cli.input);
    } else {
        cli.output = with_h5_extension(cli.output);
    }

    return true;
}


int32_t species_2_zaid(std::string ion_species){
    int32_t ion_zaid;
    if (ion_species == "d") {
        ion_zaid = 1002;
    } else if (ion_species == "t") {
        ion_zaid = 1003;
    } else if (ion_species == "a") {
        ion_zaid = 2004;
    }
    return ion_zaid;
}

double species_2_mass(std::string ion_species){
    double ion_mass;
    if (ion_species == "d") {
        ion_mass = 3.34358e-24; // [g/atom]
    } else if (ion_species == "t") {
        ion_mass = 5.00736e-24;
    } else if (ion_species == "a") {
        ion_mass = 9.10938291e-28;
    }
    return ion_mass;
}

double species_2_molar_mass(std::string ion_species){
    double ion_mass;
    if (ion_species == "d") {
        ion_mass = 2.013553214; // g/mol (deuterium)
    } else if (ion_species == "t") {
        ion_mass = 3.016049;    // g/mol (tritium)
    } else if (ion_species == "a") {
        ion_mass = 4.002602;    // g/mol (helium-4)
    } else if (ion_species == "p") {
        ion_mass = 1.007276;    // g/mol (proton)
    } else {
        std::cerr << "WARN: Unknown ion species \"" << ion_species
                  << "\". Returning 0 molar mass.\n";
    }
    return ion_mass;
}


double compute_mixture_molar_mass(const std::vector<std::string>& species,
                                  const std::vector<double>& densities)
{
    if (species.size() != densities.size() || species.empty()) {
        std::cerr << "ERROR: compute_mixture_molar_mass: species and density vectors mismatch.\n";
        return 0.0;
    }

    double sum_rho = 0.0;
    double sum_rho_over_M = 0.0;

    for (size_t i = 0; i < species.size(); ++i) {
        const double rho_i = densities[i];
        const double M_i   = species_2_molar_mass(species[i]);

        if (M_i <= 0.0) continue; // skip unknowns
        sum_rho += rho_i;
        sum_rho_over_M += rho_i / M_i;
    }

    if (sum_rho_over_M == 0.0) {
        std::cerr << "ERROR: compute_mixture_molar_mass: zero denominator.\n";
        return 0.0;
    }

    return sum_rho / sum_rho_over_M; // g/mol
}