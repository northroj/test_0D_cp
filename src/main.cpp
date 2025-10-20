// main.cpp
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <filesystem>
#include <cstdlib>
#include <memory>
#include <fstream>
#include <sstream>
#include <array>
#include <cctype>
#include <algorithm>

// Toolchain check header; not used yet
#include "cdi_CPEloss/Tabular_CP_Eloss.hh"

// HDF5 C API is a C library; wrap in extern "C" for C++ builds
extern "C" {
#include <hdf5.h>
}

// Local headers
#include "utilities.h"
#include "loop.h"


// UNITS ----------------------------------------------------------------------
// temperature keV
// distance cm
// time shk
// number_density atoms/cc



// ---- Global container: garage ----------------------------------------------
struct Garage {
    // Categories
    std::vector<std::shared_ptr<Material>> materials;
    std::vector<std::shared_ptr<Surface>>  surfaces;
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
};

Garage garage;


// ------------------------ Input parsing --------------------------------------

static std::string trim(std::string s) {
    auto isws = [](unsigned char c){ return std::isspace(c); };
    while (!s.empty() && isws(s.front())) s.erase(s.begin());
    while (!s.empty() && isws(s.back()))  s.pop_back();
    return s;
}

static std::vector<std::string> split_ws(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> tok;
    std::string t;
    while (iss >> t) tok.push_back(t);
    return tok;
}

// Parse the simple INI-like file described by the user and populate `garage`.
// Also constructs one Material for the [materials] block and pushes it to garage.materials.
static bool parse_input_file(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "ERROR: Cannot open input file: " << path << "\n";
        return false;
    }

    enum class Section { None, Materials, Source, Tallies, Settings };
    Section sec = Section::None;

    // We'll build a single Material instance from the [materials] block
    auto mat = std::make_shared<Material>();

    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;

        // section headers
        if (line.front() == '[' && line.back() == ']') {
            std::string name = line.substr(1, line.size()-2);
            if      (name == "materials") sec = Section::Materials;
            else if (name == "source")    sec = Section::Source;
            else if (name == "tallies")   sec = Section::Tallies;
            else if (name == "settings")  sec = Section::Settings;
            else                          sec = Section::None;
            continue;
        }

        // content lines
        auto tok = split_ws(line);
        if (tok.empty()) continue;

        switch (sec) {
        case Section::Materials: {
            // ion_temperature 25
            // electron_temperature 25
            // particles d t a
            // number_densities 1e-20 1e-20 0.0
            if (tok[0] == "ion_temperature" && tok.size() == 2) {
                mat->ion_temperature = std::stod(tok[1]);
            } else if (tok[0] == "electron_temperature" && tok.size() == 2) {
                mat->electron_temperature = std::stod(tok[1]);
            } else if (tok[0] == "particles" && tok.size() >= 2) {
                mat->species.assign(tok.begin()+1, tok.end());
            } else if (tok[0] == "number_densities" && tok.size() >= 2) {
                mat->number_densities.clear();
                mat->number_densities.reserve(tok.size()-1);
                for (size_t i = 1; i < tok.size(); ++i) {
                    mat->number_densities.push_back(std::stod(tok[i]));
                }
            } else {
                std::cerr << "WARN: Unrecognized materials line: " << line << "\n";
            }
        } break;

        case Section::Source: {
            // particle a
            // point 0.0 0.0 0.0
            // time 0.0
            if (tok[0] == "particle" && tok.size() == 2) {
                garage.source_particle = tok[1];
            } else if (tok[0] == "point" && tok.size() == 4) {
                garage.source_point[0] = std::stod(tok[1]);
                garage.source_point[1] = std::stod(tok[2]);
                garage.source_point[2] = std::stod(tok[3]);
            } else if (tok[0] == "time" && tok.size() == 2) {
                garage.source_time = std::stod(tok[1]);
            } else {
                std::cerr << "WARN: Unrecognized source line: " << line << "\n";
            }
        } break;

        case Section::Tallies: {
            // standard_tallies 1
            if (tok[0] == "standard_tallies" && tok.size() == 2) {
                garage.standard_tallies = std::stoi(tok[1]);
            } else {
                std::cerr << "WARN: Unrecognized tallies line: " << line << "\n";
            }
        } break;

        case Section::Settings: {
            // num_particles 100
            // num_t_steps 10
            // t_step_size 1
            if (tok[0] == "num_particles" && tok.size() == 2) {
                garage.num_particles = std::stoi(tok[1]);
            } else if (tok[0] == "num_t_steps" && tok.size() == 2) {
                garage.num_t_steps = std::stoi(tok[1]);
            } else if (tok[0] == "t_step_size" && tok.size() == 2) {
                garage.t_step_size = std::stod(tok[1]);
            } else {
                std::cerr << "WARN: Unrecognized settings line: " << line << "\n";
            }            
        } break;

        case Section::None:
        default:
            std::cerr << "WARN: Content outside recognized section: " << line << "\n";
            break;
        }
    }

    // Basic consistency check for materials
    if (!mat->species.empty() || !mat->number_densities.empty()
        || mat->ion_temperature != 0.0 || mat->electron_temperature != 0.0)
    {
        if (mat->species.size() != mat->number_densities.size()) {
            std::cerr << "ERROR: materials: species count (" << mat->species.size()
                      << ") != number_densities count (" << mat->number_densities.size() << ")\n";
            return false;
        }
        garage.materials.push_back(std::move(mat));
    }

    return true;
}

static void plan_particle_capacity_after_parse() {
    // clamp to >= 0 and convert to size_t
    const size_t n = static_cast<size_t>(std::max(garage.num_particles, 0));
    const size_t headroom = n + n / 2; // 1.5x

    garage.active_bank.reserve(n);     // likely equals the active working set
    garage.census_bank.reserve(headroom);
    garage.secondary_bank.reserve(headroom);
}


bool write_output(const std::filesystem::path& out) {
    hid_t file = H5Fcreate(out.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        std::cerr << "ERROR: Failed to create HDF5 file: " << out << "\n";
        return false;
    }

    auto mkgroup = [&](const char* name) -> bool {
        hid_t g = H5Gcreate2(file, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (g < 0) return false;
        if (H5Gclose(g) < 0) return false;
        return true;
    };

    bool ok = true;
    ok = ok && mkgroup("/input");
    ok = ok && mkgroup("/tallies");
    ok = ok && mkgroup("/details");

    if (H5Fclose(file) < 0) {
        std::cerr << "ERROR: Failed to close HDF5 file: " << out << "\n";
        return false;
    }
    if (!ok) {
        std::cerr << "ERROR: Failed creating one or more HDF5 groups in: " << out << "\n";
        return false;
    }
    return true;
}




// ---- main -------------------------------------------------------------------

int main(int argc, char** argv) try {
    using clock_t = std::chrono::steady_clock;

    std::cout << "Preparing simulation\n";

    // Start total timer
    const auto t0_total = clock_t::now();

    // Parse CLI
    Cli cli;
    if (!parse_args(argc, argv, cli)) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    
    // Parse input file (fills garage + creates one Material)
    if (!parse_input_file(cli.input)) {
        return EXIT_FAILURE;
    }
    // allocate particle banks
    plan_particle_capacity_after_parse();

    std::cout << "Running <insert_code_name> simulation\n";
    const auto t0_sim = clock_t::now();
    simulate(cli.input);
    const auto t1_sim = clock_t::now();

    // Write output HDF5 skeleton
    if (!write_output(cli.output)) {
        return EXIT_FAILURE;
    }

    // End total timer
    const auto t1_total = clock_t::now();

    // Report timings
    const double sim_s   = std::chrono::duration<double>(t1_sim   - t0_sim  ).count();
    const double total_s = std::chrono::duration<double>(t1_total - t0_total).count();

    std::cout.setf(std::ios::fixed);
    std::cout.precision(6);
    std::cout << "Simulation time: " << sim_s   << " s\n";
    std::cout << "Total time     : " << total_s << " s\n";

    return EXIT_SUCCESS;
}
catch (const std::exception& e) {
    std::cerr << "Unhandled exception: " << e.what() << "\n";
    return EXIT_FAILURE;
}
catch (...) {
    std::cerr << "Unhandled non-standard exception.\n";
    return EXIT_FAILURE;
}