#pragma once
#include <filesystem>
#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>

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
    
    double weight = 1.0; //
    double energy_weight = 1.0; // keV

    double x = 0.0; // cm
    double y = 0.0;
    double z = 0.0;
    int zone = 0;

    double t = 0.0;  // shk

    double energy = 0.0; // MeV
    double speed = 0.0; // cm/shk

    int id = -1;

    Particle() = default;

    // Easy initializer constructor
    Particle(const std::string &s,
             double x_, double y_, double z_,
             double t_,
             double energy_,
             double weight_ = 1.0);
};

// -------------------- Generic N-D array wrapper --------------------
struct NDArray {
    std::vector<size_t> dims;     // e.g., {S, T, E}
    std::vector<size_t> strides;  // row-major
    std::vector<double> data;     // flat storage

    // Compute row-major strides from dims
    static std::vector<size_t> make_strides(const std::vector<size_t>& d);

    // Allocate/resize storage and set all values to init_value
    void resize(const std::vector<size_t>& d, double init_value = 0.0);

    // Convert a multi-index (same length as dims) to flat index
    size_t flat_index(const std::vector<size_t>& idx) const;

    // Element access by multi-index
    double& at(const std::vector<size_t>& idx);
    const double& at(const std::vector<size_t>& idx) const;

    size_t size() const { return data.size(); }
};

// -------------------- Tally --------------------
class Tally {
public:
    std::string              tally_name;     // e.g., "test_tally_1"
    std::vector<std::string> species;        // labels (order matters)
    std::vector<double>      energy_bins;    // bin edges (>=2)
    std::vector<double>      time_bins;      // bin edges (>=2)
    std::string              tally_category; // Options: [csd_energy_loss, ]

    // built on finalize()
    std::unordered_map<std::string, int> species_index; // label -> S index
    NDArray counts;                                      // dims = {S, T, E}

    // Ensure defaults, build species_index and counts (call after parsing)
    void finalize();

    // Add a contribution to the appropriate [species, time, energy] bin
    // returns false if out of range or species unknown
    bool add(const std::string& sp, double time, double energy, double contribution = 1.0);

    bool add_smear(const std::string& sp,
               double t0, double t1,
               double e0, double e1,
               double contribution = 1.0);

    // Retrieve the value in the counts array for a given species, time, and energy
    double retrieve(const std::string& sp, double time, double energy) const;

    // initialize everything at once, returns false if inputs are invalid
    bool init(std::string name,
              std::vector<std::string> species_labels,
              std::vector<double> energy_edges,
              std::vector<double> time_edges,
              std::string category);

    // factory that returns a ready-to-use Tally (throws on invalid input)
    static Tally Make(std::string name,
                      std::vector<std::string> species_labels,
                      std::vector<double> energy_edges,
                      std::vector<double> time_edges,
                      std::string category);

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

int species_2_z(std::string ion_species);

double species_2_mass(std::string ion_species);

double species_2_molar_mass(std::string ion_species);

double compute_mixture_molar_mass(const std::vector<std::string>& species,
                                  const std::vector<double>& densities);


// -------------------- Binning helper --------------------
// Given monotonically increasing bin *edges* (length N), return the bin index
// in [0..N-2] for value x. If x is exactly the last edge, return N-2.
// Returns -1 if x is out of range.
int bin_index_from_edges(const std::vector<double>& edges, double x);

// Optional alias if you prefer the singular name used in some of your notes.
inline int bin_index_from_edge(const std::vector<double>& edges, double x) {
    return bin_index_from_edges(edges, x);
}

