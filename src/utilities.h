#pragma once
#include <filesystem>
#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <cmath>
#include <atomic>

#include <Random123/philox.h>
#include <Random123/uniform.hpp>



// -------------------- R123Rng --------------------
class R123Rng {
public:
    typedef r123::Philox4x32 RNG;
    typedef RNG::ctr_type ctr_type;
    typedef RNG::key_type key_type;

private:
    RNG rng;
    ctr_type ctr;
    key_type key;
    ctr_type buffer;
    int buf_index = 4; // consume all initially (forces first refill)

public:
    R123Rng(); // default constructor
    R123Rng(uint64_t seed, uint64_t stream = 0);

    // Draw a uniform double in [0,1)
    double uniform();

    // Manually increment 128-bit counter
    void increment();

    // Reset/seed the RNG
    void reseed(uint64_t seed, uint64_t stream = 0);

    // Fork a new RNG with modified stream ID (for particles)
    R123Rng fork(uint64_t stream_offset) const;
};


// ---------- classes ----------
class Material {
public:
    int mat_id = -1;
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


struct Direction {
    double ux = 0.0;
    double uy = 0.0;
    double uz = 1.0; 

    inline void normalize() {
        double n = std::sqrt(ux * ux + uy * uy + uz * uz);
        if (n == 0.0) {
            ux = 0.0;
            uy = 0.0;
            uz = 1.0;
            return;
        }
        double inv = 1.0 / n;
        ux *= inv;
        uy *= inv;
        uz *= inv;
    }

    static Direction from_mu_phi(double mu, double phi) {
        double sinTheta = std::sqrt(std::max(0.0, 1.0 - mu * mu));
        Direction d;
        d.ux = sinTheta * std::cos(phi);
        d.uy = sinTheta * std::sin(phi);
        d.uz = mu;
        return d;
    }

    // Rotate current direction by polar angle (cosTheta) and azimuth (phi)
    // using a local orthonormal basis {u_hat, a_hat, b_hat}.
    inline void rotate(double cosTheta, double phi) {
        double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
        double cphi = std::cos(phi);
        double sphi = std::sin(phi);

        double ux0 = ux, uy0 = uy, uz0 = uz;

        // Build a_hat perpendicular to u_hat
        double ax, ay, az;
        if (std::fabs(uz0) < 0.999999) {
            // cross(u, ẑ) normalized
            ax = uy0;
            ay = -ux0;
            az = 0.0;
            double an = std::sqrt(ax * ax + ay * ay);
            ax /= an;
            ay /= an;
        } else {
            // u ≈ ±ẑ → pick x̂ to avoid singularity
            ax = 0.0;
            ay = (uz0 > 0.0) ? 1.0 : -1.0;
            az = 0.0;
        }

        // b_hat = cross(u_hat, a_hat)
        double bx = uy0 * az - uz0 * ay;
        double by = uz0 * ax - ux0 * az;
        double bz = ux0 * ay - uy0 * ax;

        // New direction
        double nx = cosTheta * ux0 + sinTheta * (cphi * ax + sphi * bx);
        double ny = cosTheta * uy0 + sinTheta * (cphi * ay + sphi * by);
        double nz = cosTheta * uz0 + sinTheta * (cphi * az + sphi * bz);

        ux = nx;
        uy = ny;
        uz = nz;

        // Normalize occasionally to mitigate round-off drift
        normalize();
    }
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

    //direction
    Direction dir;

    double t = 0.0;  // shk

    double energy = 0.0; // MeV
    double speed = 0.0; // cm/shk

    R123Rng rng;

    int id = -1;

    Particle() = default;

    // Easy initializer constructor
    Particle(const std::string &s,
             double x_, double y_, double z_,
             double t_,
             double energy_,
             double weight_ = 1.0);

    void start_particle_rng();
};

// -------------------- Generic N-D array wrapper --------------------
struct NDArray {
    std::vector<size_t> dims;     // e.g., {S, T, E}
    std::vector<size_t> strides;  // row-major
    std::vector<double> data;     // flat storage

    static std::vector<size_t> make_strides(const std::vector<size_t>& d);
    void resize(const std::vector<size_t>& d, double init_value = 0.0);
    size_t flat_index(const std::vector<size_t>& idx) const;
    double& at(const std::vector<size_t>& idx);
    const double& at(const std::vector<size_t>& idx) const;

    size_t size() const { return data.size(); }
};

// -------------------- Generic tally dimension descriptions --------

// A numeric binned dimension: e.g. "time", "energy", "x", "y".
struct TallyDim {
    std::string       name;   // label, e.g. "time"
    std::vector<double> edges; // bin edges (>= 2)
};

// A single coordinate along a named dimension: {"energy", 5.0}
struct DimCoord {
    std::string name;
    double      value;
};

// A straight-line segment along a named dimension: {"energy", e0, e1}
struct DimRange {
    std::string name;
    double      v0;
    double      v1;
};

// -------------------- Tally --------------------
class Tally {
public:
    std::string              tally_name;     // e.g., "test_tally_1"
    std::vector<std::string> species;        // species labels (first dimension)
    std::string              tally_category; // free-form label

    std::vector<bool> problem_defined = {false, false, false, false}; // time, x, y, z

    // Generic numeric binned dimensions (time, energy, x, y, ...)
    std::vector<TallyDim> dims;              // order = NDArray dimension order (after species)

    // Built on finalize()
    std::unordered_map<std::string, int> species_index; // label -> species index
    std::unordered_map<std::string, int> dim_index;     // dim name -> dim slot (0..dims.size()-1)
    NDArray counts;                                     // dims = {S, N0, N1, ...}

    // Build indices and allocate counts
    void finalize();

    // Return the list of dimension labels (excluding species).
    std::vector<std::string> dimension_labels() const;

    // Add a contribution at a point in all dimensions.
    // coords gives values for each named dimension (time, energy, x, ...).
    bool add(const std::string& sp,
             const std::vector<DimCoord>& coords,
             double contribution = 1.0);

    // Add a contribution smeared along a straight segment in ALL dimensions.
    // segments must specify v0,v1 for every numeric dimension in this tally.
    bool add_smear(const std::string& sp,
                   const std::vector<DimRange>& segments,
                   double contribution = 1.0);

    // Retrieve the value at a point (same interface as add, but read-only).
    double retrieve(const std::string& sp,
                    const std::vector<DimCoord>& coords) const;

    // Initialize everything at once, returns false if inputs are invalid.
    bool init(std::string name,
              std::vector<std::string> species_labels,
              std::vector<TallyDim>   dimensions,
              std::string category);

    // Factory that returns a ready-to-use Tally (throws on invalid input)
    static Tally Make(std::string name,
                      std::vector<std::string> species_labels,
                      std::vector<TallyDim>   dimensions,
                      std::string category);
};


// ------------------------ Geometry ------------------------------------------

class MeshCell {
public:
    int cell_id;

    Material cell_material;

    std::vector<double> surface_bounds; // [-x, +x, -y, +y, -z, +z] position of boundaries
    std::vector<double> boundary_conditions; 

    MeshCell() = default;
    explicit MeshCell(int id) : cell_id(id) {}

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

double particle_energy_2_speed(double energy, double mass);


// -------------------- Binning helper --------------------
// Given monotonically increasing bin *edges* (length N), return the bin index
// in [0..N-2] for value x. If x is exactly the last edge, return N-2.
// Returns -1 if x is out of range.
int bin_index_from_edges(const std::vector<double>& edges, double x);

// Optional alias if you prefer the singular name used in some of your notes.
inline int bin_index_from_edge(const std::vector<double>& edges, double x) {
    return bin_index_from_edges(edges, x);
}

