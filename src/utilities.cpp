#include "utilities.h"
#include "garage.h"

#include <iostream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <utility>
#include <algorithm>
#include <stdexcept>

// Draco includes
#include "units/PhysicalConstantsSI.hh"

// HDF5 C API is a C library; wrap in extern "C" for C++ builds
extern "C" {
#include <hdf5.h>
}



// -------------------- R123Rng Implementation --------------------

R123Rng::R123Rng() {
    reseed(0, 0);
}

R123Rng::R123Rng(uint64_t seed, uint64_t stream) {
    reseed(seed, stream);
}

void R123Rng::reseed(uint64_t seed, uint64_t stream) {
    key = key_type{{static_cast<uint32_t>(seed & 0xffffffffu),
                    static_cast<uint32_t>(seed >> 32)}};
    ctr = ctr_type{{static_cast<uint32_t>(stream & 0xffffffffu),
                    static_cast<uint32_t>(stream >> 32),
                    0u, 0u}};
    buf_index = 4;
}

void R123Rng::increment() {
    if (++ctr[0] == 0u) {
        if (++ctr[1] == 0u) {
            if (++ctr[2] == 0u) {
                ++ctr[3];
            }
        }
    }
}

double R123Rng::uniform() {
    if (buf_index >= 4) {
        buffer = rng(ctr, key);
        buf_index = 0;
        increment();
    }
    return r123::u01<double>(buffer[buf_index++]);
}

R123Rng R123Rng::fork(uint64_t stream_offset) const {
    // create a copy with different stream ID (e.g., particle ID)
    R123Rng new_rng;
    new_rng.key = this->key;
    new_rng.ctr = this->ctr;
    new_rng.ctr[2] = static_cast<uint32_t>(stream_offset & 0xffffffffu);
    new_rng.ctr[3] = static_cast<uint32_t>(stream_offset >> 32);
    new_rng.buf_index = 4;
    return new_rng;
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
    speed = std::sqrt(2.0 * (energy_*1e-3 * rtt_units::electronChargeSI * 1e6) / (projectile_mass * 1.0e-3)) * 1.0e-8 * 1.0e2; // cm/shk
    id = 1; // TODO: make this unique
}

void Particle::start_particle_rng() {

    // Derive a unique stream based on total particles created or ID
    uint64_t stream_id = static_cast<uint64_t>(garage.total_particles_created + id);
    rng = garage.rng.fork(stream_id);
}


void Tally::finalize() {
    // Build species lookup
    species_index.clear();
    for (int i = 0; i < static_cast<int>(species.size()); ++i) {
        species_index[species[i]] = i;
    }

    // Build dimension name -> index map and ensure each has valid edges
    dim_index.clear();
    for (int d = 0; d < static_cast<int>(dims.size()); ++d) {
        // Default numeric bin edges if not provided
        if (dims[d].edges.size() < 2) {
            dims[d].edges.clear();
            dims[d].edges.push_back(0.0);
            dims[d].edges.push_back(1e20);
        }
        dim_index[dims[d].name] = d;
    }

    // Build NDArray shape: {S, N0, N1, ...}
    std::vector<size_t> nd_dims;
    nd_dims.reserve(1 + dims.size());
    nd_dims.push_back(species.size());
    for (size_t d = 0; d < dims.size(); ++d) {
        nd_dims.push_back(dims[d].edges.size() - 1); // number of bins
    }

    counts.resize(nd_dims, 0.0);
}


std::vector<std::string> Tally::dimension_labels() const {
    std::vector<std::string> labels;
    labels.reserve(dims.size());
    for (size_t d = 0; d < dims.size(); ++d) {
        labels.push_back(dims[d].name);
    }
    return labels;
}

bool Tally::add(const std::string& sp,
                const std::vector<DimCoord>& coords,
                double contribution)
{
    // 1) Find species index
    std::unordered_map<std::string,int>::const_iterator it_s = species_index.find(sp);
    if (it_s == species_index.end()) return false;

    // 2) Build a quick lookup from dim name -> coordinate value
    //    (coords.size() is usually small, so a linear search would also be fine.)
    // For simplicity, use linear search.
    std::vector<size_t> idx;
    idx.reserve(1 + dims.size());
    idx.push_back(static_cast<size_t>(it_s->second)); // species index is first

    // 3) For each numeric dimension, find the coordinate and bin it
    for (size_t d = 0; d < dims.size(); ++d) {
        const std::string& dim_name = dims[d].name;

        bool   found = false;
        double value = 0.0;
        for (size_t i = 0; i < coords.size(); ++i) {
            if (coords[i].name == dim_name) {
                value = coords[i].value;
                found = true;
                break;
            }
        }
        if (!found) return false; // missing coordinate for this dimension

        int b = bin_index_from_edges(dims[d].edges, value);
        if (b < 0) return false; // out of range

        idx.push_back(static_cast<size_t>(b));
    }

    // 4) Accumulate
    counts.at(idx) += contribution;
    return true;
}

// -----------------------------------------------------------------------------
// Helper: add all param breakpoints s in (0,1) where a0 + (a1-a0)*s == edge
// -----------------------------------------------------------------------------
inline void push_breakpoints(double a0,
                             double a1,
                             const std::vector<double>& edges,
                             std::vector<double>& svec)
{
    const double da = a1 - a0;
    if (std::abs(da) < 1e-300) return; // no variation along this axis this step

    for (double edge : edges) {
        const double s = (edge - a0) / da;
        // Accept only strict interior crossings; endpoints are already in {0,1}
        if (s > 0.0 && s < 1.0)
            svec.push_back(s);
    }
}

// -----------------------------------------------------------------------------
// Helper: sort and unique with tolerance
// -----------------------------------------------------------------------------
inline void sort_unique_with_tol(std::vector<double>& v, double tol = 1e-12)
{
    std::sort(v.begin(), v.end());
    std::vector<double> u;
    u.reserve(v.size());
    for (double x : v) {
        if (u.empty() || std::abs(x - u.back()) > tol)
            u.push_back(x);
    }
    v.swap(u);
}

bool Tally::add_smear(const std::string& sp,
                      const std::vector<DimRange>& segments,
                      double contribution)
{
    // 1) Validate species
    std::unordered_map<std::string,int>::const_iterator it_s = species_index.find(sp);
    if (it_s == species_index.end()) return false;
    const size_t s_idx = static_cast<size_t>(it_s->second);

    const size_t D = dims.size();

    // Map each dimension index to its (v0,v1)
    std::vector<double> a0(D, 0.0);
    std::vector<double> a1(D, 0.0);
    std::vector<bool>   have(D, false);

    // Consume segments: ignore any whose name is NOT a dimension in this tally.
    for (size_t i = 0; i < segments.size(); ++i) {
        const DimRange& seg = segments[i];

        std::unordered_map<std::string,int>::const_iterator it_d = dim_index.find(seg.name);
        if (it_d == dim_index.end()) {
            // Dimension name not present in this tally; ignore it.
            continue;
        }

        const size_t d = static_cast<size_t>(it_d->second);
        if (have[d]) {
            // Duplicate specification for the same dimension -> ambiguous, treat as error.
            return false;
        }

        a0[d]   = seg.v0;
        a1[d]   = seg.v1;
        have[d] = true;
    }

    // Ensure all tally dimensions were specified at least once
    for (size_t d = 0; d < D; ++d) {
        if (!have[d]) return false;
    }

    // 2) Quick reject if entire segment is outside all bins in any dimension
    for (size_t d = 0; d < D; ++d) {
        const double vmin = std::min(a0[d], a1[d]);
        const double vmax = std::max(a0[d], a1[d]);

        const std::vector<double>& edges = dims[d].edges;
        if (vmax <= edges.front() || vmin >= edges.back()) {
            return false;
        }
    }

    // 3) Build list of parametric breakpoints s ∈ [0,1]
    std::vector<double> sbreaks;
    sbreaks.reserve(2 + D * 8); // rough guess
    sbreaks.push_back(0.0);
    sbreaks.push_back(1.0);

    for (size_t d = 0; d < D; ++d) {
        push_breakpoints(a0[d], a1[d], dims[d].edges, sbreaks);
    }
    sort_unique_with_tol(sbreaks);

    // 4) Walk sub-segments and add fractional contributions
    bool any = false;

    for (size_t i = 0; i + 1 < sbreaks.size(); ++i) {
        const double sL   = sbreaks[i];
        const double sR   = sbreaks[i + 1];
        const double frac = sR - sL;
        if (frac <= 0.0) continue;

        const double sm = 0.5 * (sL + sR);

        std::vector<size_t> idx;
        idx.reserve(1 + D);
        idx.push_back(s_idx);

        bool valid = true;
        for (size_t d = 0; d < D; ++d) {
            const double v  = a0[d] + (a1[d] - a0[d]) * sm;
            const int    bin = bin_index_from_edges(dims[d].edges, v);
            if (bin < 0) { valid = false; break; }
            idx.push_back(static_cast<size_t>(bin));
        }

        if (!valid) continue;

        counts.at(idx) += contribution * frac;
        any = true;
    }

    return any;
}

double Tally::retrieve(const std::string& sp,
                       const std::vector<DimCoord>& coords) const
{
    // 1) Find species index
    std::unordered_map<std::string,int>::const_iterator it_s = species_index.find(sp);
    if (it_s == species_index.end()) return 0.0;

    std::vector<size_t> idx;
    idx.reserve(1 + dims.size());
    idx.push_back(static_cast<size_t>(it_s->second));

    // 2) For each dimension, find the coordinate and bin index
    for (size_t d = 0; d < dims.size(); ++d) {
        const std::string& dim_name = dims[d].name;

        bool   found = false;
        double value = 0.0;
        for (size_t i = 0; i < coords.size(); ++i) {
            if (coords[i].name == dim_name) {
                value = coords[i].value;
                found = true;
                break;
            }
        }
        if (!found) return 0.0; // missing coordinate

        const int b = bin_index_from_edges(dims[d].edges, value);
        if (b < 0) return 0.0; // out of range

        idx.push_back(static_cast<size_t>(b));
    }

    // 3) Access NDArray
    return counts.at(idx);
}

bool Tally::init(std::string name,
                 std::vector<std::string> species_labels,
                 std::vector<TallyDim>   dimensions,
                 std::string category)
{
    if (species_labels.empty()) return false;

    tally_name     = std::move(name);
    species        = std::move(species_labels);
    dims           = std::move(dimensions);
    tally_category = std::move(category);

    finalize();
    return true;
}

Tally Tally::Make(std::string name,
                  std::vector<std::string> species_labels,
                  std::vector<TallyDim>   dimensions,
                  std::string category)
{
    if (species_labels.empty()) {
        throw std::invalid_argument("Tally::Make: species list must not be empty.");
    }

    Tally t;
    t.tally_name     = std::move(name);
    t.species        = std::move(species_labels);
    t.dims           = std::move(dimensions);
    t.tally_category = std::move(category);
    t.finalize();
    return t;
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

int species_2_z(std::string ion_species){
    int ion_z;
    if (ion_species == "d") {
        ion_z = 1;
    } else if (ion_species == "t") {
        ion_z = 1;
    } else if (ion_species == "a") {
        ion_z = 2;
    }
    return ion_z;
}

double species_2_mass(std::string ion_species){
    double ion_mass;
    if (ion_species == "d") {
        ion_mass = 3.34358e-24; // [g/atom]
    } else if (ion_species == "t") {
        ion_mass = 5.00736e-24;
    } else if (ion_species == "a") {
        ion_mass = 6.64465e-24;
    } else if (ion_species == "e") {
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


double particle_energy_2_speed(double energy, double mass) {
    double speed = std::sqrt(2.0 * (energy*1e-3 * rtt_units::electronChargeSI * 1e6) / (mass * 1.0e-3)) * 1.0e-8 * 1.0e2; // cm/shk

    return speed;
}


// -------------------- NDArray impl --------------------
std::vector<size_t> NDArray::make_strides(const std::vector<size_t>& d) {
    std::vector<size_t> s(d.size(), 1);
    if (d.empty()) return s;
    for (int i = static_cast<int>(d.size()) - 2; i >= 0; --i) {
        s[i] = s[i + 1] * d[i + 1];
    }
    return s;
}

void NDArray::resize(const std::vector<size_t>& d, double init_value) {
    dims = d;
    strides = make_strides(dims);
    size_t total = 1;
    for (auto v : dims) total *= v;
    data.assign(total, init_value);
}

size_t NDArray::flat_index(const std::vector<size_t>& idx) const {
    size_t off = 0;
    for (size_t i = 0; i < idx.size(); ++i) {
        off += idx[i] * strides[i];
    }
    return off;
}

double& NDArray::at(const std::vector<size_t>& idx) {
    return data[flat_index(idx)];
}

const double& NDArray::at(const std::vector<size_t>& idx) const {
    return data[flat_index(idx)];
}

// -------------------- Binning helper impl --------------------
int bin_index_from_edges(const std::vector<double>& edges, double x) {
    const size_t N = edges.size();
    if (N < 2) return -1;

    // Require x within [edges.front(), edges.back()]
    if (x < edges.front()) return -1;
    if (x > edges.back())  return -1;

    // Include exact upper edge in the last bin
    if (x == edges.back()) return static_cast<int>(N) - 2;

    // First edge strictly greater than x, step back one
    auto it = std::upper_bound(edges.begin(), edges.end(), x);
    if (it == edges.begin()) return -1; // shouldn't happen due to earlier check
    return static_cast<int>(std::distance(edges.begin(), it)) - 1;
}