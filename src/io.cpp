// io.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>

#include "utilities.h"
#include "garage.h"
#include "io.h"

// HDF5 C API is a C library; wrap in extern "C" for C++ builds
extern "C" {
#include <hdf5.h>
}


// ----  (file-scope) ----
struct PaintBlock {
    int ix0, ix1, iy0, iy1, iz0, iz1;
    int mat_id;
};
static bool g_have_default_fill = false;
static int  g_default_fill = -1;
static std::vector<PaintBlock> g_paints;


inline bool mesh_ready(const Garage& g) {
    return g.nx>0 && g.ny>0 && g.nz>0 && g.mesh_cells.size()==static_cast<size_t>(g.nx*g.ny*g.nz);
}

inline size_t mesh_flatten(const Garage& g, int ix, int iy, int iz) {
    // row-major with x fastest
    return static_cast<size_t>(ix)
         + static_cast<size_t>(iy) * g.mesh_strides[1]
         + static_cast<size_t>(iz) * g.mesh_strides[2];
}

inline void mesh_unflatten(const Garage& g, size_t id, int &ix, int &iy, int &iz) {
    iz = static_cast<int>(id / g.mesh_strides[2]);
    size_t r = id % g.mesh_strides[2];
    iy = static_cast<int>(r / g.mesh_strides[1]);
    ix = static_cast<int>(r % g.mesh_strides[1]);
}

inline bool mesh_in_bounds(const Garage& g, int ix, int iy, int iz) {
    return (ix>=0 && ix<g.nx) && (iy>=0 && iy<g.ny) && (iz>=0 && iz<g.nz);
}

// neighbor: axis a∈{0:x,1:y,2:z}, step s∈{-1,+1}; returns -1 if off-domain
inline int mesh_neighbor_id(const Garage& g, int cell_id, int a, int s) {
    int ix,iy,iz;
    mesh_unflatten(g, static_cast<size_t>(cell_id), ix, iy, iz);
    if (a==0) ix += s; else if (a==1) iy += s; else iz += s;
    if (!mesh_in_bounds(g, ix, iy, iz)) return -1;
    return static_cast<int>(mesh_flatten(g, ix, iy, iz));
}

// ------------------------ Parsing helpers (file-local) ----------------------

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

static std::vector<double> make_bins(double lo, double hi, int divs, bool logspace) {
    if (divs < 0) throw std::runtime_error("dividers must be >= 0");
    if (hi <= lo) throw std::runtime_error("upper bound must be > lower bound");
    const int n = divs + 2; // edges = dividers + 2 (inclusive endpoints)

    std::vector<double> edges;
    edges.reserve(n);

    if (!logspace) {
        const double step = (hi - lo) / (n - 1);
        for (int i = 0; i < n; ++i) edges.push_back(lo + step * i);
        return edges;
    }

    if (lo <= 0.0 || hi <= 0.0)
        throw std::runtime_error("log bins require lo>0 and hi>0");

    const double ratio = std::pow(hi / lo, 1.0 / (n - 1));
    double v = lo;
    for (int i = 0; i < n; ++i) {
        edges.push_back(v);
        v *= ratio;
    }
    edges.back() = hi; // avoid last-step floating error
    return edges;
}

static bool parse_bins_spec(const std::vector<std::string>& tok,
                            std::vector<double>& out_edges) {
    // expected: <key> <lo> <hi> <divs> <lin|log>
    if (tok.size() != 5) return false;
    double lo = std::stod(tok[1]);
    double hi = std::stod(tok[2]);
    int    dv = std::stoi(tok[3]);

    std::string mode = tok[4];
    for (auto &c : mode) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    bool logspace = (mode == "log" || mode == "log10" || mode == "logspace");

    out_edges = make_bins(lo, hi, dv, logspace);
    return true;
}

static bool is_number(const std::string &s) {
    char* end=nullptr;
    std::strtod(s.c_str(), &end);
    return end && *end=='\0';
}

// Accept either:
//   x <lo> <hi> <divs> <lin|log>
// or
//   x v0 v1 v2 ... vN   (N>=1; strictly increasing)
static bool parse_axis_edges(const std::vector<std::string>& tok, std::vector<double>& out_edges) {
    if (tok.size() < 3) return false;

    // Heuristic: if there are exactly 5 tokens and the 5th is lin/log, treat as bins-spec.
    if (tok.size() == 5) {
        std::string last = tok[4];
        for (auto &c : last) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (last=="lin" || last=="log" || last=="log10" || last=="logspace") {
            return parse_bins_spec(tok, out_edges);
        }
    }

    // Otherwise treat as explicit edge list: x v0 v1 ...
    out_edges.clear();
    out_edges.reserve(tok.size()-1);
    for (size_t i=1;i<tok.size();++i) {
        if (!is_number(tok[i])) return false;
        out_edges.push_back(std::stod(tok[i]));
    }
    if (out_edges.size() < 2) return false;
    // verify strictly increasing
    for (size_t i=1;i<out_edges.size();++i) {
        if (!(out_edges[i] > out_edges[i-1])) {
            std::cerr << "ERROR: geometry axis edges must be strictly increasing.\n";
            return false;
        }
    }
    return true;
}

bool build_cartesian_mesh_after_parse() {
    // Validate edges
    if (garage.x_bin_bounds.size()<2 || garage.y_bin_bounds.size()<2 || garage.z_bin_bounds.size()<2) {
        std::cerr << "ERROR: geometry requires x,y,z with at least two edges each.\n";
        return false;
    }
    auto check_increasing = [](const std::vector<double>& e, const char* name)->bool{
        for (size_t i=1;i<e.size();++i) if (!(e[i]>e[i-1])) {
            std::cerr << "ERROR: " << name << " edges must be strictly increasing.\n";
            return false;
        }
        return true;
    };
    if (!check_increasing(garage.x_bin_bounds,"x") ||
        !check_increasing(garage.y_bin_bounds,"y") ||
        !check_increasing(garage.z_bin_bounds,"z")) return false;

    garage.nx = static_cast<int>(garage.x_bin_bounds.size()) - 1;
    garage.ny = static_cast<int>(garage.y_bin_bounds.size()) - 1;
    garage.nz = static_cast<int>(garage.z_bin_bounds.size()) - 1;

    const size_t ncell = static_cast<size_t>(garage.nx)*garage.ny*garage.nz;
    garage.mesh_strides[0] = 1;
    garage.mesh_strides[1] = static_cast<size_t>(garage.nx);
    garage.mesh_strides[2] = static_cast<size_t>(garage.nx)*garage.ny;

    for (int iz=0; iz<garage.nz; ++iz) {
      for (int iy=0; iy<garage.ny; ++iy) {
        for (int ix=0; ix<garage.nx; ++ix) {
            const size_t id = mesh_flatten(garage, ix, iy, iz);
            MeshCell cell(static_cast<int>(id));
            cell.surface_bounds.resize(6);
            cell.boundary_conditions.resize(6, -1);

            cell.surface_bounds[0] = garage.x_bin_bounds[ix];
            cell.surface_bounds[1] = garage.x_bin_bounds[ix+1];
            cell.surface_bounds[2] = garage.y_bin_bounds[iy];
            cell.surface_bounds[3] = garage.y_bin_bounds[iy+1];
            cell.surface_bounds[4] = garage.z_bin_bounds[iz];
            cell.surface_bounds[5] = garage.z_bin_bounds[iz+1];

            // -x, +x
            cell.boundary_conditions[0] = mesh_in_bounds(garage, ix-1,iy,iz) ? static_cast<int>(mesh_flatten(garage, ix-1,iy,iz)) : -1;
            cell.boundary_conditions[1] = mesh_in_bounds(garage, ix+1,iy,iz) ? static_cast<int>(mesh_flatten(garage, ix+1,iy,iz)) : -1;
            // -y, +y
            cell.boundary_conditions[2] = mesh_in_bounds(garage, ix,iy-1,iz) ? static_cast<int>(mesh_flatten(garage, ix,iy-1,iz)) : -1;
            cell.boundary_conditions[3] = mesh_in_bounds(garage, ix,iy+1,iz) ? static_cast<int>(mesh_flatten(garage, ix,iy+1,iz)) : -1;
            // -z, +z
            cell.boundary_conditions[4] = mesh_in_bounds(garage, ix,iy,iz-1) ? static_cast<int>(mesh_flatten(garage, ix,iy,iz-1)) : -1;
            cell.boundary_conditions[5] = mesh_in_bounds(garage, ix,iy,iz+1) ? static_cast<int>(mesh_flatten(garage, ix,iy,iz+1)) : -1;

            garage.mesh_cells[id] = std::move(cell);
        }
      }
    }

    // painting loops: same change
    for (const auto &pb : g_paints) {
        int ix0 = std::clamp(pb.ix0, 0, garage.nx);
        int ix1 = std::clamp(pb.ix1, 0, garage.nx);
        int iy0 = std::clamp(pb.iy0, 0, garage.ny);
        int iy1 = std::clamp(pb.iy1, 0, garage.ny);
        int iz0 = std::clamp(pb.iz0, 0, garage.nz);
        int iz1 = std::clamp(pb.iz1, 0, garage.nz);
        for (int iz=iz0; iz<iz1; ++iz)
          for (int iy=iy0; iy<iy1; ++iy)
            for (int ix=ix0; ix<ix1; ++ix) {
                size_t id = mesh_flatten(garage, ix, iy, iz);
                garage.mesh_mat_ids[id] = pb.mat_id;
            }
    }


    // Optionally resolve cell_material for each MeshCell from mat_id table
    // If you prefer to keep only mat_id per cell, you can skip this resolve and
    // look up materials by mat_id at usage sites.
    auto find_material_by_id = [&](int mid)->const Material*{
        for (const auto &m : garage.materials) {
            // you stored mat_id as a double; cast to int safely
            int im = static_cast<int>(std::llround(m.mat_id));
            if (im == mid) return &m;
        }
        return nullptr;
    };
    for (size_t id=0; id<ncell; ++id) {
        int mid = garage.mesh_mat_ids[id];
        if (mid >= 0) {
            if (const Material* mp = find_material_by_id(mid)) {
                garage.mesh_cells[id].cell_material = *mp; // copy now; or store pointer if you prefer
            } else {
                std::cerr << "WARN: mat_id " << mid << " not found in [materials]; cell " << id << " left default.\n";
            }
        }
    }

    // clear paint staging for next parse
    g_have_default_fill = false; g_default_fill = -1; g_paints.clear();

    return true;
}

// ----------------------------- Input parsing --------------------------------

bool parse_input_file(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "ERROR: Cannot open input file: " << path << "\n";
        return false;
    }

    enum class Section { None, Materials, Source, Tallies, Geometry, Settings };
    Section sec = Section::None;

    Material mat;         bool mat_active   = false;
    Tally    tally;       bool tally_active = false;
    std::vector<TallyDim> tally_dims;

    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;

        // section headers
        if (line.front() == '[' && line.back() == ']') {
            std::string name = line.substr(1, line.size()-2);

            // Finalize currently active material if leaving [materials]
            if (sec == Section::Materials && mat_active) {
                if (mat.species.size() != mat.densities.size()) {
                    std::cerr << "ERROR: materials: species count (" << mat.species.size()
                              << ") != densities count (" << mat.densities.size() << ")\n";
                    return false;
                }
                garage.materials.push_back(std::move(mat));
                mat = Material{};
                mat_active = false;
            }

            // Finalize currently active tally if leaving [tallies]
            if (sec == Section::Tallies && tally_active) {
                if (tally.tally_name.empty())
                    std::cerr << "WARN: tally without name; will still push\n";
                tally.dims = std::move(tally_dims);
                tally.finalize();
                garage.tallies.push_back(std::move(tally));
                tally = Tally{};
                tally_dims.clear();
                tally_active = false;
            }

            if      (name == "materials") sec = Section::Materials;
            else if (name == "source")    sec = Section::Source;
            else if (name == "tallies")   sec = Section::Tallies;
            else if (name == "geometry")  sec = Section::Geometry;
            else if (name == "settings")  sec = Section::Settings;
            else                          sec = Section::None;
            continue;
        }

        // content lines
        auto tok = split_ws(line);
        if (tok.empty()) continue;

        try {
            switch (sec) {
            case Section::Materials: {
                if (tok[0] == "ion_temperature" && tok.size() == 2) {
                    mat.ion_temperature = std::stod(tok[1]);
                    mat_active = true;
                } else if (tok[0] == "electron_temperature" && tok.size() == 2) {
                    mat.electron_temperature = std::stod(tok[1]);
                    mat_active = true;
                } else if (tok[0] == "particles" && tok.size() >= 2) {
                    mat.species.assign(tok.begin() + 1, tok.end());
                    mat_active = true;
                } else if (tok[0] == "densities" && tok.size() >= 2) {
                    mat.densities.clear();
                    mat.densities.reserve(tok.size() - 1);
                    for (size_t i = 1; i < tok.size(); ++i)
                        mat.densities.push_back(std::stod(tok[i]));
                    mat_active = true;
                } else if (tok[0] == "mat_id" && tok.size() == 2) {
                    mat.mat_id = std::stod(tok[1]);
                } else {
                    std::cerr << "WARN: Unrecognized materials line: " << line << "\n";
                }
            } break;

            case Section::Source: {
                if (tok[0] == "particle" && tok.size() == 2) {
                    garage.source_particle = tok[1];
                } else if (tok[0] == "point" && tok.size() == 4) {
                    garage.source_point[0] = std::stod(tok[1]);
                    garage.source_point[1] = std::stod(tok[2]);
                    garage.source_point[2] = std::stod(tok[3]);
                } else if (tok[0] == "time" && tok.size() == 2) {
                    garage.source_time = std::stod(tok[1]);
                } else if (tok[0] == "energy" && tok.size() == 2) {
                    garage.source_energy = std::stod(tok[1]);
                } else if (tok[0] == "strength" && tok.size() == 2) {
                    garage.source_strength = std::stod(tok[1]);
                } else if (tok[0] == "direction" && tok.size() == 4) {
                    garage.source_direction[0] = std::stod(tok[1]);
                    garage.source_direction[1] = std::stod(tok[2]);
                    garage.source_direction[2] = std::stod(tok[3]);
                } else {
                    std::cerr << "WARN: Unrecognized source line: " << line << "\n";
                }
            } break;

            case Section::Tallies: {
                if (tok[0] == "name" && tok.size() == 2) {
                    tally.tally_name = tok[1];
                    tally_active = true;
                } else if (tok[0] == "particles" && tok.size() >= 2) {
                    tally.species.assign(tok.begin() + 1, tok.end());
                    tally_active = true;
                } else if (tok[0] == "category" && tok.size() == 2) {
                    tally.tally_category = tok[1];
                    tally_active = true;
                } else if (tok[0] == "energy_bins") {
                    TallyDim energy_dim;
                    energy_dim.name = "energy";
                    if (!parse_bins_spec(tok, energy_dim.edges)) {
                        std::cerr << "ERROR: bad energy_bins spec: " << line << "\n";
                        return false;
                    }
                    tally_active = true;
                    tally_dims.push_back(energy_dim);
                } else if (tok[0] == "time_bins") {
                    TallyDim time_dim;
                    time_dim.name = "time";
                    if (tok.size() == 2 && tok[1] == "timesteps") {
                        time_dim.edges = {-1e20, -1e20};
                        tally.problem_defined[0] = true;
                    } else if (!parse_bins_spec(tok, time_dim.edges)) {
                        std::cerr << "ERROR: bad time_bins spec: " << line << "\n";
                        return false;
                    }
                    tally_active = true;
                    tally_dims.push_back(time_dim);
                } else if (tok[0] == "x_bins") {
                    TallyDim x_dim;
                    x_dim.name = "x";
                    if (tok.size() == 2 && tok[1] == "cells") {
                        x_dim.edges = {-1e20, -1e20};
                        tally.problem_defined[1] = true;
                    } else if (!parse_bins_spec(tok, x_dim.edges)) {
                        std::cerr << "ERROR: bad x_bins spec: " << line << "\n";
                        return false;
                    }
                    tally_dims.push_back(x_dim);
                } else if (tok[0] == "y_bins") {
                    TallyDim y_dim;
                    y_dim.name = "y";
                    if (tok.size() == 2 && tok[1] == "cells") {
                        y_dim.edges = {-1e20, -1e20};
                        tally.problem_defined[2] = true;
                    } else if (!parse_bins_spec(tok, y_dim.edges)) {
                        std::cerr << "ERROR: bad y_bins spec: " << line << "\n";
                        return false;
                    }
                    tally_dims.push_back(y_dim);
                } else if (tok[0] == "z_bins") {
                    TallyDim z_dim;
                    z_dim.name = "z";
                    if (tok.size() == 2 && tok[1] == "cells") {
                        z_dim.edges = {-1e20, -1e20};
                        tally.problem_defined[3] = true;
                    } else if (!parse_bins_spec(tok, z_dim.edges)) {
                        std::cerr << "ERROR: bad z_bins spec: " << line << "\n";
                        return false;
                    }
                    tally_dims.push_back(z_dim);
                } else {
                    std::cerr << "WARN: Unrecognized tallies line: " << line << "\n";
                }
            } break;

            case Section::Settings: {
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

            case Section::Geometry: {
                if (tok[0] == "x") {
                    if (!parse_axis_edges(tok, garage.x_bin_bounds)) { std::cerr << "ERROR: bad x edges spec: " << line << "\n"; return false; }
                } else if (tok[0] == "y") {
                    if (!parse_axis_edges(tok, garage.y_bin_bounds)) { std::cerr << "ERROR: bad y edges spec: " << line << "\n"; return false; }
                } else if (tok[0] == "z") {
                    if (!parse_axis_edges(tok, garage.z_bin_bounds)) { std::cerr << "ERROR: bad z edges spec: " << line << "\n"; return false; }
                } else if (tok[0] == "fill" && tok.size()==2) {
                    g_have_default_fill = true;
                    g_default_fill = std::stoi(tok[1]);
                } else if (tok[0] == "block" && tok.size()==8) {
                    PaintBlock pb;
                    pb.ix0 = std::stoi(tok[1]); pb.ix1 = std::stoi(tok[2]);
                    pb.iy0 = std::stoi(tok[3]); pb.iy1 = std::stoi(tok[4]);
                    pb.iz0 = std::stoi(tok[5]); pb.iz1 = std::stoi(tok[6]);
                    pb.mat_id = std::stoi(tok[7]);
                    if (!(pb.ix0<=pb.ix1 && pb.iy0<=pb.iy1 && pb.iz0<=pb.iz1)) {
                        std::cerr << "ERROR: block ranges must be non-decreasing.\n";
                        return false;
                    }
                    g_paints.push_back(pb);
                } else {
                    std::cerr << "WARN: Unrecognized geometry line: " << line << "\n";
                }
            } break;

            case Section::None:
            default:
                std::cerr << "WARN: Content outside recognized section: " << line << "\n";
                break;
            }
        } catch (const std::exception& e) {
            std::cerr << "ERROR: parse error on line: \"" << line << "\" : " << e.what() << "\n";
            return false;
        }
    }

    // Finalize any active material at EOF
    if (!mat.species.empty() || !mat.densities.empty()
        || mat.ion_temperature != 0.0 || mat.electron_temperature != 0.0)
    {
        if (mat.species.size() != mat.densities.size()) {
            std::cerr << "ERROR: materials: species count (" << mat.species.size()
                      << ") != densities count (" << mat.densities.size() << ")\n";
            return false;
        }
        garage.materials.push_back(std::move(mat));
    }

    // Finalize any active tally at EOF
    if (!tally.tally_name.empty() || !tally.species.empty())
    {
        tally.dims = std::move(tally_dims);
        tally.finalize();
        garage.tallies.push_back(std::move(tally));
    }

    return true;
}

// Reserve capacity for particle banks based on num_particles.
void plan_particle_capacity_after_parse() {
    const size_t n = static_cast<size_t>(std::max(garage.num_particles, 0));
    const size_t headroom = n + n / 2; // 1.5x

    garage.active_bank.reserve(n);
    garage.census_bank.reserve(headroom);
    garage.secondary_bank.reserve(headroom);
}

// ------------------------------ HDF5 helpers --------------------------------

static bool h5_ok(herr_t st) { return st >= 0; }

static bool h5_mkgroup(hid_t parent, const char* name) {
    hid_t g = H5Gcreate2(parent, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (g < 0) return false;
    return h5_ok(H5Gclose(g));
}

static hid_t h5_open_group(hid_t parent, const char* name) {
    return H5Gopen2(parent, name, H5P_DEFAULT);
}

static hid_t h5_create_group(hid_t parent, const char* name) {
    return H5Gcreate2(parent, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

static bool h5_write_scalar_double(hid_t parent, const char* name, double v) {
    hsize_t dims[1] = {1};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    if (space < 0) return false;
    hid_t dset = H5Dcreate2(parent, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    bool ok = (dset >= 0) && h5_ok(H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v));
    if (dset >= 0) H5Dclose(dset);
    H5Sclose(space);
    return ok;
}

static bool h5_write_scalar_int(hid_t parent, const char* name, int v) {
    hsize_t dims[1] = {1};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    if (space < 0) return false;
    hid_t dset = H5Dcreate2(parent, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    bool ok = (dset >= 0) && h5_ok(H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v));
    if (dset >= 0) H5Dclose(dset);
    H5Sclose(space);
    return ok;
}

static bool h5_write_string(hid_t parent, const char* name, const std::string& s) {
    hid_t type = H5Tcopy(H5T_C_S1);
    if (type < 0) return false;
    if (!h5_ok(H5Tset_size(type, H5T_VARIABLE))) { H5Tclose(type); return false; }

    hid_t space = H5Screate(H5S_SCALAR);
    if (space < 0) { H5Tclose(type); return false; }

    hid_t dset = H5Dcreate2(parent, name, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) { H5Sclose(space); H5Tclose(type); return false; }

    const char* ptr = s.c_str();
    bool ok = h5_ok(H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ptr));

    H5Dclose(dset);
    H5Sclose(space);
    H5Tclose(type);
    return ok;
}

static bool h5_write_strvec(hid_t parent, const char* name, const std::vector<std::string>& vec) {
    hsize_t dims[1] = { static_cast<hsize_t>(vec.size()) };
    hid_t space = H5Screate_simple(1, dims, nullptr);
    if (space < 0) return false;

    hid_t type = H5Tcopy(H5T_C_S1);
    if (type < 0) { H5Sclose(space); return false; }
    if (!h5_ok(H5Tset_size(type, H5T_VARIABLE))) { H5Tclose(type); H5Sclose(space); return false; }

    hid_t dset = H5Dcreate2(parent, name, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) { H5Tclose(type); H5Sclose(space); return false; }

    std::vector<const char*> cvec; cvec.reserve(vec.size());
    for (auto& s : vec) cvec.push_back(s.c_str());

    bool ok = h5_ok(H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, cvec.data()));

    H5Dclose(dset);
    H5Tclose(type);
    H5Sclose(space);
    return ok;
}

static bool h5_write_vec_double(hid_t parent, const char* name, const std::vector<double>& v) {
    hsize_t dims[1] = { static_cast<hsize_t>(v.size()) };
    hid_t space = H5Screate_simple(1, dims, nullptr);
    if (space < 0) return false;
    hid_t dset = H5Dcreate2(parent, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) { H5Sclose(space); return false; }
    bool ok = h5_ok(H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()));
    H5Dclose(dset);
    H5Sclose(space);
    return ok;
}

static bool h5_write_vec_int(hid_t parent, const char* name, const std::vector<int>& v) {
    hsize_t dims[1] = { static_cast<hsize_t>(v.size()) };
    hid_t space = H5Screate_simple(1, dims, nullptr);
    if (space < 0) return false;
    hid_t dset = H5Dcreate2(parent, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) { H5Sclose(space); return false; }
    bool ok = h5_ok(H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()));
    H5Dclose(dset);
    H5Sclose(space);
    return ok;
}

static bool h5_write_point3(hid_t parent, const char* name, const std::array<double,3>& p) {
    hsize_t dims[1] = { 3 };
    hid_t space = H5Screate_simple(1, dims, nullptr);
    if (space < 0) return false;
    hid_t dset = H5Dcreate2(parent, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) { H5Sclose(space); return false; }
    bool ok = h5_ok(H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, p.data()));
    H5Dclose(dset);
    H5Sclose(space);
    return ok;
}

static bool h5_write_ndarray_counts(hid_t parent, const char* name, const NDArray& arr) {
    std::vector<hsize_t> hd(arr.dims.begin(), arr.dims.end());
    hid_t space = H5Screate_simple(static_cast<int>(hd.size()), hd.data(), nullptr);
    if (space < 0) return false;

    hid_t dset = H5Dcreate2(parent, name, H5T_NATIVE_DOUBLE, space,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) { H5Sclose(space); return false; }

    bool ok = h5_ok(H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr.data.data()));

    H5Dclose(dset);
    H5Sclose(space);
    return ok;
}

// Write generic tally dimensions:
// - dim_names: ["time", "energy", "x", ...]
// - for each dim: "<name>_edges" dataset with its bin edges
static bool h5_write_tally_dims(hid_t parent, const std::vector<TallyDim>& dims_vec) {
    bool ok = true;

    // 1) Write the list of dimension names
    std::vector<std::string> names;
    names.reserve(dims_vec.size());
    for (size_t i = 0; i < dims_vec.size(); ++i) {
        names.push_back(dims_vec[i].name);
    }
    ok = ok && h5_write_strvec(parent, "dim_names", names);

    // 2) For each dimension, write its bin edges as "<name>_edges"
    for (size_t i = 0; i < dims_vec.size(); ++i) {
        const TallyDim& d = dims_vec[i];
        std::string dset_name = d.name + "_edges";
        ok = ok && h5_write_vec_double(parent, dset_name.c_str(), d.edges);
    }

    return ok;
}

static bool write_tally_block(hid_t parent_group, const std::vector<Tally>& tallies) {
    bool ok_local = true;

    for (const Tally& t : tallies) {
        hid_t g_one = h5_create_group(parent_group, t.tally_name.c_str());

        // Basic metadata
        ok_local = ok_local && h5_write_string    (g_one, "tally_name",     t.tally_name);
        ok_local = ok_local && h5_write_string    (g_one, "tally_category", t.tally_category);
        ok_local = ok_local && h5_write_strvec    (g_one, "species",        t.species);

        // NEW: generic dimensions
        ok_local = ok_local && h5_write_tally_dims(g_one, t.dims);

        // Species index map (optional, but you had it before; still valid)
        {
            std::vector<std::string> keys;
            std::vector<int>         vals;
            keys.reserve(t.species_index.size());
            vals.reserve(t.species_index.size());
            for (const auto& kv : t.species_index) {
                keys.push_back(kv.first);
                vals.push_back(kv.second);
            }
            ok_local = ok_local && h5_write_strvec(g_one, "species_index_keys",   keys);
            ok_local = ok_local && h5_write_vec_int(g_one, "species_index_values", vals);
        }

        // Counts array (same as before; NDArray now has shape {S, dim0_bins, dim1_bins, ...})
        ok_local = ok_local && h5_write_ndarray_counts(g_one, "counts", t.counts);

        H5Gclose(g_one);
        if (!ok_local) break;
    }

    return ok_local;
}

// ------------------------------ Output writer -------------------------------

bool write_output(const std::filesystem::path& out, double sim_s, double total_s) {
    hid_t file = H5Fcreate(out.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        std::cerr << "ERROR: Failed to create HDF5 file: " << out << "\n";
        return false;
    }

    bool ok = true;

    // Top-level groups
    ok = ok && h5_mkgroup(file, "/input");
    ok = ok && h5_mkgroup(file, "/tallies");
    ok = ok && h5_mkgroup(file, "/details");
    if (!ok) {
        std::cerr << "ERROR: Failed to create top-level groups.\n";
        H5Fclose(file);
        return false;
    }

    // ----- details -----
    {
        hid_t g = h5_open_group(file, "/details");
        ok = ok && h5_write_scalar_double(g, "simulation_time", sim_s);
        ok = ok && h5_write_scalar_double(g, "total_time",      total_s);
        H5Gclose(g);
    }

    // ----- input/settings and input/source -----
    {
        hid_t g_input = h5_open_group(file, "/input");

        // /input/settings
        hid_t g_sets = h5_create_group(g_input, "settings");
        ok = ok && h5_write_scalar_int(g_sets,    "num_particles", garage.num_particles);
        ok = ok && h5_write_scalar_int(g_sets,    "num_t_steps",   garage.num_t_steps);
        ok = ok && h5_write_scalar_double(g_sets, "t_step_size",   garage.t_step_size);
        H5Gclose(g_sets);

        // /input/source
        hid_t g_src = h5_create_group(g_input, "source");
        ok = ok && h5_write_string(g_src,       "particle", garage.source_particle);
        ok = ok && h5_write_point3(g_src,       "point",    garage.source_point);
        ok = ok && h5_write_scalar_double(g_src,"time",     garage.source_time);
        ok = ok && h5_write_scalar_double(g_src,"energy",   garage.source_energy);
        H5Gclose(g_src);

        H5Gclose(g_input);
    }

    // ----- tallies -----
    {
        hid_t g_tallies = h5_open_group(file, "/tallies");

        // Standard tallies
        hid_t g_std  = h5_create_group(g_tallies, "standard_tallies");
        ok = ok && write_tally_block(g_std, garage.standard_tallies);
        H5Gclose(g_std);

        // User tallies
        hid_t g_user = h5_create_group(g_tallies, "user_tallies");
        ok = ok && write_tally_block(g_user, garage.tallies);
        H5Gclose(g_user);

        H5Gclose(g_tallies);
    }

    // Close file
    if (!h5_ok(H5Fclose(file))) {
        std::cerr << "ERROR: Failed to close HDF5 file: " << out << "\n";
        return false;
    }
    if (!ok) {
        std::cerr << "ERROR: Failed writing one or more datasets/groups to: " << out << "\n";
        return false;
    }
    return true;
}