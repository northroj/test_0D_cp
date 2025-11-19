#include "loop.h"
#include "garage.h"
#include "utilities.h"

#include <algorithm>
#include <iterator>
#include <random>
#include <utility>
#include <iostream>
#include <numeric>
#include <cstdint>
#include <atomic>
#include <unordered_set>

#include <Random123/philox.h>
#include <Random123/uniform.hpp>

//Draco includes
#include "cdi_CPEloss/Analytic_CP_Eloss.hh"
#include "cdi_CPEloss/Analytic_Eloss_Model.hh"
#include "cdi_CPEloss/Analytic_KP_Alpha_Eloss_Model.hh"
#include "cdi_CPEloss/Analytic_Spitzer_Eloss_Model.hh"
#include "cdi/CPCommon.hh"
#include "units/PhysicalConstantsSI.hh"
#include "units/PhysicalConstexprs.hh"


// Move-append all of src into dst, then clear src.
static void move_append(std::vector<Particle>& dst, std::vector<Particle>& src) {
    dst.reserve(dst.size() + src.size());
    std::move(src.begin(), src.end(), std::back_inserter(dst));
    src.clear();
}

// Global-ish RNG for this TU (seeded once)
static std::mt19937& rng() {
    static std::mt19937 gen{std::random_device{}()};
    return gen;
}



void simulate() {

    // set up time steps
    time_step_setup("Uniform"); // TODO: add more type in the future

    // set up problem defined binning
    for (auto &t : garage.tallies) {
        if (t.problem_defined[0] == true) {
            for (auto &dim : t.dims) {
                if (dim.name == "time") {
                    dim.edges = garage.time_step_bins;
                }
            }
            t.finalize(); //redo this to make sure the tally dimensions are correct
        }
        if (t.problem_defined[1] == true) {
            for (auto &dim : t.dims) {
                if (dim.name == "x") {
                    dim.edges = garage.x_bin_bounds;
                }
            }
            t.finalize(); //redo this to make sure the tally dimensions are correct
        }
        if (t.problem_defined[2] == true) {
            for (auto &dim : t.dims) {
                if (dim.name == "y") {
                    dim.edges = garage.y_bin_bounds;
                }
            }
            t.finalize(); //redo this to make sure the tally dimensions are correct
        }
        if (t.problem_defined[3] == true) {
            for (auto &dim : t.dims) {
                if (dim.name == "z") {
                    dim.edges = garage.z_bin_bounds;
                }
            }
            t.finalize(); //redo this to make sure the tally dimensions are correct
        }
    }

    // set up problem defined dimensions for tallies
    TallyDim timestep_dim;
    timestep_dim.name = "time";
    timestep_dim.edges = garage.time_step_bins;
    TallyDim x_dim;
    TallyDim y_dim;
    TallyDim z_dim;
    x_dim.name = "x";
    y_dim.name = "y";
    z_dim.name = "z";
    x_dim.edges = garage.x_bin_bounds;
    y_dim.edges = garage.y_bin_bounds;
    z_dim.edges = garage.z_bin_bounds;
    // get a list of all species
    std::unordered_set<std::string> unique_species;
    for (auto &m : garage.materials) {
        for (auto &sp : m.species) {
            unique_species.insert(sp);
        }
    }
    std::vector<std::string> species_dim(unique_species.begin(), unique_species.end());
    // Initialize standard tallies
    Tally ion_temperature_time = Tally::Make("ion_temperature_time", {"all"}, {timestep_dim, x_dim, y_dim, z_dim}, "ion_temperature_time" );
    Tally electron_temperature_time = Tally::Make("electron_temperature_time", {"e"}, {timestep_dim, x_dim, y_dim, z_dim}, "electron_temperature_time" );
    Tally ion_density_time = Tally::Make("ion_density_time", species_dim, {timestep_dim, x_dim, y_dim, z_dim}, "ion_density_time"); // TODO: more than one material
    Tally energy_dump_ion = Tally::Make("energy_dump_ion", {"all"}, {timestep_dim, x_dim, y_dim, z_dim}, "energy_dump_ion");
    Tally energy_dump_electron = Tally::Make("energy_dump_electron", {"e"}, {timestep_dim, x_dim, y_dim, z_dim}, "energy_dump_electron");

    garage.standard_tallies.push_back(std::move(ion_temperature_time));
    garage.standard_tallies.push_back(std::move(electron_temperature_time));
    garage.standard_tallies.push_back(std::move(ion_density_time));
    garage.standard_tallies.push_back(std::move(energy_dump_ion));
    garage.standard_tallies.push_back(std::move(energy_dump_electron));

    // Loop over timesteps
    for(int t_it = 0; t_it < garage.num_t_steps; ++t_it) {
        garage.current_time_step = t_it;
        simulate_timestep(t_it);
    }
    
}

void simulate_timestep(int t_it) {

    // calculate the end of timestep time
    double time_census = garage.time_step_bins[t_it+1]; // shk
    // pre timestep time
    double time_start = garage.time_step_bins[t_it]; //shk

    source_particles(time_start, time_census);
    //population_control();
    
    // process active bank; when it empties, refill from secondary_bank
    auto& active    = garage.active_bank;
    auto& secondary = garage.secondary_bank;
    // while active + secondary banks have particles
    while (!active.empty() || !secondary.empty()) {
        // While we still have active particles, process them
        while (!active.empty()) {
            // Take from the back for O(1) remove
            Particle p = std::move(active.back());
            active.pop_back();

            transport_particle(p, time_census);
            // NOTE: transport_particle may push into census_bank or secondary_bank
        }

        // If active emptied but we have secondaries, move them in and continue
        if (!secondary.empty()) {
            move_append(active, secondary);
        }
    }

    // pick a value in the timestep
    double timestep_value = (garage.time_step_bins[t_it+1] - garage.time_step_bins[t_it])  /2.0  + garage.time_step_bins[t_it];

    // update material temperatures
    for (auto &cell : garage.mesh_cells) {
        // get points in the cell for tallying
        double x_mid = (cell.surface_bounds[1] - cell.surface_bounds[0]) / 2 + cell.surface_bounds[0];
        double y_mid = (cell.surface_bounds[3] - cell.surface_bounds[2]) / 2 + cell.surface_bounds[2];
        double z_mid = (cell.surface_bounds[5] - cell.surface_bounds[4]) / 2 + cell.surface_bounds[4];
        // get energy deposition from standard tallies
        double ion_energy_dep = garage.standard_tallies[3].retrieve("all", {{"time", timestep_value}, {"x", x_mid}, {"y", y_mid}, {"z", z_mid}});
        double electron_energy_dep = garage.standard_tallies[4].retrieve("e", {{"time", timestep_value}, {"x", x_mid}, {"y", y_mid}, {"z", z_mid}});
        // calculate number density
        double n_ion = 0.0;
        double n_electron = 0.0;
        Material& local_material = cell.cell_material;
        int num_background_species = local_material.species.size();
        for (int ion_it = 0; ion_it < num_background_species; ++ion_it) {
            std::string ion_species = local_material.species[ion_it];
            int ion_z = species_2_z(ion_species);
            double ion_molar_mass = species_2_z(ion_species);
            n_ion += local_material.densities[ion_it] * rtt_units::AVOGADRO / ion_molar_mass;
            n_electron += local_material.densities[ion_it] * rtt_units::AVOGADRO / ion_molar_mass * ion_z;
        }
        n_electron = 3.011e24; // FIX: kludge
        // update assuming ideal monatomic gas
        double dtemp_ion = 2.0/3.0 * ion_energy_dep / n_ion;
        double dtemp_electron = 2.0/3.0 * electron_energy_dep / n_electron;
        if (n_ion == 0.0 || n_electron == 0.0) {
            std::cout << "Error divide by zero: ion/electron density " << std::endl;
        }
        local_material.ion_temperature += dtemp_ion;
        local_material.electron_temperature += dtemp_electron;
        //std::cout << "temp update: " << ion_energy_dep << " " << electron_energy_dep << " " << n_ion << " " << n_electron << " " << dtemp_ion << " " << dtemp_electron << std::endl;

        
        // update standard tallies
        // ion temperature
        garage.standard_tallies[0].add("all", {{"time", timestep_value}, {"x", x_mid}, {"y", y_mid}, {"z", z_mid}}, cell.cell_material.ion_temperature);
        // electron temperature
        garage.standard_tallies[1].add("e", {{"time", timestep_value}, {"x", x_mid}, {"y", y_mid}, {"z", z_mid}}, cell.cell_material.electron_temperature);
        // ion_densities
        for (int species_it = 0; species_it < cell.cell_material.species.size(); ++species_it){
            std::string species = cell.cell_material.species[species_it];
            garage.standard_tallies[2].add(species, {{"time", timestep_value}, {"x", x_mid}, {"y", y_mid}, {"z", z_mid}}, cell.cell_material.densities[species_it]);
        } 
    } // end mesh_cell loop

    // update user tallies
    const std::string specified_category = "average_particle_energy";
    for (auto &t : garage.tallies) {
        if (t.tally_category != specified_category) continue;

        // get average particle energy in the census bank by species
        const auto avgE = average_energy_by_species(garage.census_bank);

        // For each species the tally cares about, if we have that species in the bank,
        // add the *average energy* as the tally amount.
        for (const auto &sp : t.species) {
            auto it = avgE.find(sp);
            if (it == avgE.end()) continue;

            const double avg_energy = it->second;
            const double energy_coordinate = 0.1; 

            // Add the average energy as the "count"/value for this (species, time) bin.
            t.add(sp, {{"time", timestep_value}}, avg_energy);
        }
    }
 
        
    // add the census_bank into the active bank for the next time step
    move_append(garage.active_bank, garage.census_bank);
    
    if (garage.num_t_steps < 101 || (t_it + 1) % 10 == 0) {
        std::cout << "Finished with timestep: " << t_it + 1 << std::endl;
    }
    
}

void transport_particle(Particle& p, double time_census) {

    // energy loss for CSD
    float csd_step = 0.01;

    int zone_index = p.zone;
    Material local_material = garage.mesh_cells[zone_index].cell_material;
    int particle_active = 2; // 2 = alive, 1 = to census, 0 = killed

    while ( particle_active == 2 ){

        //std::cout << "Particle energy keV: " << p.energy << "  time shk: " << p.t << std::endl;
        //std::cout << "Particle position cm: " << p.x << " " << p.y << " " << p.z << std::endl;
    
        // evaluate stopping power dE/dx and dE/dt
        double dedt_electron = 0.0;  // keV/shk
        double dedx_electron = 0.0;  // keV/cm
        double dedt_ion = 0.0;
        double dedx_ion = 0.0;
        //kp_alpha_csd(p); // kp not implemented fully
        spitzer_csd(p, dedt_electron, dedx_electron, dedt_ion, dedx_ion);

        //std::cout << "spitzer parameters: " << dedt_electron << " " << dedx_electron << " " << dedt_ion << " " << dedx_ion << std::endl;

        // time for energy loss
        double t_eloss = (csd_step * p.energy) / (dedt_electron + dedt_ion); // shk
        // time to census
        double t_remaining = time_census - p.t;

        //std::cout << "time parameters: " << t_eloss << " " << t_remaining << std::endl;

        double initial_p_t = p.t;

        // distance to boundary
        BoundaryCross boundary_cross;
        double dist_boundary = 1e20;
        if (distance_to_boundary(p, boundary_cross)) {
            dist_boundary = boundary_cross.distance;
        } else {
            garage.lost_particles++;
            particle_active = 0;
        }
        
        double dist_census = t_remaining * (dedt_electron + dedt_ion) / (dedx_electron + dedx_ion);
        double dist_csd = csd_step * p.energy / (dedx_electron + dedx_ion);

        std::vector<double> distances_to_event = {dist_boundary, dist_csd, dist_census};
        auto it = std::min_element(distances_to_event.begin(), distances_to_event.end());
        double smallest_distance = *it;
        size_t event_index = std::distance(distances_to_event.begin(), it);

        double eloss_total = 0.0;
        if (event_index == 0) {  // boundary event
            double csd_energy_loss = smallest_distance * (dedx_electron + dedx_ion);
            double time_to_boundary = csd_energy_loss / (dedt_electron + dedt_ion);
            eloss_total = csd_energy_loss;
            double boundary_tolerance = 1.0 + 1e-20;
            if (boundary_cross.next_zone == -1) { // vacuum boundary
                //std::cout << "vaccuum boundary at: " << p.x << " " << p.y << " " << p.z << std::endl;
                particle_active = 0;
            } else if (boundary_cross.next_zone == -2) { // reflective boundary
                boundary_tolerance = 1.0 - 1e-10;
                // reverse particle direction
                if (boundary_cross.face == 0 || boundary_cross.face == 1) { // hits x face
                    p.dir.ux = - p.dir.ux;
                } else if (boundary_cross.face == 2 || boundary_cross.face == 3) { // hits y face
                    p.dir.uy = - p.dir.uy;
                } else if (boundary_cross.face == 4 || boundary_cross.face == 5) { // hits z face
                    p.dir.uz = - p.dir.uz;
                }
            } else {
                p.zone = boundary_cross.next_zone;
            }
            // Update particle
            // position
            p.x += smallest_distance * p.dir.ux * boundary_tolerance;
            p.y += smallest_distance * p.dir.uy * boundary_tolerance;
            p.z += smallest_distance * p.dir.uz * boundary_tolerance;
            p.t += time_to_boundary;
            
        } else if (event_index == 1) { // Full csd step
            eloss_total = csd_step * p.energy;
            // update particle
            // position
            p.x += smallest_distance * p.dir.ux;
            p.y += smallest_distance * p.dir.uy;
            p.z += smallest_distance * p.dir.uz;
            p.t += t_eloss; // time

        } else if (event_index == 2) { // census event
            eloss_total = t_remaining * (dedt_electron + dedt_ion);
            // Update particle
            particle_active = 1; // trigger census with 1
            // position
            p.x += smallest_distance * p.dir.ux;
            p.y += smallest_distance * p.dir.uy;
            p.z += smallest_distance * p.dir.uz;
            p.t += t_remaining; // time
        }

        /*
        double eloss_total = 0.0;
        if (t_remaining > t_eloss){ // if the full step can be taken
            //std::cout << "made a complete step" << std::endl;
            // distance for energy loss
            eloss_total = csd_step * p.energy;
            double dist_eloss = eloss_total / (dedx_electron + dedx_ion); // cm
            p.t += t_eloss;
        } else {
            //std::cout << "hit census" << std::endl;
            eloss_total = t_remaining * (dedt_electron + dedt_ion);
            double dist_eloss = eloss_total / (dedx_electron + dedx_ion); // cm
            particle_active = 1;
            p.t += t_remaining;
        }
        */

        // calculate energy loss
        double eloss_electron = eloss_total * (dedt_electron / (dedt_electron + dedt_ion));
        double eloss_ion = eloss_total * (dedt_ion / (dedt_electron + dedt_ion));
        //std::cout << "energy loss: " << eloss_total << std::endl;

        // Tally csd energy loss
        //std::cout << "tally stuff: " << eloss_total << " " << p.weight << " " << initial_p_t << " " << p.energy << std::endl; 
        std::string specified_category = "csd_energy_loss";
        for (auto &t : garage.tallies) {
            if (t.tally_category == specified_category) {
                // Check if this tally tracks this species
                if (std::find(t.species.begin(), t.species.end(), p.species) != t.species.end()) {
                    //std::cout << "here tally" << std::endl;
                    t.add_smear(p.species, {{"time", initial_p_t, p.t}, {"energy", p.energy, p.energy-eloss_total}}, eloss_total*p.weight);
                }
            }
        }
        // ion energy dep
        double ion_eloss_fraction = dedt_ion / (dedt_ion + dedt_electron);
        garage.standard_tallies[3].add("all", {{"time", p.t}, {"x", p.x}, {"y", p.y}, {"z", p.z}}, eloss_total*ion_eloss_fraction*p.weight);
        // electron energy dep
        garage.standard_tallies[4].add("e", {{"time", p.t}, {"x", p.x}, {"y", p.y}, {"z", p.z}}, eloss_total*(1-ion_eloss_fraction)*p.weight);

        // adjust particle energy and speed
        p.energy -= eloss_total;
        double particle_mass = species_2_mass(p.species);
        p.speed = std::sqrt(2.0 * (p.energy*1e-3 * rtt_units::electronChargeSI * 1e6) / (particle_mass * 1.0e-3)) * 1.0e-8 * 1.0e2;
        
        // kill particle if thermalized
        if (p.energy < 1.5*local_material.ion_temperature){
            //std::cout << "thermalized particle" << std::endl;
            particle_active = 0; // kill the particle
            for (int species_it = 0; species_it < local_material.species.size(); ++species_it) {
                if (local_material.species[species_it] == p.species) {
                    double background_modification = p.weight * species_2_mass(p.species);
                    garage.mesh_cells[zone_index].cell_material.densities[species_it] += background_modification;
                }
            }
        }
    }

    // if a particle reaches census, add it to the census_bank
    if (particle_active == 1) {
        //std::cout << "censused particle" << std::endl;
        garage.census_bank.emplace_back(std::move(p));
    }
    // if a secondary is created, add it to the secondary_bank
    
}

static inline bool valid_zone(int zone) {
    return zone >= 0 && zone < static_cast<int>(garage.mesh_cells.size());
}

bool distance_to_boundary(Particle& p, BoundaryCross& out) {
    out.distance  = 0.0;
    out.face      = -1;
    out.next_zone = -1;

    //std::cout << "particle position: " << p.x << " " << p.y << " " << p.z << std::endl;

    if (!valid_zone(p.zone)) return false;

    //std::cout << "dist_2_col: valid zone" << std::endl;

    const MeshCell& cell = garage.mesh_cells[p.zone];
    if (cell.surface_bounds.size() < 6 || cell.boundary_conditions.size() < 6) return false;

    //std::cout << "dist_2_col: valid faces" << std::endl;

    const double ux = p.dir.ux, uy = p.dir.uy, uz = p.dir.uz;

    const double xlo = cell.surface_bounds[0];
    const double xhi = cell.surface_bounds[1];
    const double ylo = cell.surface_bounds[2];
    const double yhi = cell.surface_bounds[3];
    const double zlo = cell.surface_bounds[4];
    const double zhi = cell.surface_bounds[5];

    // Tiny tolerance to avoid “zero step” stalling when starting exactly on a plane.
    const double eps_len = 1e-20;
    const double INF = std::numeric_limits<double>::infinity();

    // coincidence nudging for particles on faces
    nudge_into_cell(p.x, xlo, xhi, eps_len);
    nudge_into_cell(p.y, ylo, yhi, eps_len);
    nudge_into_cell(p.z, zlo, zhi, eps_len);

    //std::cout << "mesh cell boundaries: " << xlo << " " << xhi << " " << ylo << " " << yhi << " " << zlo << " " << zhi << std::endl;
    //std::cout << "particle position: " << p.x << " " << p.y << " " << p.z << std::endl;

    // Time to planes along each axis; use half-open cells [lo, hi), but
    // allow exact max-edge to go to the +face.
    double tx = INF, ty = INF, tz = INF;
    int f_x = -1, f_y = -1, f_z = -1;

    if (ux > 0.0)        { tx = (xhi - p.x) / ux; f_x = 1; }  // +x face (index 1)
    else if (ux < 0.0)   { tx = (xlo - p.x) / ux; f_x = 0; }  // -x face (index 0)

    if (uy > 0.0)        { ty = (yhi - p.y) / uy; f_y = 3; }  // +y face (index 3)
    else if (uy < 0.0)   { ty = (ylo - p.y) / uy; f_y = 2; }  // -y face (index 2)

    if (uz > 0.0)        { tz = (zhi - p.z) / uz; f_z = 5; }  // +z face (index 5)
    else if (uz < 0.0)   { tz = (zlo - p.z) / uz; f_z = 4; }  // -z face (index 4)

    // Only positive forward intersections count.
    if (!(tx > eps_len)) tx = INF;
    if (!(ty > eps_len)) ty = INF;
    if (!(tz > eps_len)) tz = INF;

    // Pick the smallest positive distance; tie-breaker order x < y < z.
    double tmin = tx; int face = f_x;
    if (ty < tmin) { tmin = ty; face = f_y; }
    if (tz < tmin) { tmin = tz; face = f_z; }

    if (!std::isfinite(tmin)) return false; // no forward hit (e.g., dir==0 or already out)

    //std::cout << "dist_2_col: valid forward hit" << std::endl;

    // Pull neighbor / boundary code from the cell
    // Note: your MeshCell.boundary_conditions was shown as vector<double>;
    // cast to int to interpret cell IDs / codes.
    int next_zone = -1;
    {
        double bc_val = cell.boundary_conditions[face];
        // Clamp to int safely
        if (bc_val > static_cast<double>(std::numeric_limits<int>::max()))
            bc_val = static_cast<double>(std::numeric_limits<int>::max());
        if (bc_val < static_cast<double>(std::numeric_limits<int>::min()))
            bc_val = static_cast<double>(std::numeric_limits<int>::min());
        next_zone = static_cast<int>(bc_val);
    }

    out.distance  = tmin;
    out.face      = face;
    out.next_zone = next_zone;
    return true;
}

void nudge_into_cell(double &coord, double lo, double hi, double raw_eps) {
    // Keep the nudge well within the cell if the cell is tiny
    double eps = raw_eps;
    if (hi - lo < eps) {
        eps = (hi - lo) / 20;
    }

    // Very slightly outside on the low side
    if (coord < lo && (lo - coord) <= eps) {
        coord = lo + eps;
    }
    // Very slightly outside on the high side
    else if (coord > hi && (coord - hi) <= eps) {
        coord = hi - eps;
    }
    // Very close inside the low side
    else if (coord - lo < eps) {
        coord = lo + eps;
    }
    // Very close inside the high side
    else if (hi - coord < eps) {
        coord = hi - eps;
    }
};


void spitzer_csd(Particle& p, double& dedt_electron, double& dedx_electron, double& dedt_ion, double& dedx_ion){

    int zone_index = p.zone;
    Material local_material = garage.mesh_cells[zone_index].cell_material;
    int num_background_species = local_material.species.size();

    //std::cout << "local material species and density sizes: " << num_background_species << " " << local_material.densities.size() << std::endl;

    // projectile parameters
    std::string projectile_species = p.species;
    double projectile_mass = species_2_mass(projectile_species);
    int32_t projectile_zaid = species_2_zaid(projectile_species);
    double projectile_energy = p.energy *1e-3 * rtt_units::electronChargeSI * 1e6; // J (the 1e-3 is to convert keV to MeV)
    double projectile_speed = p.speed;

    double background_molar_mass = compute_mixture_molar_mass(local_material.species, local_material.densities);
    garage.mesh_cells[zone_index].cell_material.ion_molar_mass = background_molar_mass;

    double total_ion_density = 0; //  g/cc
    double electron_density = 0; // g/cc * Z

    rtt_units::PhysicalConstexprs<rtt_units::CGS> pc;
    double keV = pc.eV()*1e3;

    for (int ion_it = 0; ion_it < num_background_species; ++ion_it) {
        std::string ion_species = local_material.species[ion_it];
        int32_t ion_zaid = species_2_zaid(ion_species);
        double ion_mass = species_2_mass(ion_species); // g/atom

        double ion_density = local_material.densities[ion_it]; // g/cc
        if (ion_density <= 1e-5){ // 1e-5 g/cc threshold, might need to be adjusted
            break;
        }
        total_ion_density += ion_density;

        electron_density += ion_density * species_2_z(ion_species);

        double ion_number_density = ion_density / background_molar_mass * rtt_units::AVOGADRO;

        auto const ion_model(std::make_shared<rtt_cdi_cpeloss::Analytic_Spitzer_Eloss_Model>(
        ion_zaid, ion_mass, projectile_zaid, projectile_mass));

        rtt_cdi_cpeloss::Analytic_CP_Eloss const eloss_ion(ion_model, rtt_cdi::CPModelAngleCutoff::NONE);

        double dedt_temp = eloss_ion.getEloss(local_material.ion_temperature, ion_number_density, projectile_speed);

        dedt_ion += dedt_temp * ion_density;
        dedx_ion += dedt_temp * ion_density / projectile_speed; //(projectile_speed*1e8 * 1e-8/keV); // keV/cm?

    }

    dedt_ion /= total_ion_density;
    dedx_ion /= total_ion_density;

    // electron slowing down
    double electron_number_density = electron_density / background_molar_mass * rtt_units::AVOGADRO; // atoms/cc
    electron_number_density = 3.011e24; // FIX: kludge
    std::string electron_species = "e";
    double electron_mass = species_2_mass(electron_species);
    auto const electron_model(std::make_shared<rtt_cdi_cpeloss::Analytic_Spitzer_Eloss_Model>(
        -1, electron_mass, projectile_zaid, projectile_mass
    ));
    rtt_cdi_cpeloss::Analytic_CP_Eloss const eloss_electron(electron_model, rtt_cdi::CPModelAngleCutoff::NONE);
    dedt_electron = eloss_electron.getEloss(local_material.electron_temperature, electron_number_density, projectile_speed); // keV/shk
    dedx_electron = dedt_electron / projectile_speed; //(projectile_speed*1e8 * 1e-8/keV);  // keV/cm?

}


std::unordered_map<std::string, double>
average_energy_by_species(const std::vector<Particle>& bank)
{
    std::unordered_map<std::string, double> sum_wE;
    std::unordered_map<std::string, double> sum_w;

    for (const auto &p : bank) {
        sum_wE[p.species] += p.weight * p.energy;
        sum_w[p.species]  += p.weight;
    }

    std::unordered_map<std::string, double> avg;
    avg.reserve(sum_w.size());
    for (const auto &kv : sum_w) {
        const std::string &sp = kv.first;
        const double w = kv.second;
        if (w > 0.0) {
            avg[sp] = sum_wE[sp] / w;
        }
    }
    return avg;
}


void source_particles(double time_start, double time_census) {
    
    double particle_weight = (garage.source_strength/garage.num_particles);

    if (garage.source_time >= time_start && garage.source_time < time_census){ // this only works for point source in time
        //std::cout << "sourced particle: " << time_start << " " << time_census << " " << garage.source_time << std::endl;
        for (int source_it = 0; source_it < garage.num_particles; ++source_it){
            Particle p;
            p.species = garage.source_particle;
            // location
            p.x = garage.source_point[0];
            p.y = garage.source_point[1];
            p.z = garage.source_point[2];
            // direction
            p.dir.ux = garage.source_direction[0];
            p.dir.uy = garage.source_direction[1];
            p.dir.uz = garage.source_direction[2];
            p.dir.normalize();
            // other
            p.t = garage.source_time;
            p.energy = garage.source_energy;
            p.weight = particle_weight;
            double projectile_mass = species_2_mass(p.species);
            p.speed = particle_energy_2_speed(p.energy, projectile_mass);
            p.zone = locate_cell_id(p.x, p.y, p.z);
            p.start_particle_rng();
            p.id = garage.total_particles_created;
            garage.total_particles_created++;
            
            garage.active_bank.emplace_back(p);
        }
    }
    
}


static double total_weight(const std::vector<Particle>& v) {
    return std::accumulate(v.begin(), v.end(), 0.0,
        [](double s, const Particle& p){ return s + p.weight; });
}


void population_control() {
    auto& v = garage.active_bank;
    const int target_i = std::max(garage.num_particles, 0);
    const std::size_t target = static_cast<std::size_t>(target_i);

    if (v.size() == target) return;
    
    const double W_before = total_weight(v);

    auto& gen = rng();

    if (v.size() > target) {
        // Randomly eliminate down to target (Fisher�Yates tail shrink)
        // Shuffle a suffix into the tail, then resize.
        // Implementation: do a partial shuffle for k = v.size()-target steps.
        const std::size_t n = v.size();
        std::uniform_int_distribution<std::size_t> pick(0, n - 1);
        // Simpler, fine for now: shuffle entire vector and resize
        std::shuffle(v.begin(), v.end(), gen);
        v.resize(target);
    } else {
        // Duplicate random particles until we reach target
        std::uniform_int_distribution<std::size_t> pick(0, v.empty() ? 0 : v.size() - 1);
        while (v.size() < target) {
            const std::size_t idx = pick(gen);
            v.push_back(v[idx]); // copy-duplicate
            // update distribution upper bound as size grows
            pick = std::uniform_int_distribution<std::size_t>(0, v.size() - 1);
        }
    }
    
    double W_after = total_weight(v);
    const double scale = W_before / W_after;
    for (auto& p : v) p.weight *= scale;
}


void time_step_setup(std::string mode){
    if (mode == "Uniform"){
        double step_size = garage.t_step_size;
        int num_steps = garage.num_t_steps;
        double starting_point = 0.0;

                // Create a vector of size num_steps + 1
        std::vector<double> time_bins(num_steps + 1);

        // Fill it with uniformly spaced values
        for (int i = 0; i <= num_steps; ++i) {
            time_bins[i] = starting_point + i * step_size;
        }

        garage.time_step_bins = time_bins;
    }
}



// Find axis index i such that value ∈ [edges[i], edges[i+1])
// Special-case: if value == edges.back(), return last cell (size-2).
// Returns -1 if value is outside [edges.front(), edges.back()].
static inline int find_axis_index(const std::vector<double> &edges, double value) {
    const size_t n = edges.size();
    if (n < 2) return -1; // no cells
    if (value < edges.front() || value > edges.back()) return -1;
    if (value == edges.back()) return static_cast<int>(n) - 2;

    auto it = std::upper_bound(edges.begin(), edges.end(), value);
    int i = static_cast<int>(it - edges.begin()) - 1;
    // i is now in [0, n-2] if value is inside
    if (i < 0 || i >= static_cast<int>(n) - 1) return -1;
    return i;
}

bool locate_cell_indices(double x, double y, double z, int &ix, int &iy, int &iz) {
    // Require a built mesh
    if (garage.x_bin_bounds.size() < 2 ||
        garage.y_bin_bounds.size() < 2 ||
        garage.z_bin_bounds.size() < 2 ||
        garage.nx <= 0 || garage.ny <= 0 || garage.nz <= 0) {
        return false;
    }

    ix = find_axis_index(garage.x_bin_bounds, x);
    if (ix < 0) return false;
    iy = find_axis_index(garage.y_bin_bounds, y);
    if (iy < 0) return false;
    iz = find_axis_index(garage.z_bin_bounds, z);
    if (iz < 0) return false;

    return true;
}

int locate_cell_id(double x, double y, double z) {
    int ix, iy, iz;
    if (!locate_cell_indices(x, y, z, ix, iy, iz)) return -1;

    // Row-major with x fastest:
    // id = ix + iy*nx + iz*(nx*ny)
    const int nx = garage.nx;
    const int ny = garage.ny;
    // (nx, ny, nz were set during mesh build)
    return ix + iy * nx + iz * (nx * ny);
}