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

    // Loop over timesteps
  
    for(int t_it = 0; t_it < garage.num_t_steps; ++t_it) {
        garage.current_time_step = t_it;
        simulate_timestep(t_it);
    }
    
}

void simulate_timestep(int t_it) {
    source_particles(t_it);
    population_control();

    // calculate the end of timestep time
    double time_census = (t_it+1) * garage.t_step_size; // shk
    
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
 
        
    // add the census_bank into the active bank for the next time step
    move_append(garage.active_bank, garage.census_bank);
    
    if (garage.num_t_steps < 200 || (t_it + 1) % 10 == 0) {
        std::cout << "Finished with timestep: " << t_it + 1 << std::endl;
    }
    
}

void transport_particle(Particle& p, double time_census) {

    // energy loss for CSD
    float csd_step = 0.01;

    int zone_index = p.zone;
    Material local_material = garage.materials[zone_index];
    int particle_active = 2; // 2 = alive, 1 = to census, 0 = killed

    while ( particle_active == 2 ){
    
        // evaluate stopping power dE/dx and dE/dt
        double dedt_electron = 0.0;  // keV/shk
        double dedx_electron = 0.0;  // keV/cm
        double dedt_ion = 0.0;
        double dedx_ion = 0.0;
        //kp_alpha_csd(p); // kp not implemented fully
        spitzer_csd(p, dedt_electron, dedx_electron, dedt_ion, dedx_ion);

        // time for energy loss
        double t_eloss = (csd_step * p.energy * 1e3 /*MeV -> keV*/ ) / (dedt_electron + dedt_ion); // shk
        // time to census
        double t_remaining = time_census - p.t;

        double eloss_total = 0.0;
        if (t_remaining > t_eloss){ // if the full step can be taken
            // distance for energy loss
            eloss_total = csd_step * p.energy * 1e3 /*MeV -> keV*/;
            double dist_eloss = eloss_total / (dedx_electron + dedx_ion); // cm
        } else {
            eloss_total = t_remaining * (dedt_electron + dedt_ion);
            double dist_eloss = eloss_total / (dedx_electron + dedx_ion); // cm
            particle_active = 1;
        }

        // calculate energy loss
        double eloss_electron = eloss_total * (dedt_electron / (dedt_electron + dedt_ion));
        double eloss_ion = eloss_total * (dedt_ion / (dedt_electron + dedt_ion));

        p.energy -= eloss_total;
        if (p.energy*1e3 < 1.5*local_material.ion_temperature){
            particle_active = 0; // kill the particle if it thermalizes
            // TODO: add the particle to the background
        }
    }

    // if a particle reaches census, add it to the census_bank
    if (particle_active == 1) {
        garage.census_bank.emplace_back(std::move(p));
    }
    // if a secondary is created, add it to the secondary_bank
    
}


void kp_alpha_csd(Particle& p) {

    int zone_index = p.zone;
    
    Material local_material = garage.materials[zone_index];
    int num_background_species = local_material.species.size();

    // Single background species, whatever is listed first in the input
    std::string ion_species = local_material.species[0];
    int32_t ion_zaid = species_2_zaid(ion_species);
    double ion_mass = species_2_mass(ion_species);

    // projectile values
    std::string projectile_species = p.species;
    int32_t projectile_zaid = species_2_zaid(projectile_species);
    double projectile_mass = species_2_mass(projectile_species);

    // KP for alphas slowing down in d + t
    auto const model_in(std::make_shared<rtt_cdi_cpeloss::Analytic_KP_Alpha_Eloss_Model>(
      ion_zaid, ion_mass, projectile_zaid, projectile_mass));
    rtt_cdi_cpeloss::Analytic_CP_Eloss const eloss_mod(model_in, rtt_cdi::CPModelAngleCutoff::NONE);
    
    double eloss_coeff = eloss_mod.getEloss(
        local_material.ion_temperature, // keV
        local_material.densities[0] / ion_mass, // [g/cc]/[g/atom] = atoms/cc
        1.0 // cm/shk
    );


}

void spitzer_csd(Particle& p, double& dedt_electron, double& dedx_electron, double& dedt_ion, double& dedx_ion){

    int zone_index = p.zone;
    Material local_material = garage.materials[zone_index];
    int num_background_species = local_material.species.size();

    // projectile parameters
    std::string projectile_species = p.species;
    double projectile_mass = species_2_mass(projectile_species);
    int32_t projectile_zaid = species_2_zaid(projectile_species);
    double projectile_energy = p.energy * rtt_units::electronChargeSI * 1e6; // J
    double projectile_speed = p.speed;

    double background_molar_mass = compute_mixture_molar_mass(local_material.species, local_material.densities);
    garage.materials[zone_index].ion_molar_mass = background_molar_mass;

    double total_ion_density = 0; //  g/cc

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

        double ion_number_density = ion_density / background_molar_mass * rtt_units::AVOGADRO;

        auto const ion_model(std::make_shared<rtt_cdi_cpeloss::Analytic_Spitzer_Eloss_Model>(
        ion_zaid, ion_mass, projectile_zaid, projectile_mass));

        rtt_cdi_cpeloss::Analytic_CP_Eloss const eloss_ion(ion_model, rtt_cdi::CPModelAngleCutoff::NONE);

        double dedt_temp = eloss_ion.getEloss(local_material.ion_temperature, ion_number_density, projectile_speed);

        dedt_ion += dedt_temp;
        dedx_ion += dedt_temp / (projectile_speed*1e8 * 1e-8/keV); // keV/cm?

    }

    // electron slowing down
    double electron_number_density = total_ion_density / background_molar_mass * rtt_units::AVOGADRO; // atoms/cc
    double electron_mass = 9.1093837e-28;
    auto const electron_model(std::make_shared<rtt_cdi_cpeloss::Analytic_Spitzer_Eloss_Model>(
        -1, electron_mass, projectile_zaid, projectile_mass
    ));
    rtt_cdi_cpeloss::Analytic_CP_Eloss const eloss_electron(electron_model, rtt_cdi::CPModelAngleCutoff::NONE);
    dedt_electron = eloss_electron.getEloss(local_material.electron_temperature, electron_number_density, projectile_speed); // keV/shk
    dedx_electron = dedt_electron / (projectile_speed*1e8 * 1e-8/keV);  // keV/cm?

}


void source_particles(int t_it) {

    garage.active_bank.emplace_back(
        garage.source_particle,
        garage.source_point[0],
        garage.source_point[1],
        garage.source_point[2],
        garage.source_time,
        garage.source_energy,
        1.0
    );
    
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
