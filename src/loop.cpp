#include "loop.h"
#include "garage.h"
#include "utilities.h"

#include <algorithm>
#include <iterator>
#include <random>
#include <utility>
#include <iostream>


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

            transport_particle(p);
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

void transport_particle(Particle& p) {

    // if a particle reaches census, add it to the census_bank
    if (true) {
        garage.census_bank.emplace_back(std::move(p));
    }
    // if a secondary is created, add it to the secondary_bank
    
}

void source_particles(int t_it) {

    garage.active_bank.emplace_back(
        garage.garage.source_particle,
        garage.source_point[0],
        garage.source_point[1],
        garage.source_point[2],
        garage.source_time,
        garage.source_energy,
    );
    
}

void population_control() {
    auto& v = garage.active_bank;
    const int target_i = std::max(garage.num_particles, 0);
    const std::size_t target = static_cast<std::size_t>(target_i);

    if (v.size() == target) return;

    auto& gen = rng();

    if (v.size() > target) {
        // Randomly eliminate down to target (Fisher–Yates tail shrink)
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
}
