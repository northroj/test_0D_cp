// main.cpp
#include <iostream>
#include <chrono>
#include <filesystem>
#include <cstdlib>

#include "utilities.h"
#include "loop.h"
#include "garage.h"
#include "io.h"    

// ---- Global container: garage ----------------------------------------------
Garage garage;

int main(int argc, char** argv) try {
    using clock_t = std::chrono::steady_clock;

    std::cout << "Preparing simulation\n";

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
    simulate();
    const auto t1_sim = clock_t::now();

    const auto t1_total = clock_t::now();

    const double sim_s   = std::chrono::duration<double>(t1_sim   - t0_sim  ).count();
    const double total_s = std::chrono::duration<double>(t1_total - t0_total).count();

    std::cout.setf(std::ios::fixed);
    std::cout.precision(6);
    std::cout << "Simulation time: " << sim_s   << " s\n";
    std::cout << "Total time     : " << total_s << " s\n";

    if (!write_output(cli.output, sim_s, total_s)) {
        return EXIT_FAILURE;
    }

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