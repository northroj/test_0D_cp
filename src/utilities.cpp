#include "utilities.h"

#include <iostream>
#include <cstdlib>
#include <string>

// HDF5 C API is a C library; wrap in extern "C" for C++ builds
extern "C" {
#include <hdf5.h>
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
