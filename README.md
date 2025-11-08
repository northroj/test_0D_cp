# galvanize

<p align="center">
  <img src="docs/assets/0076Golem-Alola.png" alt="alolan golem" width="200"/>
</p>

Infant code for charged particle transport using Draco


## Basic setup

1) Install Draco (from their github, lanl/draco, or using spack. Current instructions are for using a custom install of tag draco-7_21_0)

2) Clone galvanize:

        git clone https://github.com/northroj/galvanize.git

3) Modify CMakeLists.txt to properly point to draco

4) Install galvanize

        From the home directory:
        cmake -S . -B build -DCMAKE_PREFIX_PATH=path/to/draco_install
        cmake --build build


## Basic usage

Short example input deck in examples/example_basic.txt

To run:

        cd examples
        ../build/galvanize -i example_basic.txt

It will run and create an output hdf5 file called example_basic.h5 (use the -o flag to specify an output name)



## Input explanation

Check out https://northroj.github.io/galvanize/ for the user documentation on input formatting (work in progress)


