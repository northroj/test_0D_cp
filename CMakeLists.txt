cmake_minimum_required(VERSION 3.18)
project(galvanize LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the installed Draco package
set(draco_DIR "$ENV{HOME}/lifft/draco_install/lib/cmake/draco")
find_package(draco CONFIG REQUIRED)

find_package(MPI  REQUIRED COMPONENTS C CXX Fortran)
find_package(GSL  REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C CXX)

# Helper: link an executable to many Draco targets if they exist in this install
function(link_if_exists exe)
  foreach(cand IN LISTS ARGN)
    if(TARGET ${cand})
      message(STATUS "Linking ${exe} to Draco target: ${cand}")
      target_link_libraries(${exe} PRIVATE ${cand})
    else()
      message(STATUS "Draco target not present (skipping): ${cand}")
    endif()
  endforeach()
endfunction()

add_executable(galvanize
  src/main.cpp
  src/loop.cpp
  src/utilities.cpp
  src/io.cpp
)

if (TARGET hdf5::hdf5)
  target_link_libraries(galvanize PRIVATE hdf5::hdf5)
elseif (TARGET HDF5::HDF5)
  target_link_libraries(galvanize PRIVATE HDF5::HDF5)
else()
  # Legacy FindHDF5 variables
  target_include_directories(galvanize PRIVATE ${HDF5_INCLUDE_DIRS})
  if (HDF5_C_LIBRARIES)
    target_link_libraries(galvanize PRIVATE ${HDF5_C_LIBRARIES})
  else()
    target_link_libraries(galvanize PRIVATE ${HDF5_LIBRARIES})
  endif()
endif()

# Link against all Draco library targets present in the file (skip Exe_*).
# (Order here doesn�t matter�CMake handles transitive links via imported target metadata.)
link_if_exists(galvanize
  coverage_config
  Lib_dsxx
  Lib_c4
  Lib_c4_fc4
  Lib_cdi
  Lib_linear
  Lib_memory
  Lib_mesh_element
  Lib_units
  Lib_cdi_CPEloss
  Lib_cdi_ipcress
  Lib_diagnostics
  Lib_kde
  Lib_meshReaders
  Lib_min
  Lib_norms
  Lib_parser
  Lib_rng
  Lib_roots
  Lib_viz
  Lib_predict
  Lib_lightweight_quadrature
  Lib_cdi_analytic
  Lib_cdi_ndi
  Lib_cdi_cpt
  Lib_RTT_Format_Reader
  Lib_special_functions
  Lib_compton_tools
  Lib_mesh
  Lib_quadrature
  Lib_quadrature_fquadrature
)

# (We intentionally do NOT link imported executables:)
# Exe_xthi Exe_ythi Exe_Ipcress_Interpreter Exe_draco_info Exe_CSKrw

# Optional: set an rpath so your built binary can find Draco�s .so without LD_LIBRARY_PATH
# (adjust lib vs lib64 to match your install)
set_target_properties(galvanize PROPERTIES
  BUILD_RPATH "$ENV{HOME}/lifft/draco_install/lib;$ENV{HOME}/lifft/draco_install/lib64"
  INSTALL_RPATH "$ENV{HOME}/lifft/draco_install/lib;$ENV{HOME}/lifft/draco_install/lib64"
)

# If you later install this project:
# install(TARGETS galvanize RUNTIME DESTINATION bin)