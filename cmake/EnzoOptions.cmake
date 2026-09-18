# cmake/EnzoOptions.cmake
# Configurable options for Enzo build system

# Parameters
set(ENZO_MAX_SUBGRIDS "100000" CACHE STRING "Maximum number of subgrids")
set(ENZO_MAX_BARYONS "30" CACHE STRING "Maximum number of baryon fields")
set(ENZO_MAX_TASKS_PER_NODE "8" CACHE STRING "Maximum tasks per node")
set(ENZO_MEMORY_POOL_SIZE "100000" CACHE STRING "Memory pool size")

# Numeric Precisions
set(ENZO_PRECISION "64" CACHE STRING "Baryon floating point precision (32, 64)")
set(ENZO_PARTICLES "64" CACHE STRING "Particle position precision (32, 64, 128)")
set(ENZO_INTEGERS "64" CACHE STRING "Integer size (32, 64)")
set(ENZO_PARTICLE_IDS "64" CACHE STRING "Particle ID integer size (32, 64)")
set(ENZO_INITS_PRECISION "64" CACHE STRING "Inits precision (32, 64)")
set(ENZO_IO_PRECISION "32" CACHE STRING "IO precision (32, 64)")

# Feature Flags
option(ENZO_USE_MPI "Enable MPI support" ON)
option(ENZO_TASKMAP "Enable taskmap" OFF)
option(ENZO_PACKED_AMR "Enable packed AMR" ON)
option(ENZO_PACKED_MEM "Enable packed memory with packed AMR" OFF)
option(ENZO_USE_LCAPERF "Enable LCAPERF profiling" OFF)
option(ENZO_USE_PAPI "Enable PAPI hardware counters with LCAPERF" OFF)
option(ENZO_USE_PYTHON "Enable embedded Python interpreter" OFF)
option(ENZO_USE_LIBYT "Enable libyt in situ analysis" OFF)
option(ENZO_LIBYT_INTERACTIVE "Enable libyt interactive prompt" OFF)
option(ENZO_LIBYT_RELOAD "Enable libyt reloading script feature" OFF)
option(ENZO_LIBYT_JUPYTER "Enable libyt Jupyter Notebook feature" OFF)
option(ENZO_NEW_PROBLEM_TYPES "Enable new problem initializers" OFF)
option(ENZO_OOC_BOUNDARY "Enable out-of-core boundary conditions" OFF)
option(ENZO_ACCELERATION_BOUNDARY "Enable setting acceleration boundary" ON)
option(ENZO_TESTING "Enable hooks for test suites" OFF)
option(ENZO_PHOTON "Enable adaptive ray tracing for radiative transfer" ON)
option(ENZO_USE_HYPRE "Enable HYPRE bindings for implicit solvers" OFF)
option(ENZO_EMISSIVITY "Enable emissivity field allocation" OFF)
option(ENZO_NEW_GRID_IO "Enable new simpler grid IO routines" ON)
option(ENZO_FAST_SIB "Enable fast sibling locator" ON)
option(ENZO_USE_HDF4 "Use HDF4 instead of HDF5" OFF)
option(ENZO_BITWISE_IDENTICALITY "Use blocking potential solves for bitwise identicality" OFF)
option(ENZO_USE_CUDA "Enable CUDA GPU acceleration" OFF)
option(ENZO_GRAVITY_4S "Use 4th order gravity solver" OFF)
option(ENZO_USE_GRACKLE "Enable Grackle chemistry and cooling library" OFF)
option(ENZO_ENZO_PERFORMANCE "Enable performance and timing measurements" ON)
option(ENZO_LOG2ALLOC "Enable power of 2 block size allocations" OFF)
option(ENZO_USE_UUID "Enable UUID functionality" ON)

# Assemble preprocessor definitions into global variable ENZO_COMPILE_DEFINITIONS
set(ENZO_COMPILE_DEFINITIONS
    LINUX
    H5_USE_16_API
    HAVE_UNISTD_H
    HAVE_SYS_TIME_H
    HAVE_STDLIB_H
    __max_subgrids=${ENZO_MAX_SUBGRIDS}
    __max_baryons=${ENZO_MAX_BARYONS}
    __max_cpu_per_node=${ENZO_MAX_TASKS_PER_NODE}
    __memory_pool_size=${ENZO_MEMORY_POOL_SIZE}
)

if (ENZO_PRECISION EQUAL 32)
    list(APPEND ENZO_COMPILE_DEFINITIONS CONFIG_BFLOAT_4)
elseif (ENZO_PRECISION EQUAL 64)
    list(APPEND ENZO_COMPILE_DEFINITIONS CONFIG_BFLOAT_8)
else()
    message(FATAL_ERROR "Invalid ENZO_PRECISION value '${ENZO_PRECISION}'. Must be 32 or 64.")
endif()

if (ENZO_PARTICLES EQUAL 32)
    list(APPEND ENZO_COMPILE_DEFINITIONS CONFIG_PFLOAT_4)
elseif (ENZO_PARTICLES EQUAL 64)
    list(APPEND ENZO_COMPILE_DEFINITIONS CONFIG_PFLOAT_8)
elseif (ENZO_PARTICLES EQUAL 128)
    list(APPEND ENZO_COMPILE_DEFINITIONS CONFIG_PFLOAT_16)
else()
    message(FATAL_ERROR "Invalid ENZO_PARTICLES value '${ENZO_PARTICLES}'. Must be 32, 64, or 128.")
endif()

if (ENZO_INTEGERS EQUAL 32)
    list(APPEND ENZO_COMPILE_DEFINITIONS SMALL_INTS)
elseif (ENZO_INTEGERS EQUAL 64)
    list(APPEND ENZO_COMPILE_DEFINITIONS LARGE_INTS)
else()
    message(FATAL_ERROR "Invalid ENZO_INTEGERS value '${ENZO_INTEGERS}'. Must be 32 or 64.")
endif()

if (ENZO_PARTICLE_IDS EQUAL 32)
    list(APPEND ENZO_COMPILE_DEFINITIONS CONFIG_PINT_4)
elseif (ENZO_PARTICLE_IDS EQUAL 64)
    list(APPEND ENZO_COMPILE_DEFINITIONS CONFIG_PINT_8)
else()
    message(FATAL_ERROR "Invalid ENZO_PARTICLE_IDS value '${ENZO_PARTICLE_IDS}'. Must be 32 or 64.")
endif()

if (ENZO_INITS_PRECISION EQUAL 32)
    list(APPEND ENZO_COMPILE_DEFINITIONS INITS32)
elseif (ENZO_INITS_PRECISION EQUAL 64)
    list(APPEND ENZO_COMPILE_DEFINITIONS INITS64)
else()
    message(FATAL_ERROR "Invalid ENZO_INITS_PRECISION value '${ENZO_INITS_PRECISION}'. Must be 32 or 64.")
endif()

if (ENZO_IO_PRECISION EQUAL 32)
    list(APPEND ENZO_COMPILE_DEFINITIONS IO_32)
elseif (ENZO_IO_PRECISION EQUAL 64)
    list(APPEND ENZO_COMPILE_DEFINITIONS IO_64)
else()
    message(FATAL_ERROR "Invalid ENZO_IO_PRECISION value '${ENZO_IO_PRECISION}'. Must be 32 or 64.")
endif()

if (ENZO_USE_MPI)
    list(APPEND ENZO_COMPILE_DEFINITIONS USE_MPI)
endif()

if (ENZO_TASKMAP)
    list(APPEND ENZO_COMPILE_DEFINITIONS TASKMAP ENABLE_TASKMAP)
endif()

if (ENZO_PACKED_AMR)
    list(APPEND ENZO_COMPILE_DEFINITIONS USE_HDF5_GROUPS)
endif()

if (ENZO_PACKED_MEM)
    list(APPEND ENZO_COMPILE_DEFINITIONS USE_HDF5_OUTPUT_BUFFERING)
endif()

if (ENZO_USE_LCAPERF)
    list(APPEND ENZO_COMPILE_DEFINITIONS USE_LCAPERF)
endif()

if (ENZO_USE_PAPI)
    list(APPEND ENZO_COMPILE_DEFINITIONS USE_PAPI)
endif()

if (ENZO_USE_PYTHON)
    list(APPEND ENZO_COMPILE_DEFINITIONS USE_PYTHON)
endif()

if (ENZO_USE_LIBYT)
    list(APPEND ENZO_COMPILE_DEFINITIONS USE_LIBYT)
    if (ENZO_LIBYT_INTERACTIVE)
        list(APPEND ENZO_COMPILE_DEFINITIONS USE_LIBYT_INTERACTIVE)
    endif()
    if (ENZO_LIBYT_RELOAD)
        list(APPEND ENZO_COMPILE_DEFINITIONS USE_LIBYT_RELOAD)
    endif()
    if (ENZO_LIBYT_JUPYTER)
        list(APPEND ENZO_COMPILE_DEFINITIONS USE_LIBYT_JUPYTER)
    endif()
endif()

if (ENZO_NEW_PROBLEM_TYPES)
    list(APPEND ENZO_COMPILE_DEFINITIONS NEW_PROBLEM_TYPES)
endif()

if (ENZO_OOC_BOUNDARY)
    list(APPEND ENZO_COMPILE_DEFINITIONS OOC_BOUNDARY)
endif()

if (ENZO_ACCELERATION_BOUNDARY)
    list(APPEND ENZO_COMPILE_DEFINITIONS SAB)
endif()

if (ENZO_TESTING)
    list(APPEND ENZO_COMPILE_DEFINITIONS CONFIG_TESTING)
endif()

if (ENZO_PHOTON)
    list(APPEND ENZO_COMPILE_DEFINITIONS TRANSFER)
endif()

if (ENZO_USE_HYPRE)
    list(APPEND ENZO_COMPILE_DEFINITIONS USE_HYPRE)
endif()

if (ENZO_EMISSIVITY)
    list(APPEND ENZO_COMPILE_DEFINITIONS EMISSIVITY)
endif()

if (ENZO_NEW_GRID_IO)
    list(APPEND ENZO_COMPILE_DEFINITIONS NEW_GRID_IO)
endif()

if (ENZO_FAST_SIB)
    list(APPEND ENZO_COMPILE_DEFINITIONS FAST_SIB)
endif()

if (ENZO_USE_HDF4)
    list(APPEND ENZO_COMPILE_DEFINITIONS USE_HDF4)
endif()

if (ENZO_BITWISE_IDENTICALITY)
    list(APPEND ENZO_COMPILE_DEFINITIONS BITWISE_IDENTICALITY)
endif()

if (ENZO_USE_CUDA)
    list(APPEND ENZO_COMPILE_DEFINITIONS ECUDA)
endif()

if (ENZO_GRAVITY_4S)
    list(APPEND ENZO_COMPILE_DEFINITIONS GRAVITY_4S)
endif()

if (ENZO_ENZO_PERFORMANCE)
    list(APPEND ENZO_COMPILE_DEFINITIONS ENZO_PERFORMANCE)
endif()

if (ENZO_USE_GRACKLE)
    list(APPEND ENZO_COMPILE_DEFINITIONS USE_GRACKLE)
endif()

if (ENZO_LOG2ALLOC)
    list(APPEND ENZO_COMPILE_DEFINITIONS USE_LOG2ALLOC)
endif()

if (ENZO_USE_UUID)
    list(APPEND ENZO_COMPILE_DEFINITIONS USE_UUID)
endif()
