cmake_minimum_required(VERSION 3.12)
project(pi_calculator C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Set the path of the CMake module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Search for dependency libraries
find_package(GMP REQUIRED)
find_package(OpenMP REQUIRED)

# Compiler optimization
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(MSVC)
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O2 /fp:fast /GL /openmp")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LTCG")
    else()
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -ffast-math -march=native -flto -fopenmp")
    endif()
endif()

# Static compilation options
option(BUILD_STATIC "Build as static executable" ON)
if(BUILD_STATIC)
    if(MSVC)
        # Use /MT for static runtime linking
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
        set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
    else()
        # Use GCC/Clang specific static linking flags
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
    endif()
endif()

# Add include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Add executable file
add_executable(pi_calculator
    src/pi.c
    main.c
)

# Add include directories
target_include_directories(pi_calculator PRIVATE include)

# Link libraries
target_link_libraries(pi_calculator
    PRIVATE
    GMP::GMP
    OpenMP::OpenMP_C
)
