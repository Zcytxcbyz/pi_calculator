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
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O2 /openmp")
    else()
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -fopenmp")
    endif()
endif()

# LTO support
option(ENABLE_LTO "Enable Link Time Optimization" ON)
if(ENABLE_LTO)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    if(MSVC)
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /GL")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LTCG")
    else()
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -flto")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
    endif()
endif()

# SIMD support
option(ENABLE_SIMD "Enable SIMD instructions (AVX2/SSE4)" ON)
if(ENABLE_SIMD AND CMAKE_BUILD_TYPE STREQUAL "Release")
    if(MSVC)
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /arch:AVX2")
    else()
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mavx2")
    endif()
endif()

# Enable Fast Math
option(ENABLE_FAST_MATH "Enable Fast Math" ON)
if(ENABLE_FAST_MATH)
    if(MSVC)
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /fp:fast")
    else()
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -ffast-math")
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
