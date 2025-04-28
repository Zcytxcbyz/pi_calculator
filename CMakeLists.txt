cmake_minimum_required(VERSION 3.12)
project(pi_calculator C)

# Basic configuration
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set the path of the CMake module
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Dependency lookups
find_package(GMP REQUIRED)
find_package(OpenMP REQUIRED)

# User options configuration
option(ENABLE_SIMD "Enable SIMD instructions (AVX2/SSE4)" ON)
option(ENABLE_LTO "Enable Link Time Optimization (LTO)" ON)
option(BUILD_STATIC "Build as static executable" OFF)
option(ENABLE_CACHE "Enable cache for large calculations" ON)
option(ENABLE_BLOCK_FACTORIAL "Enable block factorial optimization" ON)

# Compiler optimization configuration
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    # Basic optimization
    if(MSVC)
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O2 /DNDEBUG /fp:fast /fopenmp")
    else()
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -DNDEBUG -ffast-math -fopenmp")
    endif()

    # SIMD instruction set optimization
    if(ENABLE_SIMD)
        if(MSVC)
            set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /arch:AVX2")
        else()
        include(CheckCCompilerFlag)
            check_c_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
            if(COMPILER_SUPPORTS_MARCH_NATIVE)
                set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native")
            else()
                set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mavx2 -mbmi2 -mfma")
            endif()
        endif()
    endif()

    # Link-Time Optimization (LTO)
    if(ENABLE_LTO)
        include(CheckIPOSupported)
        check_ipo_supported(RESULT LTO_SUPPORTED OUTPUT output)
        if(LTO_SUPPORTED)
            set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
            if(MSVC)
                set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /GL")
                set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
            else()
                set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -flto")
                set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
            endif()
        else()
            message(WARNING "LTO not supported: ${output}")
        endif()
    endif()
else()
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DDEBUG -fopenmp")
    # Debug mode configuration
    add_compile_definitions(DEBUG)
endif()

# Static link configuration
if(BUILD_STATIC)
    if(MSVC)
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
        set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
    else()
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
    endif()
endif()

# Cache configuration
if(ENABLE_CACHE)
    add_compile_definitions(ENABLE_CACHE)
endif()

# Block factorial optimization configuration
if(ENABLE_BLOCK_FACTORIAL)
    add_compile_definitions(ENABLE_BLOCK_FACTORIAL)
endif()

# Executable file configuration
add_executable(pi_calculator
    src/pi.c
    main.c
)

# Include directories
target_include_directories(pi_calculator PRIVATE
    include
)

# Link libraries
target_link_libraries(pi_calculator PRIVATE
    GMP::GMP
    OpenMP::OpenMP_C
)

# Define macros for copying dependent libraries
function(copy_windows_dependencies target)
    if(WIN32 AND CMAKE_C_COMPILER_ID STREQUAL "GNU")  # Only handle MinGW environments
        # Retrieve the bin directory of MinGW (inferred from compiler path)
        get_filename_component(MINGW_BIN_DIR ${CMAKE_C_COMPILER} DIRECTORY)
        get_filename_component(MINGW_BIN_DIR ${MINGW_BIN_DIR}/../bin ABSOLUTE)

        # List of DLLs that need to be copied
        set(REQUIRED_DLLS
            "${MINGW_BIN_DIR}/libgmp-10.dll"   # GMP dynamic library
            "${MINGW_BIN_DIR}/libgcc_s_seh-1.dll"
            "${MINGW_BIN_DIR}/libgomp-1.dll"
            "${MINGW_BIN_DIR}/libwinpthread-1.dll"
        )

        # Add copy command
        foreach(DLL ${REQUIRED_DLLS})
            add_custom_command(TARGET ${target} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${DLL}"
                "$<TARGET_FILE_DIR:${target}>"
                COMMENT "Copying ${DLL} to output directory..."
            )
        endforeach()
    endif()
endfunction()

if(NOT BUILD_STATIC)
    # Call Windows dependency copy macro
    copy_windows_dependencies(pi_calculator)
endif()

# Installation rules
install(TARGETS pi_calculator DESTINATION bin)
