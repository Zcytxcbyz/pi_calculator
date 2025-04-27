# cmake/FindGMP.cmake
find_path(GMP_INCLUDE_DIR gmp.h) # Find the GMP header file
find_library(GMP_LIBRARY NAMES gmp) # Find the GMP library

include(FindPackageHandleStandardArgs) # Include the standard arguments handling module
find_package_handle_standard_args(GMP DEFAULT_MSG GMP_INCLUDE_DIR GMP_LIBRARY) # Handle standard arguments for the package

# Check if the GMP library was found
if(GMP_FOUND)
    message(STATUS "Found GMP: ${GMP_LIBRARY}") # Print the found library

    # Check if the library is static or dynamic
    get_filename_component(GMP_LIBRARY_EXTENSION "${GMP_LIBRARY}" EXT)

    if(GMP_LIBRARY_EXTENSION MATCHES ".a|.lib")
        # Static library
        set(GMP_LIBRARY_TYPE STATIC)
    elseif(GMP_LIBRARY_EXTENSION MATCHES ".so|.dll")
        # Dynamic library
        set(GMP_LIBRARY_TYPE SHARED)
    else()
        message(WARNING "Unknown library type for GMP: ${GMP_LIBRARY}")
        set(GMP_LIBRARY_TYPE UNKNOWN)
    endif()

    add_library(GMP::GMP ${GMP_LIBRARY_TYPE} IMPORTED) # Create an imported target

    # Set the include directories and library location
    set_target_properties(GMP::GMP PROPERTIES
        INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
        IMPORTED_LOCATION "${GMP_LIBRARY}"
    )
else()
    message(FATAL_ERROR "GMP not found. Please install GMP or set the GMP_DIR variable.") # Error if not found
endif()
