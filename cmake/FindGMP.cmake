# cmake/FindGMP.cmake
find_path(GMP_INCLUDE_DIR gmp.h
    PATHS
    /usr/include
    /mingw64/include
)

find_library(GMP_LIBRARY NAMES gmp
    PATHS
    /usr/lib
    /mingw64/lib
)

find_library(GMPXX_LIBRARY NAMES gmpxx
    PATHS
    /usr/lib
    /mingw64/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG GMP_INCLUDE_DIR GMP_LIBRARY)

if(GMP_FOUND)
    set(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})
    set(GMP_LIBRARIES ${GMP_LIBRARY} ${GMPXX_LIBRARY})
    add_library(GMP::GMP STATIC IMPORTED)
    set_target_properties(GMP::GMP PROPERTIES IMPORTED_LOCATION "${GMP_LIBRARY}")
endif()
