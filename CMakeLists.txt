cmake_minimum_required(VERSION 3.0)

project(QuadraticDifferential C CXX Fortran)
set (QuadraticDifferential_VERSION_MAJOR 0)
set (QuadraticDifferential_VERSION_MINOR 1)

set(PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${PROJECT_ROOT}/bin CACHE PATH "" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(BLAS_LIBRARY_PATH "" CACHE PATH "")
set(CMAKE_LIBRARY_PATH "${CMAKE_LIBRARY_PATH} ${BLAS_LIBRARY_PATH}")
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
message("BLAS LIBRARIES ${BLAS_LIBRARIES}")

# use openmp
OPTION (USE_OpenMP "Use OpenMP" OFF)
if(USE_OpenMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# set vc runtime library as /MT
if(${MSVC})
    foreach(flag
            CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
            CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    if(${flag} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${flag} "${${flag}}")
    endif(${flag} MATCHES "/MD")
    endforeach(flag)
    # pass to external projects
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_RELEASE} /MTd")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
    set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} /MT")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /MT")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /libs:static")
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()
if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgfortran")
endif()

include_directories(third-party/MeshLib/core)
include_directories(third-party/eigen)

set(ARPACK_ROOT "${PROJECT_ROOT}/third-party/arpack-ng/")
set(ARPACK_BUILD_DIR "${PROJECT_ROOT}/build/arpack-ng/")
include(ExternalProject)
ExternalProject_Add(arpack
    PREFIX "${ARPACK_ROOT}"
    BINARY_DIR "${ARPACK_BUILD_DIR}"
    STAMP_DIR  "${ARPACK_BUILD_DIR}"
    SOURCE_DIR "${ARPACK_ROOT}"
    INSTALL_DIR "${ARPACK_ROOT}"
    CMAKE_CACHE_ARGS "-DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}" "-DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}"
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${ARPACK_ROOT} -DBUILD_SHARED_LIBS=OFF
        -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
        #-DCMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH}
        #-DBLA_STATIC=ON
)
add_library(libarpack STATIC IMPORTED)
add_dependencies(libarpack arpack)

if(WIN32)
    set_property(TARGET libarpack PROPERTY IMPORTED_LOCATION ${ARPACK_ROOT}/lib/arpack.lib)
endif()
if(UNIX)
    set_property(TARGET libarpack PROPERTY IMPORTED_LOCATION ${ARPACK_ROOT}/lib/libarpack.a)
endif()
#ExternalProject_Get_Property(arpack OUTPUT_NAME)

link_directories(${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})

add_executable(QuadraticDifferential
    source/main.cpp
    source/QuadraticDifferential.cpp
)
set_property(TARGET QuadraticDifferential PROPERTY CXX_STANDARD 14)
set_property(TARGET QuadraticDifferential PROPERTY CXX_STANDARD_REQUIRED ON)

target_link_libraries(QuadraticDifferential libarpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

install (TARGETS QuadraticDifferential DESTINATION ${PROJECT_ROOT}/bin)
