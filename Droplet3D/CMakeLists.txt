cmake_minimum_required (VERSION 2.6)

cmake_policy(SET CMP0042 OLD)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

project (Droplet3D)
set (MULTITRACKER_MAJOR_VERSION "1")
set (MULTITRACKER_MINOR_VERSION "0")
set (MULTITRACKER_BUGFIX_VERSION "0")
set (MULTITRACKER_VERSION_NUMBER "${MULTITRACKER_MAJOR_VERSION}.${MULTITRACKER_MINOR_VERSION}.${MULTITRACKER_BUGFIX_VERSION}")
set (MULTITRACKER_VERSION "${MULTITRACKER_VERSION_NUMBER}")
add_definitions (-DMULTITRACKER_VERSION="${MULTITRACKER_VERSION}")

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# add directory with macros
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# add cmake file
include (common)

# set options

if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions (-DDEBUG)
endif (CMAKE_BUILD_TYPE MATCHES Debug)

# Eigen library is required
find_package (Eigen REQUIRED)
include_directories (${EIGEN_INCLUDE_DIR})

option (USE_PNG "Builds in support for exporting to png" ON)
if (USE_PNG)
  find_package (PNG)
  if (PNG_FOUND)
    add_definitions (-DHAVE_PNG)
    add_definitions (${PNG_DEFINITIONS})
    include_directories (${PNG_INCLUDE_DIR})
    set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${PNG_LIBRARIES})
  endif (PNG_FOUND)
endif (USE_PNG)

find_package (LAPACK REQUIRED)
add_definitions (-DHAVE_LAPACK)
set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${LAPACK_LIBRARIES})

find_package (BLAS REQUIRED)
add_definitions (-DHAVE_BLAS)
set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${BLAS_LIBRARIES})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/LosTopos3D)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common/newsparse)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LosTopos/common/tunicate)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Apps/Droplet3D/fmmtl)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Apps/Droplet3D)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory (LosTopos)

add_subdirectory (Apps)

