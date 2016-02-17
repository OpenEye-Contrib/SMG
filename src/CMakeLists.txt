cmake_minimum_required(VERSION 2.8)

project(SMG)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build_utils/CMakeModules/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -W -Wall  -Wshadow -Wno-error=shadow -Wno-error=unused-parameter -Wno-error=unused-function -fomit-frame-pointer -ffast-math -fvisibility=hidden -fvisibility-inlines-hidden -fPIC  -Werror " )

if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^3.1" )
  message( "Hopefully, the Ubuntu 14.04 machine, with a modern compiler." )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
endif()
if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^2.6." )
  message( "It looks like a trusty old Centos 6 machine with an ancient compiler." )
endif()

find_package(OEToolkits COMPONENTS oedepict oechem oesystem oeplatform)

set(SMG_SRCS ${SMG_SOURCE_DIR}/smg.cc
${SMG_SOURCE_DIR}/spiv_nogr_bits.cc
${SMG_SOURCE_DIR}/SpivMolecule.cc)

set(SMG_INCS
${SMG_SOURCE_DIR}/SpivMolecule.H
${SMG_SOURCE_DIR}/spiv_nogr_bits.H)

set(SMG_DACLIB_SRCS
${SMG_SOURCE_DIR}/apply_daylight_arom_model_to_oemol.cc
${SMG_SOURCE_DIR}/build_time.cc
${SMG_SOURCE_DIR}/superfast_hash.cc
${SMG_SOURCE_DIR}/MurmurHash2.cc
${SMG_SOURCE_DIR}/read_smarts_file.cc
${SMG_SOURCE_DIR}/PharmPoint.cc)

set(SMG_DACLIB_INCS
${SMG_SOURCE_DIR}/crash.H
${SMG_SOURCE_DIR}/stddefs.H
${SMG_SOURCE_DIR}/FileExceptions.H
${SMG_SOURCE_DIR}/PharmPoint.H
${SMG_SOURCE_DIR}/SMARTSExceptions.H)

include_directories( SYSTEM ${OEToolkits_INCLUDE_DIRS})

set(EXECUTABLE_OUTPUT_PATH ${SMG_SOURCE_DIR}/exe_${CMAKE_BUILD_TYPE})

set(SMG_LIBS ${LIBS}
  ${OEToolkits_LIBRARIES}
  ${Boost_LIBRARIES})

add_executable(smg ${SMG_SRCS} ${SMG_DACLIB_SRCS}
  ${SMG_INCS}  ${SMG_DACLIB_INCS})
target_link_libraries(smg z ${SMG_LIBS} z pthread rt)
