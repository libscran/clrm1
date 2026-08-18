cmake_minimum_required(VERSION 3.24)

project(
    clrm1
    VERSION 0.2.1
    DESCRIPTION "Centered log-ratio minus 1"
    LANGUAGES CXX
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Library
add_library(clrm1 INTERFACE)
add_library(libscran::clrm1 ALIAS clrm1)

target_include_directories(clrm1 INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/clrm1>)
target_compile_features(clrm1 INTERFACE cxx_std_17)

# Dependencies
option(CLRM1_FETCH_EXTERN "Automatically fetch clrm1's external dependencies." ON)
if(CLRM1_FETCH_EXTERN)
    add_subdirectory(extern)
else()
    find_package(tatami_tatami 4.1.0 CONFIG REQUIRED)
    find_package(tatami_tatami_stats 3.1.0 CONFIG REQUIRED)
endif()

target_link_libraries(clrm1 INTERFACE tatami::tatami tatami::tatami_stats)

install(
    DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clrm1
)

install(
    TARGETS clrm1
    EXPORT clrm1Targets
)

install(
    EXPORT clrm1Targets
    FILE libscran_clrm1Targets.cmake
    NAMESPACE libscran::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_clrm1
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_clrm1Config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_clrm1
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_clrm1ConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion
)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/libscran_clrm1Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_clrm1ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_clrm1
)
