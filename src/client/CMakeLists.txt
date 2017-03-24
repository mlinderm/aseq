
set(target client)
message(STATUS "App ${target}")


# External libraries
set(DOCOPT_DIR "${VENDOR_SOURCE_DIR}/docopt")
add_subdirectory(${DOCOPT_DIR} "${VENDOR_BINARY_DIR}/docopt")

# Includes

include_directories(
  "${DOCOPT_DIR}"
)

include_directories(
    BEFORE
    ${CMAKE_SOURCE_DIR}/src/library/include
)


# Libraries

set(libs
  glog
  docopt
  library
)


# Compiler definitions

if (OPTION_BUILD_STATIC)
    add_definitions("-D${META_PROJECT_NAME_UPPER}_STATIC")
endif()


# Sources

set(sources
    main.cpp
    commands.hpp
    variants.cpp)


# Build executable

add_executable(${target} ${sources})

target_link_libraries(${target} ${libs})

target_compile_options(${target} PRIVATE ${DEFAULT_COMPILE_FLAGS})

set_target_properties(${target}
    PROPERTIES
    LINKER_LANGUAGE             CXX
    CXX_STANDARD                14
    CXX_STANDARD_REQUIRED       ON
    FOLDER                      "${IDE_FOLDER}"
    COMPILE_DEFINITIONS_DEBUG   "${DEFAULT_COMPILE_DEFS_DEBUG}"
    COMPILE_DEFINITIONS_RELEASE "${DEFAULT_COMPILE_DEFS_RELEASE}"
    LINK_FLAGS                  "${DEFAULT_LINKER_FLAGS}"
    LINK_FLAGS_DEBUG            "${DEFAULT_LINKER_FLAGS_DEBUG}"
    LINK_FLAGS_RELEASE          "${DEFAULT_LINKER_FLAGS_RELEASE}"
    DEBUG_POSTFIX               "d${DEBUG_POSTFIX}"
    VERSION                     "${META_VERSION}"
    OUTPUT_NAME                 "${PROJECT_NAME}"
    CLEAN_DIRECT_OUTPUT         1
)


# Deployment

install(TARGETS ${target}
    RUNTIME DESTINATION ${INSTALL_BIN}
)
