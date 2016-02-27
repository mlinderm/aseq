# External libraries
find_package(Boost 1.59.0 REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

include_directories(${GLOG_INCLUDE_DIR})

set(CPPFORMAT_DIR "${VENDOR_SOURCE_DIR}/cppformat")
add_subdirectory(${CPPFORMAT_DIR} "${VENDOR_BINARY_DIR}/cppformat")
include_directories(${CPPFORMAT_DIR})

# Generate meta-header
configure_file(aseq-version.h.in ${CMAKE_CURRENT_BINARY_DIR}/aseq-version.h)

# We do NOT want to deploy this header (at least not for binary packages)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/template/template-version.h DESTINATION ${INSTALL_INCLUDE}/template)

# Include directories
include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}
)

# Libraries
set(IDE_FOLDER "")
add_subdirectory(library)

# Applications
set(IDE_FOLDER "Applications")
add_subdirectory(client)

# Tests
set(IDE_FOLDER "Tests")
add_subdirectory(tests)
