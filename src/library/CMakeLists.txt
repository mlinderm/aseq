
set(target library)
message(STATUS "Lib ${target}")


# External libraries
set(HTSLIB_DIR "${VENDOR_SOURCE_DIR}/htslib")
add_subdirectory(${HTSLIB_DIR} "${VENDOR_BINARY_DIR}/htslib")

# e.g., find_package(THIRDPARTY REQUIRED)


# Includes

include_directories(
    ${HTSLIB_DIR}
    # e.g., ${TRDP_INCLUDE_DIR}
)

include_directories(
    BEFORE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)


# Libraries

set(libs
    glog
    cppformat
    htslib
    # e.g., ${THIRDPARTY_LIBRARY}
)


# Compiler definitions

if (OPTION_BUILD_STATIC)
    add_definitions("-D${META_PROJECT_NAME_UPPER}_STATIC")
else()
    add_definitions("-DFIBLIB_EXPORTS")
endif()


# Sources

set(include_path "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(source_path "${CMAKE_CURRENT_SOURCE_DIR}/src")

set(api_includes
    include/aseq/util/exception.hpp
    include/aseq/util/flyweight.hpp
    include/aseq/io/line.hpp
    include/aseq/io/variant.hpp
    include/aseq/io/file_format.hpp
    include/aseq/io/vcf.hpp
    include/aseq/io/reference.hpp
    include/aseq/model/contig.hpp
    include/aseq/model/allele.hpp
    include/aseq/model/variant_context.hpp
    include/aseq/algorithm/variant.hpp
    include/aseq/util/any.hpp
    include/aseq/util/attributes.hpp
)

set(sources
    src/io/line_reader.cpp
    src/io/line_writer.cpp
    src/io/variant_source.cpp
    src/io/variant_sink.cpp
    src/io/vcf.cpp
    src/io/vcf_source.cpp
    src/io/vcf_sink.cpp
    src/io/reference_source.cpp
    include/aseq/io/reference-mock.hpp
    src/model/allele.cpp
    src/model/variant_context.cpp
    src/algorithm/variant_consensus.cpp
)

# Group src files
set(header_group "Header Files (API)")
set(source_group "Source Files")
source_group_by_path(${include_path} "\\\\.h$|\\\\.hpp$" ${header_group} ${api_includes})
source_group_by_path(${source_path} "\\\\.cpp$|\\\\.c$|\\\\.h$|\\\\.hpp$" ${source_group} ${sources})


# Build library
add_library(${target} ${api_includes} ${sources})

target_link_libraries(${target} ${Boost_LIBRARIES} ${libs})

target_compile_options(${target} PRIVATE ${DEFAULT_COMPILE_FLAGS})

set_target_properties(${target}
    PROPERTIES
    LINKER_LANGUAGE             CXX
    CXX_STANDARD                14
    FOLDER                      "${IDE_FOLDER}"
    COMPILE_DEFINITIONS_DEBUG   "${DEFAULT_COMPILE_DEFS_DEBUG}"
    COMPILE_DEFINITIONS_RELEASE "${DEFAULT_COMPILE_DEFS_RELEASE}"
    LINK_FLAGS_DEBUG            "${DEFAULT_LINKER_FLAGS_DEBUG}"
    LINK_FLAGS_RELEASE          "${DEFAULT_LINKER_FLAGS_RELEASE}"
    DEBUG_POSTFIX               "d${DEBUG_POSTFIX}"
    INCLUDE_PATH	            ${include_path}
    OUTPUT_NAME                 "${PROJECT_NAME}"
)


# Deployment

# Library
install(TARGETS ${target}
    RUNTIME DESTINATION ${INSTALL_BIN}
    LIBRARY DESTINATION ${INSTALL_SHARED}
    ARCHIVE DESTINATION ${INSTALL_LIB}
)

# Header files
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/${target} DESTINATION ${INSTALL_INCLUDE})
