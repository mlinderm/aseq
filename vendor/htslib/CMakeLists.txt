set(target htslib)

set(api_includes
    htslib/bgzf.h
    htslib/faidx.h
    htslib/hfile.h
    htslib/hts.h
    htslib/hts_defs.h
    htslib/kfunc.h
    htslib/khash.h
    htslib/khash_str2int.h
    htslib/klist.h
    htslib/knetfile.h
    htslib/kseq.h
    htslib/ksort.h
    htslib/kstring.h
    htslib/regidx.h
    htslib/sam.h
    htslib/synced_bcf_reader.h
    htslib/tbx.h
    htslib/vcf.h
    htslib/vcf_sweep.h
    htslib/vcfutils.h
    )

set(sources
    kfunc.c
    knetfile.c
    kstring.c
    bgzf.c
    faidx.c
    hfile.c
    hfile_net.c
    hts.c
    md5.c
    regidx.c
    sam.c
    synced_bcf_reader.c
    vcf_sweep.c
    tbx.c
    vcf.c
    vcfutils.c
    cram/cram_codecs.c
    cram/cram_decode.c
    cram/cram_encode.c
    cram/cram_external.c
    cram/cram_index.c
    cram/cram_io.c
    cram/cram_samtools.c
    cram/cram_stats.c
    cram/files.c
    cram/mFILE.c
    cram/open_trace_file.c
    cram/pooled_alloc.c
    cram/rANS_static.c
    cram/sam_header.c
    cram/string_alloc.c
    cram/thread_pool.c
    cram/vlen.c
    cram/zfio.c
    )

configure_file(version.h.in ${CMAKE_CURRENT_BINARY_DIR}/version.h)
configure_file(config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)

find_package(ZLIB REQUIRED)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${ZLIB_INCLUDE_DIRS}
    )

add_definitions(-D_FILE_OFFSET_BITS=64)
add_library(${target} STATIC ${api_includes} ${sources})
target_compile_options(${target} PUBLIC -fPIC)
target_link_libraries(${target} ${ZLIB_LIBRARIES} )
set_target_properties(${target}
    PROPERTIES
    OUTPUT_NAME         "${target}"
    CLEAN_DIRECT_OUTPUT 1
    )

install(TARGETS ${target}
    RUNTIME DESTINATION ${INSTALL_BIN}
    LIBRARY DESTINATION ${INSTALL_SHARED}
    ARCHIVE DESTINATION ${INSTALL_LIB}
    )
