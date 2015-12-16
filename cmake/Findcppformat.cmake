# Find the cppformat include directory
# The following variables are set if cppformat is found.
#  cppformat_INCLUDE_DIR
#  cppformat_LIBRARY

include(FindPackageHandleStandardArgs)

find_path(cppformat_INCLUDE_DIR
    NAMES cppformat/format.h
    DOC "cppformat library header files")

find_library(cppformat_LIBRARY
    NAMES cppformat
    PATH_SUFFIXES lib lib64)

find_package_handle_standard_args(cppformat DEFAULT_MSG
    cppformat_INCLUDE_DIR cppformat_LIBRARY)
