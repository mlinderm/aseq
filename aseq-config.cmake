######
# This is the find-script for your project/library used by CMake if this project/libary is to be
# included in another project/library.
######

set(ASEQ_INCLUDES "")
set(ASEQ_LIBRARIES "")

# Definition of function "find" with two mandatory arguments, "LIB_NAME" and "HEADER".
macro (find LIB_NAME HEADER)

	set(HINT_PATHS ${ARGN})

	if (${LIB_NAME} STREQUAL "aseq")
		set(LIB_NAME_UPPER ASEQ)
		set(LIBNAME aseq)
	else()
		string(TOUPPER ASEQ_${LIB_NAME} LIB_NAME_UPPER)
		set(LIBNAME ${LIB_NAME})
	endif()

	find_path(
		${LIB_NAME_UPPER}_INCLUDE_DIR
		${HEADER}
		${ENV_ASEQ_DIR}/include
		${ENV_ASEQ_DIR}/source/${LIB_NAME}/include
		${ASEQ_DIR}/include
		${ASEQ_DIR}/source/${LIB_NAME}/include
		${ENV_PROGRAMFILES}/aseq/include
		/usr/include
		/usr/local/include
		/sw/include
		/opt/local/include
		DOC "The directory where ${HEADER} resides"
	)


	find_library(
		${LIB_NAME_UPPER}_LIBRARY_RELEASE
		NAMES ${LIBNAME}
		PATHS ${HINT_PATHS}
		DOC "The ${LIB_NAME} library"
	)
	find_library(
		${LIB_NAME_UPPER}_LIBRARY_DEBUG
		NAMES ${LIBNAME}d
		PATHS ${HINT_PATHS}
		DOC "The ${LIB_NAME} debug library"
	)


	if(${LIB_NAME_UPPER}_LIBRARY_RELEASE AND ${LIB_NAME_UPPER}_LIBRARY_DEBUG)
		set(${LIB_NAME_UPPER}_LIBRARY "optimized" ${${LIB_NAME_UPPER}_LIBRARY_RELEASE} "debug" ${${LIB_NAME_UPPER}_LIBRARY_DEBUG})
	elseif(${LIB_NAME_UPPER}_LIBRARY_RELEASE)
		set(${LIB_NAME_UPPER}_LIBRARY ${${LIB_NAME_UPPER}_LIBRARY_RELEASE})
	elseif(${LIB_NAME_UPPER}_LIBRARY_DEBUG)
		set(${LIB_NAME_UPPER}_LIBRARY ${${LIB_NAME_UPPER}_LIBRARY_DEBUG})
	endif()

	list(APPEND ASEQ_INCLUDES ${${LIB_NAME_UPPER}_INCLUDE_DIR})
	list(APPEND ASEQ_LIBRARIES ${${LIB_NAME_UPPER}_LIBRARY})

	# DEBUG MESSAGES
	# message("${LIB_NAME_UPPER}_INCLUDE_DIR     = ${${LIB_NAME_UPPER}_INCLUDE_DIR}")
	# message("${LIB_NAME_UPPER}_LIBRARY_RELEASE = ${${LIB_NAME_UPPER}_LIBRARY_RELEASE}")
	# message("${LIB_NAME_UPPER}_LIBRARY_DEBUG   = ${${LIB_NAME_UPPER}_LIBRARY_DEBUG}")
	# message("${LIB_NAME_UPPER}_LIBRARY         = ${${LIB_NAME_UPPER}_LIBRARY}")

endmacro(find)


# load standard CMake arguments (c.f. http://stackoverflow.com/questions/7005782/cmake-include-findpackagehandlestandardargs-cmake)
include(FindPackageHandleStandardArgs)

if(CMAKE_CURRENT_LIST_FILE)
	get_filename_component(ASEQ_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
endif()

file(TO_CMAKE_PATH "$ENV{PROGRAMFILES}" ENV_PROGRAMFILES)
file(TO_CMAKE_PATH "$ENV{ASEQ_DIR}" ENV_ASEQ_DIR)

set(LIB_PATHS
	${ASEQ_DIR}/build
	${ASEQ_DIR}/build/Release
	${ASEQ_DIR}/build/Debug
	${ASEQ_DIR}/build-release
	${ASEQ_DIR}/build-debug
	${ASEQ_DIR}/lib
	${ENV_ASEQ_DIR}/lib
	${ENV_PROGRAMFILES}/aseq/lib
	/usr/lib
	/usr/local/lib
	/sw/lib
	/opt/local/lib
	/usr/lib64
	/usr/local/lib64
	/sw/lib64
	/opt/local/lib64
)

# Find libraries
find(fiblib fiblib/fiblib_api.h ${LIB_PATHS})

if(ASEQ_FIBLIB_LIBRARY)
	# add dependencies
endif()


# DEBUG
# message("ASEQ_INCLUDES  = ${ASEQ_INCLUDES}")
# message("ASEQ_LIBRARIES = ${ASEQ_LIBRARIES}")

find_package_handle_standard_args(ASEQ DEFAULT_MSG ASEQ_LIBRARIES ASEQ_INCLUDES)
mark_as_advanced(ASEQ_FOUND)
