message(STATUS "Configuring for platform Linux/GCC.")

execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
    OUTPUT_VARIABLE GCC_VERSION)

set(LINUX_COMPILE_DEFS
    LINUX                     # Linux system
    PIC                       # Position-independent code
    _REENTRANT                # Reentrant code
)
set(DEFAULT_COMPILE_DEFS_DEBUG
    ${LINUX_COMPILE_DEFS}
    _DEBUG                    # Debug build
)
set(DEFAULT_COMPILE_DEFS_RELEASE
    ${LINUX_COMPILE_DEFS}
    NDEBUG                    # Release build
)

set(LINUX_COMPILE_FLAGS
      -pthread      # -> use pthread library
    # -no-rtti      # -> disable c++ rtti
      -pipe
      -Wall
      -Wextra
      -fPIC
      -Wreturn-type 
      -Wfloat-equal 
      -Wno-shadow
      -Wcast-align 
      -Wconversion
      -Wno-unused-parameter
)

set(DEFAULT_COMPILE_FLAGS
    ${LINUX_COMPILE_FLAGS}
    $<$<CONFIG:Debug>:   
    >
    $<$<CONFIG:Release>: 
    >
)

set(LINUX_LINKER_FLAGS 
	-pthread
)

set(LINUX_LINKER_FLAGS "-pthread")

set(DEFAULT_LINKER_FLAGS_RELEASE ${LINUX_LINKER_FLAGS})
set(DEFAULT_LINKER_FLAGS_DEBUG ${LINUX_LINKER_FLAGS})
set(DEFAULT_LINKER_FLAGS ${LINUX_LINKER_FLAGS})
