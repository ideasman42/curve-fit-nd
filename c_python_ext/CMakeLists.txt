cmake_minimum_required(VERSION 3.10)

project(curve_fit_nd_ext C)

option(WITH_DEBUG_TIME "Print timing information for curve fitting" OFF)
set(WITH_PORT "C" CACHE STRING "Implementation port to use (C, RUST)")

# avoid having empty buildtype
set(CMAKE_BUILD_TYPE_INIT Release)

set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS
	$<$<CONFIG:Debug>:DEBUG;_DEBUG>
	$<$<CONFIG:Release>:NDEBUG>
	$<$<CONFIG:MinSizeRel>:NDEBUG>
	$<$<CONFIG:RelWithDebInfo>:NDEBUG>
)

if(WITH_DEBUG_TIME)
	add_definitions(-DDEBUG_TIME)
endif()

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin CACHE INTERNAL "" FORCE )
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib CACHE INTERNAL "" FORCE )

# needed for python lib
if(CMAKE_COMPILER_IS_GNUCC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -funsigned-char")
endif()

if(MSVC)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /J")
endif()

# -----------------------------------------------------------------------------
# configure python

if(NOT DEFINED PYTHON_ROOT_DIR)
	execute_process(
		COMMAND python3 -c "import sys; print(sys.prefix)"
		OUTPUT_VARIABLE _prefix
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	set(PYTHON_ROOT_DIR "${_prefix}" CACHE PATH INTERNAL)
	unset(_prefix)
endif()

if(NOT DEFINED PYTHON_SITE_PACKAGES)
	execute_process(
		COMMAND python3 -c "import sysconfig; print(sysconfig.get_path('purelib'))"
		OUTPUT_VARIABLE _site_packages
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	set(PYTHON_SITE_PACKAGES "${_site_packages}" CACHE PATH INTERNAL)
	unset(_site_packages)
endif()

if(NOT DEFINED PYTHON_INCLUDE_DIR)
	execute_process(
		COMMAND python3 -c "import sysconfig; print(sysconfig.get_path('include'))"
		OUTPUT_VARIABLE _inc_dir
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	set(PYTHON_INCLUDE_DIR "${_inc_dir}" CACHE PATH INTERNAL)
	unset(_inc_dir)
endif()

# set standard cmake/python names, uncached
set(PYTHON_INCLUDE_DIRS "${PYTHON_INCLUDE_DIR}")


# -----------------------------------------------------------------------------
# curve_fit_nd library (C or Rust)

include_directories(
	SYSTEM
	${PYTHON_INCLUDE_DIRS}
)

if(WITH_PORT STREQUAL "RUST")
	# Rust implementation
	set(RUST_DIR "${CMAKE_SOURCE_DIR}/../rust")

	# Maximum supported dimensions (stack-allocated array size).
	set(RUST_DIMS_MAX "32" CACHE STRING "Maximum supported dimensions for curve fitting")
	if(NOT RUST_DIMS_MAX GREATER 0)
		message(FATAL_ERROR "RUST_DIMS_MAX must be greater than 0")
	endif()

	# Select debug or release build based on CMAKE_BUILD_TYPE
	if(CMAKE_BUILD_TYPE STREQUAL "Debug")
		set(RUST_BUILD_TYPE "debug")
		set(RUST_BUILD_FLAG "")
	else()
		set(RUST_BUILD_TYPE "release")
		set(RUST_BUILD_FLAG "--release")
	endif()

	set(RUST_LIB "${RUST_DIR}/target/${RUST_BUILD_TYPE}/libcurve_fit_nd.a")

	# Check that cargo is available
	find_program(CARGO_EXECUTABLE cargo)
	if(NOT CARGO_EXECUTABLE)
		message(FATAL_ERROR "cargo not found. Please install Rust: https://rustup.rs/")
	endif()

	# Rust source files for dependency tracking
	set(RUST_SOURCES
		${RUST_DIR}/src/lib.rs
		${RUST_DIR}/src/mod.rs
		${RUST_DIR}/src/curve_fit_cubic_refit.rs
		${RUST_DIR}/src/curve_fit_from_polys.rs
		${RUST_DIR}/src/curve_fit_cubic.rs
		${RUST_DIR}/src/math_vector.rs
		${RUST_DIR}/src/min_heap.rs
		${RUST_DIR}/Cargo.toml
		${RUST_DIR}/build.rs
	)

	# Build the Rust library
	# Disable FMA (fused multiply-add) to match C/GCC floating-point behavior.
	# Without this, LLVM's FMA optimizations cause tiny differences that compound
	# through Newton-Raphson iterations, leading to different split decisions.
	add_custom_command(
		OUTPUT ${RUST_LIB}
		COMMAND ${CMAKE_COMMAND} -E env DIMS_MAX=${RUST_DIMS_MAX} RUSTFLAGS=-Ctarget-feature=-fma cargo build ${RUST_BUILD_FLAG}
		WORKING_DIRECTORY ${RUST_DIR}
		DEPENDS ${RUST_SOURCES}
		COMMENT "Building Rust curve_fit_nd library (${RUST_BUILD_TYPE}, dims_max=${RUST_DIMS_MAX})"
	)

	add_custom_target(curve_fit_nd_rust DEPENDS ${RUST_LIB})

	# Import the pre-built Rust static library
	add_library(curve_fit_nd_lib STATIC IMPORTED GLOBAL)
	set_target_properties(curve_fit_nd_lib PROPERTIES
		IMPORTED_LOCATION ${RUST_LIB}
	)
	add_dependencies(curve_fit_nd_lib curve_fit_nd_rust)

	# Clean target for Rust build artifacts
	add_custom_target(clean_rust
		COMMAND cargo clean
		WORKING_DIRECTORY ${RUST_DIR}
		COMMENT "Cleaning Rust build artifacts"
	)

	# Pass DIMS_MAX to C code for runtime validation
	add_definitions(-DWITH_RUST_PORT -DRUST_DIMS_MAX=${RUST_DIMS_MAX})
else()
	# C implementation (default)
	set(SRC
		../c/intern/curve_fit_cubic.c
		../c/intern/curve_fit_cubic_refit.c

		../c/curve_fit_nd.h
		../c/intern/curve_fit_inline.h

		# generic helpers
		../c/intern/generic_heap.c

		../c/intern/generic_alloc_impl.h
		../c/intern/generic_heap.h
	)

	add_library(curve_fit_nd_lib ${SRC})
endif()


# -----------------------------------------------------------------------------
# curve_fit_nd_ext (C/Python API)

set(SRC
	curve_fit_nd_ext.c
)

include_directories(
	../c
)


# -----------------------------------------------------------------------------
# python module

add_library(curve_fit_nd_module SHARED ${SRC})

set_target_properties(
	curve_fit_nd_module
	PROPERTIES
		PREFIX ""
		OUTPUT_NAME curve_fit_nd
		LINKER_LANGUAGE C
		LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
		RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin  # only needed on windows
)

if(WIN32)
	# python modules use this
	set_target_properties(
		curve_fit_nd_module
		PROPERTIES
		SUFFIX ".pyd"
	)	
endif()

# link
target_link_libraries(curve_fit_nd_module
	curve_fit_nd_lib
)

# -----------------------------------------------------------------------------
# installation

install(
	TARGETS curve_fit_nd_module
	LIBRARY
	DESTINATION
		"${PYTHON_SITE_PACKAGES}"
)
