
# avoid having empty buildtype
set(CMAKE_BUILD_TYPE_INIT Release)

set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS
	$<$<CONFIG:Debug>:DEBUG;_DEBUG>
	$<$<CONFIG:Release>:NDEBUG>
	$<$<CONFIG:MinSizeRel>:NDEBUG>
	$<$<CONFIG:RelWithDebInfo>:NDEBUG>
)

cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0005 NEW)

cmake_minimum_required(VERSION 2.8)

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
		COMMAND python3 -c "from distutils.sysconfig import PREFIX; print(PREFIX)"
		OUTPUT_VARIABLE _prefix
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	set(PYTHON_ROOT_DIR "${_prefix}" CACHE PATH INTERNAL)
	unset(_prefix)
endif()

if(NOT DEFINED PYTHON_SITE_PACKAGES)
	execute_process(
		COMMAND python3 -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
		OUTPUT_VARIABLE _site_packages
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	set(PYTHON_SITE_PACKAGES "${_site_packages}" CACHE PATH INTERNAL)
	unset(_site_packages)
endif()

if(NOT DEFINED PYTHON_INCLUDE_DIR)
	execute_process(
		COMMAND python3 -c "from distutils.sysconfig import get_python_inc; print(get_python_inc(plat_specific=1))"
		OUTPUT_VARIABLE _inc_dir
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	set(PYTHON_INCLUDE_DIR "${_inc_dir}" CACHE PATH INTERNAL)
	unset(_inc_dir)
endif()

# set standard cmake/python names, uncached
set(PYTHON_INCLUDE_DIRS "${PYTHON_INCLUDE_DIR}")


# -----------------------------------------------------------------------------
# curve_fit_nd (C)

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

include_directories(
	SYSTEM
	${PYTHON_INCLUDE_DIRS}
)


add_library(curve_fit_nd_lib ${SRC})


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
