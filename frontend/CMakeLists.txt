
#project(frontend)


file(GLOB_RECURSE SRCS *.cpp)

include_directories(${CMAKE_SOURCE_DIR})



add_executable(frontend_libOTe  ${SRCS})

# default location for NTL
set(NTL_Dirs "${CMAKE_SOURCE_DIR}/thirdparty/linux/ntl/")
include_directories("${NTL_Dirs}/include/")
link_directories("${NTL_Dirs}/src/")

#NTL
set(NTL_Lib_Dirs "${NTL_Dirs}/src/")
set(NTL_Include_Dirs "${NTL_Dirs}/include/")
set(NTL_LIB "${NTL_Lib_Dirs}libntl.a")
if(NOT EXISTS "${NTL_LIB}")
	# we failed to find it, try using find_library
	find_library(
			NTL_LIB
			NAMES ntl
			HINTS NTL_Lib_Dirs)

	# if we still cant fint it, throw an error
	if(NOT NTL_LIB)
		Message(${NTL_LIB})
		message(FATAL_ERROR "Failed to find NTL at ${NTL_Lib_Dirs} or at system locations." )
	endif()
else()
	include_directories(${NTL_Include_Dirs})
	#message(" ----  NTL -----  ${NTL_Include_Dirs}")
endif()


find_package(xxHash CONFIG REQUIRED)

target_link_libraries(frontend_libOTe PRIVATE xxHash::xxhash libOTe_Tests tests_cryptoTools crypto.a dl ssl gmpxx gmp pthread ${NTL_LIB} linbox iml boost_system blas givaro gomp)

