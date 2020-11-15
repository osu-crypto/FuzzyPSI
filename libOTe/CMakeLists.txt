enable_language(ASM)

file(GLOB_RECURSE SRCS *.cpp)

set(SRCS "${SRCS}")

add_library(libOTe STATIC ${SRCS})

target_include_directories(libOTe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..) 
target_link_libraries(libOTe cryptoTools)

if(ENABLE_SSE)
    target_compile_options(libOTe PRIVATE -maes -msse2 -msse3 -msse4.1 -mpclmul)
endif()

if(ENABLE_SIMPLESTOT_ASM)
	target_link_libraries(libOTe SimplestOT)
    #target_compile_options(libOTe PRIVATE -fPIC -no-pie)
endif()

if(ENABLE_MR_KYBER)
	target_link_libraries(libOTe KyberOT)
endif(ENABLE_MR_KYBER)



#  coproto
############################################


if(NOT DEFINED COPROTO_DIR)
set(COPROTO_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../coproto")
endif()
if(NOT DEFINED COPROTO_INCLUDE_DIRS)
    set(COPROTO_INCLUDE_DIRS "${COPROTO_DIR}")
endif()


find_library(COPROTO_LIB NAMES coproto  HINTS "${COPROTO_DIR}/lib")

# if we cant fint it, throw an error
if(NOT COPROTO_LIB)
    message(FATAL_ERROR "Failed to find coproto at ${COPROTO_DIR}. You can specify its location via cmake with -DCOPROTO_DIR=path" )
endif()

message(STATUS "COPROTO_LIB:  ${COPROTO_LIB}")
message(STATUS "COPROTO_INC:  ${COPROTO_INCLUDE_DIRS}\n")


target_include_directories(libOTe PUBLIC ${COPROTO_INCLUDE_DIRS}) 
target_link_libraries(libOTe ${COPROTO_LIB})



#############################################
#            Install                        #
#############################################

# install library
install(TARGETS libOTe DESTINATION lib)

# install headers
install(DIRECTORY . DESTINATION include/libOTe FILES_MATCHING PATTERN "*.h")
