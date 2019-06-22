# Find a Intel® Math Kernel Library (Intel® MKL) installation and provide
# all necessary variables and macros to compile software for it.
#
# MKLROOT is required in your system
#
# we use mkl_link_tool to get the library needed depending on variables
#
# The following are set after the configuration is done:
#  MKL_FOUND        -  system has MKL
#  MKL_ROOT_DIR     -  path to the MKL base directory
#  MKL_INCLUDE      -  the MKL include directory
#  MKL_LIBRARIES    -  MKL libraries
#
#
#
# Sample usage:
#    find_package(MKL REQUIRED)
#    if (MKL_FOUND)
#        include_directories(${MKL_INCLUDE})
#        # and for each of your dependent executable/library targets:
#        target_link_libraries(<YourTarget> ${MKL_LIBRARIES})
#    endif()
#
#
# AUTHOR
# Adriano Amaricci (adriano.amaricci.AT.sissa.it)
# essentially based on previous work of. Simplified version for SciFortran
# Joan MASSICH (joan.massich-vall.AT.inria.fr).
# Alexandre GRAMFORT (Alexandre.Gramfort.AT.inria.fr)
# Théodore PAPADOPOULO (papadop.AT.inria.fr)



set(CMAKE_FIND_DEBUG_MODE 1)


# Find MKL ROOT
find_path(MKL_ROOT_DIR NAMES include/mkl.h include/mkl.fi  PATHS $ENV{MKLROOT})

# Convert symlinks to real paths

get_filename_component(MKL_ROOT_DIR ${MKL_ROOT_DIR} REALPATH)

if (NOT MKL_ROOT_DIR)
  if (MKL_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find MKL: please set environment variable {MKLROOT}")
  else()
    unset(MKL_ROOT_DIR CACHE)
  endif()

else()
  set(MKL_INCLUDE_DIR ${MKL_ROOT_DIR}/include)
  
  # set arguments to call the MKL provided tool for linking
  set(MKL_LINK_TOOL ${MKL_ROOT_DIR}/tools/mkl_link_tool)
  if (NOT EXISTS "${MKL_LINK_TOOL}")
    message(FATAL_ERROR "cannot find MKL tool: ${MKL_LINK_TOOL}")
  endif()
  
  
  set(MKL_LINK_TOOL_LIBS ${MKL_LINK_TOOL} -libs -l static)
  set(MKL_LINK_TOOL_INCS ${MKL_LINK_TOOL} -opts)
  
  
  execute_process(COMMAND  ${MKL_LINK_TOOL_LIBS}
    OUTPUT_VARIABLE MKL_LIBRARIES
    RESULT_VARIABLE COMMAND_WORKED
    TIMEOUT 2 ERROR_QUIET)
  if (NOT ${COMMAND_WORKED} EQUAL 0)
    message(FATAL_ERROR "Cannot find MKL libraries. The mkl_link_tool command executed was:\n ${MKL_LINK_TOOL_LIBS}.")
  endif()
  
  execute_process(COMMAND ${MKL_LINK_TOOL_INCS}
    OUTPUT_VARIABLE MKL_INCLUDE
    RESULT_VARIABLE COMMAND_WORKED
    TIMEOUT 2 ERROR_QUIET)
  if (NOT ${COMMAND_WORKED} EQUAL 0)
    message(FATAL_ERROR "Cannot find MKL libraries. The mkl_link_tool command executed was:\n ${MKL_LINK_TOOL_INCS}.")
  endif()
  
  
  
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE )
  
  mark_as_advanced(MKL_INCLUDE MKL_LIBRARIES MKL_ROOT_DIR)

  message(STATUS "MKL found at: ${MKL_ROOT_DIR}")

  if (CMAKE_FIND_DEBUG_MODE)
    message(STATUS "Exectuted command: ${MKL_LINK_TOOL_LIBS}; ${MKL_LINK_TOOL_INCS}")
    message(STATUS "Found MKL_LIBRARIES:${MKL_LIBRARIES}")
    message(STATUS "Found MKL_INCLUDE:${MKL_INCLUDE}")
  endif()
  
endif()
