########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2021 SA Galindo Torres                                 #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

CMAKE_MINIMUM_REQUIRED (VERSION 3.5)
PROJECT                (MECHSYS)
ENABLE_TESTING         ()

# Features                                                                         
OPTION(A_WITH_TLBM          "Compile LBM tests"              ON )
OPTION(A_WITH_TEMLBM        "Compile EMLBM tests"            ON )
OPTION(A_WITH_TADLBM        "Compile ADLBM tests"            ON )
OPTION(A_WITH_TFLBM         "Compile FLBM  tests"            ON )
OPTION(A_WITH_TDEM          "Compile DEM tests"              ON )
OPTION(A_WITH_TNN           "Compile NN tests"               ON )
OPTION(A_WITH_TSPH          "Compile SPH tests"              ON )
OPTION(A_WITH_TDFN          "Compile DFN tests"              ON )
OPTION(A_WITH_TMPM          "Compile MPM tests"              ON )
OPTION(A_WITH_TLBMMPM       "Compile LBMMPM tests "          ON )
OPTION(A_WITH_TLBMDEM       "Compile LBMDEM tests "          ON )

# needed when building Debian package
IF(NOT INSTALL_ROOT)
	SET(INSTALL_ROOT /)
ENDIF(NOT INSTALL_ROOT)

INCLUDE_DIRECTORIES (${MECHSYS_SOURCE_DIR})

### DEPENDENCIES ########################################################################################

INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindDEPS.cmake)
#INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindDEPSFL.cmake)
if(MISSING)
    MESSAGE("Missing dependencies =${MISSING}")
endif(MISSING)

### PYTHON MODULE #######################################################################################

#if(A_WITH_PYTHON_MODULE)
    #INCLUDE (FindPythonLibs)
    #if(PYTHONLIBS_FOUND)
        #INCLUDE_DIRECTORIES (${PYTHON_INCLUDE_PATH})
        #EXECUTE_PROCESS (COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
        #FIND_LIBRARY (BOOST_PYTHON NAMES boost_python-mt PATHS /usr/lib)
        #IF(NOT BOOST_PYTHON)
            #FIND_LIBRARY (BOOST_PYTHON NAMES boost_python PATHS /usr/lib)
        #ENDIF(NOT BOOST_PYTHON)
        #SET(LIBS ${LIBS} ${PYTHON_LIBRARIES} ${BOOST_PYTHON})
    #else(PYTHONLIBS_FOUND)
        #MESSAGE("To compile the Python module, boost-python must be installed")
    #endif(PYTHONLIBS_FOUND)
#endif(A_WITH_PYTHON_MODULE)
#
### FEATURES ############################################################################################

SET(ALLDIRS ) # sub-directories

IF(A_WITH_TLBM)
	SET(ALLDIRS ${ALLDIRS} tlbm)
ENDIF(A_WITH_TLBM)

IF(A_WITH_TEMLBM)
	SET(ALLDIRS ${ALLDIRS} temlbm temlbm2)
ENDIF(A_WITH_TEMLBM)

IF(A_WITH_TADLBM)
	SET(ALLDIRS ${ALLDIRS} tadlbm)
ENDIF(A_WITH_TADLBM)

IF(A_WITH_TFLBM)
	SET(ALLDIRS ${ALLDIRS} tflbm)
ENDIF(A_WITH_TFLBM)

IF(A_WITH_TDEM)
	SET(ALLDIRS ${ALLDIRS} tdem)
ENDIF(A_WITH_TDEM)

IF(A_WITH_TNN)
	SET(ALLDIRS ${ALLDIRS} tnn)
ENDIF(A_WITH_TNN)

IF(A_WITH_TSPH)
	SET(ALLDIRS ${ALLDIRS} tsph)
ENDIF(A_WITH_TSPH)

IF(A_WITH_TDFN)
	SET(ALLDIRS ${ALLDIRS} tdfn)
ENDIF(A_WITH_TDFN)

IF(A_WITH_TMPM)
	SET(ALLDIRS ${ALLDIRS} tmpm)
ENDIF(A_WITH_TMPM)

IF(A_WITH_TLBMMPM)
	SET(ALLDIRS ${ALLDIRS} tlbmmpm)
ENDIF(A_WITH_TLBMMPM)

IF(A_WITH_TLBMDEM)
	SET(ALLDIRS ${ALLDIRS} tlbmdem)
ENDIF(A_WITH_TLBMDEM)

### CUDA ARCHITECTURES ##################################################################################

# MECHSYS_CUDA_ARCH is the architecture list every CUDA target in this tree
# builds for.  Leave it empty (the default) and the A_MAKE_CUDA_NATIVE switch
# decides; set it to name architectures explicitly.  It is a normal CMake
# architecture list: numbers, "N-real", "N-virtual", or the keywords
# "native"/"all"/"all-major".
#
#   A_MAKE_CUDA_NATIVE  (cmake menu, default ON)
#     ON   build for the GPU(s) in this machine -- tuned, small, fast to build.
#          If the machine has no visible GPU (an HPC compile node) this falls
#          back to the major list below rather than failing.
#     OFF  build for all major architectures the toolkit supports -- portable,
#          runs anywhere from Pascal up, ~4x the CUDA compile time.
#
#   MECHSYS_CUDA_ARCH   (override, default empty)
#     Set this when the switch is not enough, which on a cluster it usually is
#     not: the compile node has no GPU, so the run node's architecture has to be
#     *selected* rather than detected.
#
#       cmake -DMECHSYS_CUDA_ARCH=90 ..            Hopper only
#       cmake -DMECHSYS_CUDA_ARCH="80;90" ..       fat binary, Ampere + Hopper
#       cmake -DMECHSYS_CUDA_ARCH=80-real ..       SASS only (no PTX for JIT)
#       cmake -DMECHSYS_CUDA_ARCH=90-virtual ..    PTX only
#
# Any architecture from sm_60 upwards may be named, whether or not a GPU is
# present and whether or not it belongs to this machine.  Anything below sm_60
# is dropped with a warning: the solvers call atomicAdd(double*,double), which
# arrived with Pascal, so those builds cannot work.
SET(MECHSYS_CUDA_ARCH "" CACHE STRING
    "Override the CUDA architecture list, e.g. 80, 90 or \"80;90\". Empty = use A_MAKE_CUDA_NATIVE.")

# Resolution order, first match wins:
#   1. MECHSYS_CUDA_ARCH                named explicitly by the caller
#   2. -DCMAKE_CUDA_ARCHITECTURES=...   the standard CMake setting, if given
#   3. A_MAKE_CUDA_NATIVE=ON            nvidia-smi on the build machine
#   4. otherwise                        every major arch the installed toolkit
#                                       can target
#
# This is resolved at *configure* time deliberately.  nvcc's own "native"
# resolution happens at compile time and, when no GPU is visible then, silently
# falls back to a pre-Maxwell default (sm_52) where atomicAdd(double*) does not
# exist and the build dies with ~100 errors.  Resolving here means a tree
# configured with a card present still builds correctly without one.
#
# Nothing below is tied to a particular card: on this machine step 3 yields 86,
# on an A100 host it yields 80, on an H100 it yields 90, and on a machine with
# no NVIDIA driver step 4 yields whatever the local toolkit supports.

SET(MECHSYS_CUDA_MIN_ARCH 60)

# CMake's "native"/"all"/"all-major" keywords are expanded here instead of being
# handed to nvcc, for two reasons that both bite exactly on a cluster:
#
#   * nvcc expands "all"/"all-major" starting at sm_50, and this code cannot
#     compile below sm_60, so the keyword would fail the build;
#   * nvcc resolves "native" when it runs, so on a compile node with no GPU it
#     warns "Cannot find valid GPU for '-arch=native', default arch is used",
#     falls back to sm_52 and then fails on atomicAdd(double*).  CMake records
#     the same failure as the placeholder string "No CUDA devices found.-real".
FUNCTION(MECHSYS_EXPAND_CUDA_ARCHS _out)
    SET(_expanded "")
    FOREACH(_a ${ARGN})
        IF(_a STREQUAL "native")
            LIST(APPEND _expanded ${CMAKE_CUDA_ARCHITECTURES_NATIVE})
        ELSEIF(_a STREQUAL "all-major")
            LIST(APPEND _expanded ${CMAKE_CUDA_ARCHITECTURES_ALL_MAJOR})
        ELSEIF(_a STREQUAL "all")
            LIST(APPEND _expanded ${CMAKE_CUDA_ARCHITECTURES_ALL})
        ELSE()
            LIST(APPEND _expanded "${_a}")
        ENDIF()
    ENDFOREACH()
    SET(${_out} "${_expanded}" PARENT_SCOPE)
ENDFUNCTION()

# Keep only numeric archs at or above the floor.  Entries may carry a CMake
# suffix ("80-real", "90-virtual"); anything unrecognisable -- including CMake's
# "No CUDA devices found.-real" placeholder -- is dropped.
FUNCTION(MECHSYS_FILTER_CUDA_ARCHS _out)
    SET(_kept "")
    FOREACH(_a ${ARGN})
        IF(_a MATCHES "^([0-9]+)(-real|-virtual)?$")
            IF(CMAKE_MATCH_1 GREATER_EQUAL MECHSYS_CUDA_MIN_ARCH)
                LIST(APPEND _kept "${_a}")
            ENDIF()
        ENDIF()
    ENDFOREACH()
    IF(_kept)
        LIST(REMOVE_DUPLICATES _kept)
    ENDIF()
    SET(${_out} "${_kept}" PARENT_SCOPE)
ENDFUNCTION()

# (1) the caller named the architectures explicitly -- the HPC path
IF(NOT MECHSYS_CUDA_ARCH STREQUAL "")
    SET(_mechsys_request "${MECHSYS_CUDA_ARCH}")
    MECHSYS_EXPAND_CUDA_ARCHS(_mechsys_exp ${_mechsys_request})
    MECHSYS_FILTER_CUDA_ARCHS(_mechsys_archs ${_mechsys_exp})
    LIST(LENGTH _mechsys_archs _mechsys_nkept)
    IF(_mechsys_nkept EQUAL 0)
        MESSAGE(WARNING
            "MECHSYS_CUDA_ARCH=\"${_mechsys_request}\" names no usable architecture. "
            "MechSys needs sm_60 or newer (atomicAdd(double*) arrived with Pascal). "
            "If you asked for \"native\", note it cannot be resolved on a host with no "
            "visible GPU -- name the run node's architecture instead, e.g. "
            "-DMECHSYS_CUDA_ARCH=90.  Falling back to auto-detection.")
        # Shadow the cache entry with an empty normal variable, otherwise the
        # auto-detection steps below would still see the rejected value.
        SET(MECHSYS_CUDA_ARCH "")
    ELSE()
        LIST(LENGTH _mechsys_exp _mechsys_nexp)
        IF(NOT _mechsys_nkept EQUAL _mechsys_nexp)
            MESSAGE(WARNING
                "MECHSYS_CUDA_ARCH=\"${_mechsys_request}\": architectures below sm_60 were "
                "dropped; building for \"${_mechsys_archs}\" instead.")
        ENDIF()
        SET(MECHSYS_CUDA_ARCH "${_mechsys_archs}")
    ENDIF()
ENDIF()

# (2) the caller used the standard CMake knob
IF(MECHSYS_CUDA_ARCH STREQUAL "" AND
   CMAKE_CUDA_ARCHITECTURES AND NOT CMAKE_CUDA_ARCHITECTURES STREQUAL "")
    MECHSYS_EXPAND_CUDA_ARCHS(_mechsys_exp ${CMAKE_CUDA_ARCHITECTURES})
    MECHSYS_FILTER_CUDA_ARCHS(_mechsys_archs ${_mechsys_exp})
    SET(MECHSYS_CUDA_ARCH "${_mechsys_archs}")
ENDIF()

# (3) A_MAKE_CUDA_NATIVE is ON (the default): this machine's GPU(s).
#     One line per card, "8.6" -> "86".
IF(MECHSYS_CUDA_ARCH STREQUAL "" AND A_MAKE_CUDA_NATIVE)
    FIND_PROGRAM(MECHSYS_NVIDIA_SMI nvidia-smi)
    IF(MECHSYS_NVIDIA_SMI)
        EXECUTE_PROCESS(COMMAND ${MECHSYS_NVIDIA_SMI} --query-gpu=compute_cap --format=csv,noheader
                        OUTPUT_VARIABLE _mechsys_cc
                        RESULT_VARIABLE _mechsys_rc
                        ERROR_QUIET)
        IF(_mechsys_rc EQUAL 0 AND _mechsys_cc)
            STRING(REGEX REPLACE "[\r\n]+" ";" _mechsys_cc "${_mechsys_cc}")
            SET(_mechsys_archs "")
            FOREACH(_cc ${_mechsys_cc})
                STRING(REGEX REPLACE "[^0-9]" "" _cc "${_cc}")
                IF(_cc MATCHES "^[0-9]+$")
                    LIST(APPEND _mechsys_archs "${_cc}")
                ENDIF()
            ENDFOREACH()
            # A heterogeneous multi-GPU host keeps one entry per distinct card.
            MECHSYS_FILTER_CUDA_ARCHS(_mechsys_archs ${_mechsys_archs})
            SET(MECHSYS_CUDA_ARCH "${_mechsys_archs}")
        ENDIF()
    ENDIF()
ENDIF()

# (4) A_MAKE_CUDA_NATIVE is OFF, or it is ON but this host has no visible GPU
#     (the usual case on an HPC compile node).  Build for every major
#     architecture the installed toolkit can target, which is portable: sm_60
#     SASS also serves 61/62, sm_70 serves 72/75, sm_80 serves 86/87/89, and the
#     newest PTX JITs forward.
IF(MECHSYS_CUDA_ARCH STREQUAL "")
    MECHSYS_FILTER_CUDA_ARCHS(_mechsys_archs ${CMAKE_CUDA_ARCHITECTURES_ALL_MAJOR})
    IF(NOT _mechsys_archs)
        MECHSYS_FILTER_CUDA_ARCHS(_mechsys_archs ${CMAKE_CUDA_ARCHITECTURES_ALL})
    ENDIF()
    IF(NOT _mechsys_archs)
        SET(_mechsys_archs "${MECHSYS_CUDA_MIN_ARCH}")
    ENDIF()
    SET(MECHSYS_CUDA_ARCH "${_mechsys_archs}")
ENDIF()

# Warn early if an explicitly requested arch is newer than the installed
# toolkit, rather than letting nvcc fail deep into the first kernel.
IF(CMAKE_CUDA_ARCHITECTURES_ALL)
    SET(_mechsys_known "")
    FOREACH(_a ${CMAKE_CUDA_ARCHITECTURES_ALL})
        IF(_a MATCHES "^([0-9]+)")
            LIST(APPEND _mechsys_known "${CMAKE_MATCH_1}")
        ENDIF()
    ENDFOREACH()
    FOREACH(_a ${MECHSYS_CUDA_ARCH})
        IF(_a MATCHES "^([0-9]+)")
            IF(NOT CMAKE_MATCH_1 IN_LIST _mechsys_known)
                MESSAGE(WARNING
                    "MECHSYS_CUDA_ARCH asks for sm_${CMAKE_MATCH_1}, which the installed CUDA "
                    "toolkit does not list (it supports: ${_mechsys_known}). nvcc will probably "
                    "reject it -- check the toolkit on this build host.")
            ENDIF()
        ENDIF()
    ENDFOREACH()
ENDIF()

# The CUDA_ARCHITECTURES target property needs CMake >= 3.18.  Older CMake is
# still common on clusters, so rather than abandoning the kernel (which would
# fall back to sm_52 and fail), translate the list into -gencode flags here and
# skip the property.  Provenance: sm_60 SASS also serves 61/62, sm_70 serves
# 72/75, sm_80 serves 86/87/89, and PTX JITs forward to later architectures.
FUNCTION(MECHSYS_CUDA_GENCODE_FLAGS _out)
    SET(_flags "")
    FOREACH(_a ${ARGN})
        IF(_a MATCHES "^([0-9]+)(-real|-virtual)?$")
            SET(_n "${CMAKE_MATCH_1}")
            IF(CMAKE_MATCH_2 STREQUAL "-real")
                LIST(APPEND _flags "-gencode=arch=compute_${_n},code=sm_${_n}")
            ELSEIF(CMAKE_MATCH_2 STREQUAL "-virtual")
                LIST(APPEND _flags "-gencode=arch=compute_${_n},code=compute_${_n}")
            ELSE()
                LIST(APPEND _flags "-gencode=arch=compute_${_n},code=[compute_${_n},sm_${_n}]")
            ENDIF()
        ENDIF()
    ENDFOREACH()
    SET(${_out} "${_flags}" PARENT_SCOPE)
ENDFUNCTION()

IF(CMAKE_VERSION VERSION_LESS 3.18)
    MECHSYS_CUDA_GENCODE_FLAGS(_mechsys_gencode ${MECHSYS_CUDA_ARCH})
    STRING(REPLACE ";" " " _mechsys_gencode "${_mechsys_gencode}")
    SET(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${_mechsys_gencode}")
    SET(MECHSYS_CUDA_ARCH_VIA_FLAGS TRUE)
    MESSAGE(STATUS "CMake ${CMAKE_VERSION} has no CUDA_ARCHITECTURES property; passing -gencode directly: ${_mechsys_gencode}")
ELSE()
    SET(MECHSYS_CUDA_ARCH_VIA_FLAGS FALSE)
ENDIF()

MESSAGE(STATUS "CUDA architectures: ${MECHSYS_CUDA_ARCH}")

IF(ALLDIRS)
    SUBDIRS(${ALLDIRS})
ENDIF(ALLDIRS)


