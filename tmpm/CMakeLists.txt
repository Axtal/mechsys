########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             #
# Copyright (C) 2009 Sergio Galindo                                    #
# Copyright (C) 2013 William Oquendo                                   #
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
CMAKE_MINIMUM_REQUIRED (VERSION 3.28)

# --- Force the correct CUDA compiler and architecture ---
set(CMAKE_CUDA_COMPILER "/usr/local/cuda-13.0/bin/nvcc" CACHE FILEPATH "CUDA compiler" FORCE)
#set(CMAKE_CUDA_FLAGS "-arch=sm_86" CACHE STRING "CUDA flags" FORCE)
#set(CMAKE_CUDA_ARCHITECTURES "86" CACHE STRING "CUDA architectures" FORCE)
#if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
#    set(CMAKE_CUDA_ARCHITECTURES native)
#endif()


# --- Auto-detect CUDA architecture with fallback ---
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    # Try to detect using nvidia-smi (most common)
    execute_process(
        COMMAND nvidia-smi --query-gpu=compute_cap --format=csv,noheader
        OUTPUT_VARIABLE GPU_CC
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    
    if(GPU_CC AND NOT GPU_CC STREQUAL "N/A")
        string(REPLACE "." "" GPU_CC_NUM ${GPU_CC})
        message(STATUS "Detected GPU: sm_${GPU_CC_NUM}")
        set(CMAKE_CUDA_ARCHITECTURES "${GPU_CC_NUM}" CACHE STRING "CUDA architectures" FORCE)
    else()
        message(STATUS "Could not detect GPU. Using default: sm_61")
        set(CMAKE_CUDA_ARCHITECTURES "61" CACHE STRING "CUDA architectures" FORCE)
    endif()
endif()

PROJECT (MECHSYSSANDBOX)

INCLUDE(${MECHSYSSANDBOX_SOURCE_DIR}/UseMECHSYS.cmake)


SET(PROGS
   tmpm01
   tmpm02
   tmpm03
   tmpm04
)

SET(CUDAPROGS
    tmpm_cu_01
    tmpm_cu_02
    tmpm_cu_03
    tmpm_cu_04
)

FOREACH(var ${PROGS})
    ADD_EXECUTABLE        (${var} "${var}.cpp")
    TARGET_LINK_LIBRARIES (${var} ${LIBS})
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
ENDFOREACH(var)

FOREACH(var ${CUDAPROGS})
    ADD_EXECUTABLE        (${var} "${var}.cu")
    TARGET_LINK_LIBRARIES (${var} ${CUDA_LIBRARIES})
    # Remove any explicit -arch flag; rely on global CMAKE_CUDA_FLAGS
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS} -DUSE_CUDA" LINK_FLAGS "${LFLAGS}")
    # No need for target CUDA_ARCHITECTURES; global setting handles it
ENDFOREACH(var)
