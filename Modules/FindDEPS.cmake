#####################################################################################
# MechSys - A C++ library to simulate Mechanical Systems                            #
# Copyright (C) 2010 Sergio Galindo                                                 #
#                                                                                   #
# This file is part of MechSys.                                                     #
#                                                                                   #
# MechSys is free software; you can redistribute it and/or modify it under the      #
# terms of the GNU General Public License as published by the Free Software         #
# Foundation; either version 2 of the License, or (at your option) any later        #
# version.                                                                          #
#                                                                                   #
# MechSys is distributed in the hope that it will be useful, but WITHOUT ANY        #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A   #
# PARTICULAR PURPOSE. See the GNU General Public License for more details.          #
#                                                                                   #
# You should have received a copy of the GNU General Public License along with      #
# MechSys; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, #
# Fifth Floor, Boston, MA 02110-1301, USA                                           #
#####################################################################################

# Flags
OPTION(A_MAKE_VERBOSE        "Show additional messages during compilation/linking?"                         ON )
OPTION(A_MAKE_ALL_WARNINGS   "Make with all warnings (-Wall)"                                               OFF)
OPTION(A_MAKE_DEBUG_SYMBOLS  "Make with debug symbols (-g)"                                                 OFF)
OPTION(A_MAKE_PROFILING      "Make with profiling information (-pg)"                                        OFF)
OPTION(A_MAKE_OPTIMIZED      "Make optimized (-O3)"                                                         ON )
OPTION(A_MAKE_WXW_MONO       "Use wxWidgets monolithic libraries ?"                                         OFF)
OPTION(A_MAKE_TERM_WHITEBG   "Select colors for terminal with white background ?"                           OFF)
OPTION(A_MAKE_TERM_NOCOLORS  "Don't use colors when printing to terminal ?"                                 OFF)
OPTION(A_MAKE_STDVECTOR      "Use std::vector instead of own implemenatation ?"                             ON )
OPTION(A_MAKE_CHECK_OVERLAP  "Check for maximun overlapping in DEM simulations"                             ON )
OPTION(A_MAKE_USE_GPU_DOUBLE "Use double precision numbers in GPU computation"                              ON )
OPTION(A_MAKE_IGNORE_SOLID   "Ignore deep solid cells from LBM computations"                                OFF)
OPTION(A_MAKE_USE_IBB        "Use Immersed Bounce Back instead of Gamma Method for DEM-LBM simulations"     OFF)
                                                                                   
# Options                                                                          
OPTION(A_USE_OMP            "Use OpenMP  ?"                                        ON )
OPTION(A_USE_VTK            "Use VTK ?"                                            OFF)
OPTION(A_USE_HDF5           "Use HDF5 ?"                                           ON )

#ADD_DEFINITIONS(-fmessage-length=0) # Each error message will appear on a single line; no line-wrapping will be done.
#ADD_DEFINITIONS(-std=gnu++11)                   # New C++ standard
#ADD_DEFINITIONS(-std=c++17)                      # New C++ standard
#ADD_DEFINITIONS(-fpermissive)                    # New C++ standard
ADD_DEFINITIONS(-Wno-deprecated-declarations)    # Remove depracated warnings
#INCLUDE_DIRECTORIES (${INCLUDE_DIRECTORIES} $ENV{HOME}/pkg/boost_1_59_0)
#INCLUDE_DIRECTORIES (${INCLUDE_DIRECTORIES} $ENV{MECHSYS_ROOT}/pkg/boost_1_59_0)

ENABLE_LANGUAGE(CUDA)

### FLAGS ###############################################################################################

IF(A_MAKE_VERBOSE)
	SET (CMAKE_VERBOSE_MAKEFILE TRUE)
ENDIF(A_MAKE_VERBOSE)

IF(A_MAKE_ALL_WARNINGS)
	ADD_DEFINITIONS (-Wall)
ENDIF(A_MAKE_ALL_WARNINGS)

IF(A_MAKE_DEBUG_SYMBOLS)
	ADD_DEFINITIONS (-g)
ENDIF(A_MAKE_DEBUG_SYMBOLS)

IF(A_MAKE_PROFILING)
	ADD_DEFINITIONS (-pg)
    SET (LFLAGS "${LFLAGS} -pg")
ENDIF(A_MAKE_PROFILING)

IF(A_MAKE_OPTIMIZED)
	ADD_DEFINITIONS (-O3)
	#ADD_DEFINITIONS (-Ofast)
ENDIF(A_MAKE_OPTIMIZED)

IF(A_MAKE_WXW_MONO)
    SET (WXW_COMPONENTS base core)
ELSE(A_MAKE_WXW_MONO)
    SET (WXW_COMPONENTS base core aui)
ENDIF(A_MAKE_WXW_MONO)

IF(A_MAKE_TERM_WHITEBG AND NOT A_MAKE_TERM_NOCOLORS)
    ADD_DEFINITIONS (-DTERM_WHITEBG)
ENDIF(A_MAKE_TERM_WHITEBG AND NOT A_MAKE_TERM_NOCOLORS)

IF(A_MAKE_TERM_NOCOLORS)
    ADD_DEFINITIONS (-DTERM_NOCOLORS)
ENDIF(A_MAKE_TERM_NOCOLORS)

IF(A_MAKE_STDVECTOR)
    ADD_DEFINITIONS (-DUSE_STDVECTOR)
ENDIF(A_MAKE_STDVECTOR)

IF(A_MAKE_CHECK_OVERLAP)
    ADD_DEFINITIONS (-DUSE_CHECK_OVERLAP)
ENDIF(A_MAKE_CHECK_OVERLAP)

IF(A_MAKE_USE_GPU_DOUBLE)
    ADD_DEFINITIONS (-DUSE_GPU_DOUBLE)
ENDIF(A_MAKE_USE_GPU_DOUBLE)

IF(A_MAKE_IGNORE_SOLID)
    ADD_DEFINITIONS (-DIGNORESOLID)
ENDIF(A_MAKE_IGNORE_SOLID)

IF(A_MAKE_USE_IBB)
    ADD_DEFINITIONS (-DUSE_IBB)
ENDIF(A_MAKE_USE_IBB)

### FIND DEPENDENCIES AND SET FLAGS AND LIBRARIES #######################################################

SET (FLAGS   "${FLAGS}")
SET (LIBS    ${LIBS})
SET (CUDA_LIBRARIES ${CUDA_LIBRARIES})
SET (LFLAGS  "${LFLAGS}")
SET (MISSING "")

SET (Boost_USE_STATIC_LIBS ON)
ENABLE_LANGUAGE (Fortran)

if(A_USE_VTK)
INCLUDE      (FindVTK)                                      #  1
endif(A_USE_VTK)
#FIND_PACKAGE (HDF5 COMPONENTS     HL)                       #  2
#FIND_PACKAGE (HDF5 COMPONENTS CXX HL)                       #  2
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindHDF5.cmake     ) #  2
INCLUDE      (FindBoost)                                    #  3
INCLUDE      (FindLAPACK)                                   #  4
#INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindLocLAPACK.cmake) #  4
#INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindBLAZE.cmake    ) #  5
#INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindARMA.cmake     ) #  5
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindBLITZ.cmake    ) #  5
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindGSL.cmake      ) #  6
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindVORO.cmake     ) # 7
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindTETGEN.cmake   ) # 8
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindTRIANGLE.cmake ) # 9
INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindIGRAPH.cmake   ) # 10
INCLUDE (FindOpenMP )                                       # 11

# 1
if(VTK_FOUND AND A_USE_VTK)
    ADD_DEFINITIONS (-DUSE_VTK)
    INCLUDE_DIRECTORIES (${VTK_INCLUDE_DIRS})
    SET (LIBS  ${LIBS} vtkRendering vtkHybrid) 
    SET (CUDA_LIBRARIES  ${CUDA_LIBRARIES} vtkRendering vtkHybrid)
    SET (FLAGS "${FLAGS} -DVTK_EXCLUDE_STRSTREAM_HEADERS")
else(VTK_FOUND AND A_USE_VTK)
    if(A_USE_VTK)
        SET (MISSING "${MISSING} VTK")
    endif(A_USE_VTK)
endif(VTK_FOUND AND A_USE_VTK)

# 2
if(HDF5_FOUND AND A_USE_HDF5)
    ADD_DEFINITIONS (-DH5_NO_DEPRECATED_SYMBOLS -DH5Gcreate_vers=2 -DH5Gopen_vers=2 -DUSE_HDF5)
	INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
	SET (LIBS ${LIBS} ${HDF5_LIBRARIES})
	SET (CUDA_LIBRARIES ${CUDA_LIBRARIES} ${HDF5_LIBRARIES})
else(HDF5_FOUND AND A_USE_HDF5)
    if(A_USE_HDF5)
        SET (MISSING "${MISSING} HDF5")
    endif(A_USE_HDF5)
endif(HDF5_FOUND AND A_USE_HDF5)

# 4
if(LAPACK_FOUND)
    SET (LIBS ${LIBS} ${LAPACK_LIBRARIES})
    SET (CUDA_LIBRARIES ${CUDA_LIBRARIES} ${LAPACK_LIBRARIES})
else(LAPACK_FOUND)
    INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindLocLAPACK.cmake)
    if(LocLAPACK_FOUND)
        SET (LIBS ${LIBS} ${LocLAPACK_LIBRARIES} "gfortran")
        SET (CUDA_LIBRARIES ${CUDA_LIBRARIES} ${LocLAPACK_LIBRARIES})
    else(LocLAPACK_FOUND)
        SET (MISSING "${MISSING} LaPACK")
    endif(LocLAPACK_FOUND)
endif(LAPACK_FOUND)

# 5
#if(BLAZE_FOUND)
	#INCLUDE_DIRECTORIES (${BLAZE_INCLUDE_DIRS})
#else(BLAZE_FOUND)
    #SET (MISSING "${MISSING} blaze")
#endif(BLAZE_FOUND)

# 5
#if(ARMA_FOUND)
	#INCLUDE_DIRECTORIES (${ARMA_INCLUDE_DIRS})
#else(ARMA_FOUND)
    #SET (MISSING "${MISSING} armadillo")
#endif(ARMA_FOUND)

# 5
if(BLITZ_FOUND)
	INCLUDE_DIRECTORIES (${BLITZ_INCLUDE_DIRS})
else(BLITZ_FOUND)
    SET (MISSING "${MISSING} blitz++")
endif(BLITZ_FOUND)

# 6
if(GSL_FOUND)
	INCLUDE_DIRECTORIES (${GSL_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${GSL_LIBRARIES})
	SET (CUDA_LIBRARIES ${CUDA_LIBRARIES} ${GSL_LIBRARIES})
else(GSL_FOUND)
    SET (MISSING "${MISSING} GSL")
endif(GSL_FOUND)

# 7
if(VORO_FOUND)
	INCLUDE_DIRECTORIES (${VORO_INCLUDE_DIRS})
    SET (LIBS ${LIBS} ${VORO_LIBRARIES})
    SET (CUDA_LIBRARIES ${CUDA_LIBRARIES} ${VORO_LIBRARIES})
else(VORO_FOUND)
    SET (MISSING "${MISSING} Voro++")
endif(VORO_FOUND)

# 8
if(TETGEN_FOUND)
	INCLUDE_DIRECTORIES (${TETGEN_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${TETGEN_LIBRARIES})
	SET (CUDA_LIBRARIES ${CUDA_LIBRARIES} ${TETGEN_LIBRARIES})
else(TETGEN_FOUND)
    SET (MISSING "${MISSING} Tetgen")
endif(TETGEN_FOUND)

# 9
if(TRIANGLE_FOUND)
	INCLUDE_DIRECTORIES (${TRIANGLE_INCLUDE_DIRS})
	SET (LIBS ${LIBS} ${TRIANGLE_LIBRARIES})
	SET (CUDA_LIBRARIES ${CUDA_LIBRARIES} ${TRIANGLE_LIBRARIES})
else(TRIANGLE_FOUND)
    SET (MISSING "${MISSING} Triangle")
endif(TRIANGLE_FOUND)

# 10
if(IGRAPH_FOUND)
    INCLUDE_DIRECTORIES (${IGRAPH_INCLUDE_DIRS})
    SET (LIBS ${LIBS} ${IGRAPH_LIBRARIES})
    SET (CUDA_LIBRARIES ${CUDA_LIBRARIES} ${IGRAPH_LIBRARIES})
	ADD_DEFINITIONS(-DHAS_IGRAPH)
else(IGRAPH_FOUND)
    SET (MISSING "${MISSING} IGraph")
endif(IGRAPH_FOUND)

# 11
if(OPENMP_FOUND AND A_USE_OMP)
    ADD_DEFINITIONS (-DUSE_OMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    SET(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=${OpenMP_CXX_FLAGS}") 
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}") 
else(OPENMP_FOUND AND A_USE_OMP)
    if(A_USE_OMP)
        SET (MISSING "${MISSING} OpenMP")
    endif(A_USE_OMP)
endif(OPENMP_FOUND AND A_USE_OMP)

