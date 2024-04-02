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

SET(UMFPACK_INCLUDE_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/SuiteSparse/UFconfig
  $ENV{MECHSYS_ROOT}/pkg/SuiteSparse/AMD/Include
  $ENV{MECHSYS_ROOT}/pkg/SuiteSparse/UMFPACK/Include
  $ENV{HOME}/pkg/SuiteSparse/UFconfig
  $ENV{HOME}/pkg/SuiteSparse/AMD/Include
  $ENV{HOME}/pkg/SuiteSparse/UMFPACK/Include
  /usr/include/umfpack
  /usr/include/suitesparse)

SET(UMFPACK_LIBRARY_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/SuiteSparse/AMD/Lib
  $ENV{MECHSYS_ROOT}/pkg/SuiteSparse/UMFPACK/Lib
  $ENV{HOME}/pkg/SuiteSparse/AMD/Lib
  $ENV{HOME}/pkg/SuiteSparse/UMFPACK/Lib
  /usr/lib)

FIND_PATH(UMFPACK_AMD_H      amd.h      ${UMFPACK_INCLUDE_SEARCH_PATH})
FIND_PATH(UMFPACK_UMFPACK_H  umfpack.h  ${UMFPACK_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(UMFPACK_UMFPACK NAMES umfpack PATHS ${UMFPACK_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(UMFPACK_AMD     NAMES amd     PATHS ${UMFPACK_LIBRARY_SEARCH_PATH})

SET(UMFPACK_FOUND 1)
FOREACH(var UMFPACK_UMFPACK_H UMFPACK_AMD_H UMFPACK_UMFPACK UMFPACK_AMD)
  IF(NOT ${var})
	SET(UMFPACK_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(UMFPACK_FOUND)
  SET(UMFPACK_INCLUDE_DIRS ${UMFPACK_AMD_H} ${UMFPACK_UMFPACK_H})
  SET(UMFPACK_LIBRARIES    ${UMFPACK_UMFPACK} ${UMFPACK_AMD})
ENDIF(UMFPACK_FOUND)
