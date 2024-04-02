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

SET(HDF5_INCLUDE_SEARCH_PATH
    $ENV{MECHSYS_ROOT}/pkg/szip-2.1.1/src
    $ENV{MECHSYS_ROOT}/pkg/hdf5-1.8.15-patch1/src
    $ENV{MECHSYS_ROOT}/pkg/hdf5-1.8.15-patch1/hl/src
    $ENV{HOME}/pkg/szip-2.1.1/src
    $ENV{HOME}/pkg/hdf5-1.8.15-patch1/src
    $ENV{HOME}/pkg/hdf5-1.8.15-patch1/hl/src
#   $ENV{HOME}/pkg/zlib-1.2.5
    )
  #/usr/include)

SET(HDF5_LIBRARY_SEARCH_PATH
    $ENV{MECHSYS_ROOT}/pkg/szip-2.1.1/src/.libs
    $ENV{MECHSYS_ROOT}/pkg/hdf5-1.8.15-patch1/src/.libs
    $ENV{MECHSYS_ROOT}/pkg/hdf5-1.8.15-patch1/hl/src/.libs
    $ENV{HOME}/pkg/szip-2.1.1/src/.libs
    $ENV{HOME}/pkg/hdf5-1.8.15-patch1/src/.libs
    $ENV{HOME}/pkg/hdf5-1.8.15-patch1/hl/src/.libs
#   $ENV{HOME}/pkg/zlib-1.2.5
    )
  #/usr/lib)

#MESSAGE(${HDF5_HL})

FIND_PATH(HDF5_H    hdf5.h    ${HDF5_INCLUDE_SEARCH_PATH})
FIND_PATH(HDF5_HL_H hdf5_hl.h ${HDF5_INCLUDE_SEARCH_PATH})


#MESSAGE(${HDF5_LIBRARY_SEARCH_PATH})

FIND_LIBRARY(SZ      NAMES sz      PATHS ${HDF5_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(HDF5    NAMES hdf5    PATHS ${HDF5_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(HDF5_HL NAMES hdf5_hl PATHS ${HDF5_LIBRARY_SEARCH_PATH})
#FIND_LIBRARY(LZ      NAMES z       PATHS ${HDF5_LIBRARY_SEARCH_PATH})

SET(HDF5_FOUND 1)
#FOREACH(var LZ SZ HDF5_H HDF5_HL_H HDF5 HDF5_HL)
FOREACH(var SZ HDF5_H HDF5_HL_H HDF5 HDF5_HL)
  IF(NOT ${var})
	SET(HDF5_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(HDF5_FOUND)
  SET(HDF5_INCLUDE_DIR  ${HDF5_H} ${HDF5_HL_H})
#  SET(HDF5_LIBRARIES    ${HDF5_HL} ${HDF5} ${LZ} ${SZ})
  SET(HDF5_LIBRARIES    ${HDF5_HL} ${HDF5} ${SZ})
ENDIF(HDF5_FOUND)

#MESSAGE(${HDF5_HL})
