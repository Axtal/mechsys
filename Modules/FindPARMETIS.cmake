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

SET(METIS_INCLUDE_SEARCH_PATH
    /usr/include/metis
    /usr/include/
    $ENV{MECHSYS_ROOT}/pkg/metis-5.1.0/include
    )

SET(METIS_LIBRARY_SEARCH_PATH
    /usr/lib
    $ENV{MECHSYS_ROOT}/pkg/metis-5.1.0/libmetis
    )

FIND_PATH(METIS_H metis.h ${METIS_INCLUDE_SEARCH_PATH})
FIND_LIBRARY(METIS NAMES metis PATHS ${METIS_LIBRARY_SEARCH_PATH})

SET(METIS_FOUND 1)
FOREACH(var METIS)
  IF(NOT ${var})
	SET(METIS_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(METIS_FOUND)
  SET(METIS_INCLUDE_DIRS ${METIS_H})
  SET(METIS_LIBRARIES    ${METIS})
ENDIF(METIS_FOUND)
