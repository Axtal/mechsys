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

SET(SOPLEX_INCLUDE_SEARCH_PATH
  $ENV{HOME}/pkg/soplex-1.5.0/src
  /usr/include
  /usr/local/include)

SET(SOPLEX_LIBRARY_SEARCH_PATH
  $ENV{HOME}/pkg/soplex-1.5.0/lib
  /usr/lib
  /usr/local/lib)

FIND_PATH(SOPLEX_SOPLEX_H soplex.h ${SOPLEX_INCLUDE_SEARCH_PATH})
FIND_LIBRARY(SOPLEX_SOPLEX NAMES soplex PATHS ${SOPLEX_LIBRARY_SEARCH_PATH})

SET(SOPLEX_FOUND 1)
FOREACH(var SOPLEX_SOPLEX_H SOPLEX_SOPLEX)
  IF(NOT ${var})
	SET(SOPLEX_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(SOPLEX_FOUND)
  SET(SOPLEX_INCLUDE_DIRS ${SOPLEX_SOPLEX_H})
  SET(SOPLEX_LIBRARIES    ${SOPLEX_SOPLEX})
ENDIF(SOPLEX_FOUND)
