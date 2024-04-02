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

SET(MTL_INCLUDE_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/mtl4
  $ENV{HOME}/pkg/mtl4
  /usr/include
  /usr/local/include)

FIND_PATH(MTL_MTL_H boost/numeric/mtl/mtl.hpp ${MTL_INCLUDE_SEARCH_PATH})

SET(MTL_FOUND 1)
FOREACH(var MTL_MTL_H)
  IF(NOT ${var})
	SET(MTL_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(MTL_FOUND)
  SET(MTL_INCLUDE_DIRS ${MTL_MTL_H})
ENDIF(MTL_FOUND)
