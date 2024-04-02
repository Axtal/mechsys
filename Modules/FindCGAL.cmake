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

SET(CGAL_INCLUDE_SEARCH_PATH
  /usr/include
  /usr/local/include
  $ENV{MECHSYS_ROOT}/cgal-releases-CGAL-4.10/include
  )

SET(CGAL_LIBRARY_SEARCH_PATH
  /usr/lib
  /usr/local/lib
  $ENV{MECHSYS_ROOT}/cgal-releases-CGAL-4.10/lib
  )

FIND_PATH(CGAL_CGAL_H CGAL/Exact_predicates_inexact_constructions_kernel.h ${CGAL_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(CGAL_CGAL NAMES CGAL PATHS ${CGAL_LIBRARY_SEARCH_PATH})

SET(CGAL_FOUND 1)
FOREACH(var CGAL_CGAL_H CGAL_CGAL)
  IF(NOT ${var})
	SET(CGAL_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(CGAL_FOUND)
  SET(CGAL_INCLUDE_DIRS ${CGAL_CGAL_H})
  SET(CGAL_LIBRARIES    ${CGAL_CGAL})
ENDIF(CGAL_FOUND)
