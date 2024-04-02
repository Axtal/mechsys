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

SET(GSL_INCLUDE_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/gsl-2.1
  $ENV{HOME}/pkg/gsl-2.1
  /usr/include
  /usr/local/include)

SET(GSL_LIBRARY_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/gsl-2.1/.libs
  $ENV{HOME}/pkg/gsl-2.1/.libs
  $ENV{MECHSYS_ROOT}/pkg/gsl-2.1/cblas/.libs
  $ENV{HOME}/pkg/gsl-2.1/cblas/.libs
  /usr/lib
  /usr/local/lib)

FIND_PATH(GSL_MONTE_H gsl/gsl_monte.h       ${GSL_INCLUDE_SEARCH_PATH})
FIND_PATH(GSL_INTEG_H gsl/gsl_integration.h ${GSL_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(GSL_GSL NAMES gsl PATHS ${GSL_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(GSL_CBLASGSL NAMES gslcblas PATHS ${GSL_LIBRARY_SEARCH_PATH})


SET(GSL_FOUND 1)
FOREACH(var GSL_MONTE_H GSL_INTEG_H GSL_GSL)
  IF(NOT ${var})
	SET(GSL_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(GSL_FOUND)
  SET(GSL_INCLUDE_DIRS ${GSL_INTEG_H})
  SET(GSL_LIBRARIES    ${GSL_GSL} ${GSL_CBLASGSL})
ENDIF(GSL_FOUND)
