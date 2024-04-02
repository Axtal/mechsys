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

SET(VTK_INCLUDE_PATH
  /usr/include/vtk-5.8
  /usr/local/include/vtk-5.8)

SET(VTK_LIBRARY_PATH
  /usr/lib
  /usr/local/lib/vtk-5.8)

SET(VTK_INCS
    vtkRenderWindow.h)

SET(VTK_LIBS
    MapReduceMPI
    mpistubs
    vtkalglib
    vtkCharts
    vtkCommon
    vtkDICOMParser
    vtkexoIIc
    vtkexpat
    vtkFiltering
    vtkfreetype
    vtkftgl
    vtkGenericFiltering
    vtkGeovis
    vtkGraphics
    vtkhdf5
    vtkHybrid
    vtkImaging
    vtkInfovis
    vtkIO
    vtkjpeg
    vtklibxml2
    vtkmetaio
    vtkNetCDF
    vtkNetCDF_cxx
    vtkpng
    vtkproj4
    vtkRendering
    vtksqlite
    vtksys
    vtktiff
    vtkverdict
    vtkViews
    vtkVolumeRendering
    vtkWidgets
    vtkzlib
    GL
    GLU)

SET(VTK_FOUND 1)

FOREACH(inc ${VTK_INCS})
    FIND_PATH(VTK_${inc} ${inc} ${VTK_INCLUDE_PATH})
    IF(NOT VTK_${inc})
        SET(VTK_FOUND 0)
    ELSE(NOT VTK_${inc})
        SET(VTK_INCLUDE_DIRS ${VTK_INCLUDE_DIRS} ${VTK_${inc}})
    ENDIF(NOT VTK_${inc})
ENDFOREACH(inc)

FOREACH(lib ${VTK_LIBS})
    FIND_LIBRARY(VTK_${lib} NAMES ${lib} PATHS ${VTK_LIBRARY_PATH})
    IF(NOT VTK_${lib})
        SET(VTK_FOUND 0)
    ELSE(NOT VTK_${lib})
        SET(VTK_LIBRARIES ${VTK_LIBRARIES} ${VTK_${lib}})
    ENDIF(NOT VTK_${lib})
ENDFOREACH(lib)

FOREACH(inc ${VTK_INCS})
    UNSET(VTK_${inc})
    UNSET(VTK_${inc} CACHE)
ENDFOREACH(inc)

FOREACH(lib ${VTK_LIBS})
    UNSET(VTK_${lib})
    UNSET(VTK_${lib} CACHE)
ENDFOREACH(lib)

#MESSAGE(${VTK_INCLUDE_DIRS})
#MESSAGE(${VTK_LIBRARIES})
