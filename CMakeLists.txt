########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2021 SA Galindo Torres                                 #
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

PROJECT                (MECHSYS)
CMAKE_MINIMUM_REQUIRED (VERSION 3.2)
ENABLE_TESTING         ()

# Features                                                                         
OPTION(A_WITH_TLBM          "Compile LBM tests"              ON )
OPTION(A_WITH_TEMLBM        "Compile EMLBM tests"            ON )
OPTION(A_WITH_TADLBM        "Compile ADLBM tests"            ON )
OPTION(A_WITH_TFLBM         "Compile FLBM  tests"            ON )
OPTION(A_WITH_TDEM          "Compile DEM tests"              ON )
OPTION(A_WITH_TNN           "Compile NN tests"               ON )
OPTION(A_WITH_TSPH          "Compile SPH tests"              ON )
OPTION(A_WITH_TDFN          "Compile DFN tests"              ON )
OPTION(A_WITH_TMPM          "Compile MPM tests"              ON )
OPTION(A_WITH_TLBMMPM       "Compile LBMMPM tests "          ON )
OPTION(A_WITH_TLBMDEM       "Compile LBMDEM tests "          ON )

# needed when building Debian package
IF(NOT INSTALL_ROOT)
	SET(INSTALL_ROOT /)
ENDIF(NOT INSTALL_ROOT)

INCLUDE_DIRECTORIES (${MECHSYS_SOURCE_DIR})

### DEPENDENCIES ########################################################################################

INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindDEPS.cmake)
#INCLUDE (${MECHSYS_SOURCE_DIR}/Modules/FindDEPSFL.cmake)
if(MISSING)
    MESSAGE("Missing dependencies =${MISSING}")
endif(MISSING)

### PYTHON MODULE #######################################################################################

#if(A_WITH_PYTHON_MODULE)
    #INCLUDE (FindPythonLibs)
    #if(PYTHONLIBS_FOUND)
        #INCLUDE_DIRECTORIES (${PYTHON_INCLUDE_PATH})
        #EXECUTE_PROCESS (COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
        #FIND_LIBRARY (BOOST_PYTHON NAMES boost_python-mt PATHS /usr/lib)
        #IF(NOT BOOST_PYTHON)
            #FIND_LIBRARY (BOOST_PYTHON NAMES boost_python PATHS /usr/lib)
        #ENDIF(NOT BOOST_PYTHON)
        #SET(LIBS ${LIBS} ${PYTHON_LIBRARIES} ${BOOST_PYTHON})
    #else(PYTHONLIBS_FOUND)
        #MESSAGE("To compile the Python module, boost-python must be installed")
    #endif(PYTHONLIBS_FOUND)
#endif(A_WITH_PYTHON_MODULE)
#
### FEATURES ############################################################################################

SET(ALLDIRS ) # sub-directories

IF(A_WITH_TLBM)
	SET(ALLDIRS ${ALLDIRS} tlbm)
ENDIF(A_WITH_TLBM)

IF(A_WITH_TEMLBM)
	SET(ALLDIRS ${ALLDIRS} temlbm temlbm2)
ENDIF(A_WITH_TEMLBM)

IF(A_WITH_TADLBM)
	SET(ALLDIRS ${ALLDIRS} tadlbm)
ENDIF(A_WITH_TADLBM)

IF(A_WITH_TFLBM)
	SET(ALLDIRS ${ALLDIRS} tflbm)
ENDIF(A_WITH_TFLBM)

IF(A_WITH_TDEM)
	SET(ALLDIRS ${ALLDIRS} tdem)
ENDIF(A_WITH_TDEM)

IF(A_WITH_TNN)
	SET(ALLDIRS ${ALLDIRS} tnn)
ENDIF(A_WITH_TNN)

IF(A_WITH_TSPH)
	SET(ALLDIRS ${ALLDIRS} tsph)
ENDIF(A_WITH_TSPH)

IF(A_WITH_TDFN)
	SET(ALLDIRS ${ALLDIRS} tdfn)
ENDIF(A_WITH_TDFN)

IF(A_WITH_TMPM)
	SET(ALLDIRS ${ALLDIRS} tmpm)
ENDIF(A_WITH_TMPM)

IF(A_WITH_TLBMMPM)
	SET(ALLDIRS ${ALLDIRS} tlbmmpm)
ENDIF(A_WITH_TLBMMPM)

IF(A_WITH_TLBMDEM)
	SET(ALLDIRS ${ALLDIRS} tlbmdem)
ENDIF(A_WITH_TLBMDEM)

IF(ALLDIRS)
    SUBDIRS(${ALLDIRS})
ENDIF(ALLDIRS)


