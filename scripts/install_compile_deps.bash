#!/bin/bash

set -e

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

echo
echo "*************************************************************************"
echo "*                                                                       *"
echo "* You can call this script with options to force recompiling all        *"
echo "* dependencies and/or download their code again:                        *"
echo "*                                                                       *"
echo "*                                             recompile     download    *"
echo "* Example:                                            |     |           *"
echo "*                                                     V     V           *"
echo "*   bash ~/mechsys/scripts/install_compile_deps.sh {0,1} {0,1} {NPROCS} *"
echo "*                                                                       *"
echo "* By default, a dependency will not be re-compiled if the corresponding *"
echo "* directory exists under the MECHSYS_ROOT/pkg directory.                *"
echo "*                                                                       *"
echo "* By default, -j4 (4 processors) is given to make. This can be changed  *"
echo "* with a new value of NPROCS.                                           *"
echo "*                                                                       *"
echo "*************************************************************************"

if [ -d "$MECHSYS_ROOT/mechsys" ]; then
    echo
    echo "Found: $MECHSYS_ROOT/mechsys ==> OK"
else
    echo
    echo "Directory named 'mechsys' does not exist"
    echo "Please, download 'mechsys' first:"
    echo
    echo "   hg clone http://hg.savannah.nongnu.org/hgweb/mechsys/"
    echo
    exit 1
fi

RECOMPILE=0
if [ "$#" -gt 0 ]; then
    RECOMPILE=$1
    if [ "$RECOMPILE" -lt 0 -o "$RECOMPILE" -gt 1 ]; then
        echo
        echo "The option for re-compilation must be either 0 or 1"
        echo "  $1 is invalid"
        echo
        exit 1
    fi
fi

FORCEDOWNLOAD=0
if [ "$#" -gt 1 ]; then
    FORCEDOWNLOAD=$2
    if [ "$FORCEDOWNLOAD" -lt 0 -o "$FORCEDOWNLOAD" -gt 1 ]; then
        echo
        echo "The option for downloading and compilation of additional packages must be either 0 or 1"
        echo "  $2 is invalid"
        echo
        exit 1
    fi
fi

NPROCS=8
if [ "$#" -gt 2 ]; then
    NPROCS=$3
fi

test -d $MECHSYS_ROOT/pkg || mkdir $MECHSYS_ROOT/pkg

VER_BLITZ=0.9
VER_TRIANGLE=1.6
VER_TETGEN=1.4.3
#VER_VORO=0.3.1
VER_VORO=0.4.5
VER_IGRAPH=0.8.2
#VER_SOPLEX=1.6.0
#VER_WXWIDGETS=2.9.4
VER_GSL=2.1
VER_BOOST=1_59_0
VER_LAPACK=3.9.0
VER_HDF5=1.8.15-patch1
VER_SZIP=2.1.1

error_message() {
    echo
    echo
    echo "    [1;31m Error: $1 [0m"
    echo
    echo
}

download_and_compile() {
    PKG=""
    PKG_DIR=""
    EXT=tar.gz
    LOCATION=""
    RENAME_PKG=""
    EXTRA_CONF=""
    EXTRA_CMD=""
    CONF_PRMS=""
    IS_SVN=0
    DO_PATCH=0
    DO_CONF=0
    DO_CMAKECONF=0
    DO_MAKE=1
    DO_MAKE_INST=0
    case "$1" in
        blitz)
            PKG=blitz-$VER_BLITZ
            EXT=tar.gz
            LOCATION=http://mechsys.nongnu.org/software/$PKG.$EXT 
            #LOCATION=http://downloads.sourceforge.net/project/blitz/blitz/Blitz%2B%2B%20$VER_BLITZ/blitz-$VER_BLITZ.$EXT 
            DO_CONF=1
            DO_PATCH=1
            ;;
        triangle)
            PKG=triangle$VER_TRIANGLE
            LOCATION=http://mechsys.nongnu.org/software/$PKG.$EXT
            DO_PATCH=1
            ;;
        tetgen)
            PKG=tetgen$VER_TETGEN
            LOCATION=http://mechsys.nongnu.org/software/$PKG.$EXT
            DO_PATCH=1
            ;;
        voro)
            PKG=voro++-$VER_VORO
            LOCATION=http://mechsys.nongnu.org/software/voro++-$VER_VORO.$EXT
            #LOCATION=http://math.lbl.gov/voro++/download/dir/voro++-$VER_VORO.$EXT
            #DO_PATCH=1
            DO_MAKE=1
            ;;
        igraph)
            PKG=igraph-$VER_IGRAPH
            LOCATION=http://mechsys.nongnu.org/software/$PKG.$EXT
            #LOCATION=http://sourceforge.net/projects/igraph/files/C%20library/$VER_IGRAPH/$PKG.$EXT
            CONF_PRMS="--enable-tls"
            DO_CONF=1
            ;;
        #soplex)
            #PKG=soplex-$VER_SOPLEX
            #EXT=tgz
            #LOCATION=http://soplex.zib.de/download/$PKG.$EXT
            #NPROCS=1 # it seems soples cannot be compiled in parallel...
            #;;
        #wxwidgets)
            #PKG=wxWidgets-$VER_WXWIDGETS
            #EXT=tar.bz2
            #LOCATION=http://downloads.sourceforge.net/project/wxwindows/$VER_WXWIDGETS/wxWidgets-$VER_WXWIDGETS.$EXT
            #DO_CONF=1
            #DO_MAKE_INST=1
            #CONF_PRMS="-enable-monolithic --disable-compat26 --disable-compat28 --with-opengl"
            #;;
        gsl)
            PKG=gsl-$VER_GSL
            EXT=tar.gz
            LOCATION=ftp://ftp.gnu.org/gnu/gsl/$PKG.$EXT
            DO_CONF=1
            DO_MAKE=1
            ;;
        boost)
            PKG=boost_$VER_BOOST
            EXT=tar.bz2
            LOCATION=http://internode.dl.sourceforge.net/project/boost/boost/1.59.0/$PKG.$EXT
            DO_CONF=1
            DO_MAKE=1
            ;;
        lapack)
            PKG=lapack-$VER_LAPACK
            EXT=tgz
            LOCATION=http://www.netlib.org/lapack/$PKG.$EXT
            DO_CONF=0
            DO_CMAKECONF=1
            DO_MAKE=1
            ;;
        hdf5)
            PKG=hdf5-$VER_HDF5
            EXT=tar.gz
            LOCATION=http://mechsys.nongnu.org/software/$PKG.$EXT
            #LOCATION=https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-$VER_HDF5/src/$PKG.$EXT
            DO_CONF=1
            DO_MAKE=1
            ;;
        szip)
            PKG=szip-$VER_SZIP
            EXT=tar.gz
            LOCATION=http://mechsys.nongnu.org/software/$PKG.$EXT
            #LOCATION=https://support.hdfgroup.org/ftp/lib-external/szip/$VER_SZIP/src/szip-$VER_SZIP.tar.gz
            DO_CONF=1
            DO_MAKE=1
            ;;
        *)
            error_message "download_and_compile: __Internal_error__"
            exit 1
            ;;
    esac
    echo
    echo "********************************** ${1} ********************************"

    # change into the packages directory
    cd $MECHSYS_ROOT/pkg

    # package filename and directory
    PKG_FILENAME=$PKG.$EXT
    if [ -z "$PKG_DIR" ]; then PKG_DIR=$PKG; fi

    # check for package that must be existing (cannot be downloaded)
    if [ -z "$LOCATION" ]; then
        if [ ! -e "$PKG_FILENAME" ]; then
            error_message "Please download <$PKG_FILENAME> first"
            return
        fi
    fi

    # (re)compile or return (erasing existing package) ?
    if [ "$IS_SVN" -eq 0 ]; then
        if [ -d "$MECHSYS_ROOT/pkg/$PKG_DIR" ]; then
            if [ "$RECOMPILE" -eq 1   -o   "$FORCEDOWNLOAD" -eq 1 ]; then
                echo "    Erasing existing $PKG_DIR"
                rm -rf $MECHSYS_ROOT/pkg/$PKG_DIR
            else
                echo "    Using existing $PKG_DIR"
                return
            fi
        fi
    else
        if [ -d "$MECHSYS_ROOT/pkg/$PKG_DIR" ]; then
            if [ "$RECOMPILE" -eq 1   -o   "$FORCEDOWNLOAD" -eq 1 ]; then
                echo "    Updating existing $PKG SVN repository"
                cd $PKG_DIR
                svn up
                cd $MECHSYS_ROOT/pkg
            else
                echo "    Using existing $PKG SVN repository in $PKG_DIR"
                return
            fi
        fi
    fi

    # download package
    if [ "$IS_SVN" -eq 0 ]; then
        if [ "$FORCEDOWNLOAD" -eq 1   -o   ! -e "$PKG_FILENAME" ]; then
            if [ -e "$PKG_FILENAME" ]; then
                echo "    Removing existing <$PKG_FILENAME>"
                rm $PKG_FILENAME
            fi
            echo "    Downloading <$PKG_FILENAME>"
            wget $LOCATION
            # rename package
            if [ ! -z "$RENAME_PKG" ]; then
                echo "        . . . renaming package . . . . . . ."
                $RENAME_PKG
            fi
        fi
    else
        if [ ! -d "$MECHSYS_ROOT/pkg/$PKG_DIR" ]; then
            echo "    Downloading new $PKG SVN repository"
            svn co $LOCATION $PKG
        fi
    fi

    # uncompress package
    if [ "$IS_SVN" -eq 0 ]; then
        echo "        . . . uncompressing . . ."
        if [ "$EXT" = "tar.bz2" ]; then
            tar xjf $PKG_FILENAME
        else
            tar xzf $PKG_FILENAME
        fi
    fi

    # change into the package directory
    cd $PKG_DIR

    # patch
    if [ "$DO_PATCH" -eq 1 ]; then
        echo "        . . . patching . . ."
        bash $MECHSYS_ROOT/mechsys/patches/${1}/do_patch.bash
    fi

    # configure
    if [ "$DO_CONF" -eq 1 ]; then
        echo "        . . . configuring . . ."
        ./configure $CONF_PRMS 2> /dev/null
    fi

    # extra configuration
    if [ ! -z "$EXTRA_CONF" ]; then
        echo "        . . . extra configuration . . . . . ."
        $EXTRA_CONF
    fi

    # cmake configuration
    if [ "$DO_CMAKECONF" -eq 1 ]; then
        echo "        . . . configuring using cmake . . ."
        cd ..
        mkdir $PKG-build
        cd $PKG-build
        cmake ../$PKG_DIR 2> /dev/null
        #cd ../$PKG_DIR
    fi

    # compilation
    if [ "$DO_MAKE" -eq 1 ]; then
        echo "        . . . compiling . . ."
        make -j$NPROCS > /dev/null 2> /dev/null
    fi

    # make install
    if [ "$DO_MAKE_INST" -eq 1 ]; then
        echo "        . . . installing . . ."
        sudo make install > /dev/null 2> /dev/null
        echo "           . . ldconfig . ."
        sudo ldconfig
    fi

    # execute specific command
    if [ ! -z "$EXTRA_CMD" ]; then
        echo "        . . . command . . . . . ."
        $EXTRA_CMD
    fi

    # finished
    echo "        . . . finished . . . . . "
}

download_and_compile blitz
download_and_compile triangle
download_and_compile tetgen
download_and_compile voro
download_and_compile igraph
#download_and_compile soplex
#download_and_compile wxwidgets
download_and_compile gsl
#download_and_compile boost
download_and_compile lapack
download_and_compile hdf5
download_and_compile szip

echo
echo "Finished ###################################################################"
echo
