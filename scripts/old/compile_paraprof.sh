#!/bin/bash

set -e

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

echo
echo "****************************************************************************"
echo "* You can call this script with an option to force recompiling everything  *"
echo "*                                                                          *"
echo "* Example:                                                                 *"
echo "*   sh $MECHSYS_ROOT/mechsys/scripts/do_compile_deps.sh {0,1}              *"
echo "*                                                                          *"
echo "* By default, the code will not be compiled if this was done before.       *"
echo "****************************************************************************"

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
        echo "The option for re-compilation must be either 1 or 0. ($1 is invalid)"
        echo
        exit 1
    fi
fi

test -d $MECHSYS_ROOT/pkg || mkdir $MECHSYS_ROOT/pkg

PDT_VER=3.16
TAU_VER=2.19.2

download_and_compile() {
    DO_CONF=1
    DO_MAKE=1
    DO_LINST=1
    CONF_PRMS=""
    case "$1" in
        pdt)
            PKG=pdt.tgz
            PKG_DIR=pdtoolkit-$PDT_VER
            LOCATION=http://tau.uoregon.edu/$PKG
            ;;
        tau)
            PKG=tau.tgz
            PKG_DIR=tau-$TAU_VER
            LOCATION=http://tau.uoregon.edu/$PKG
            CONF_PRMS="-pdt=$HOME/pkg/pdtoolkit-$PDT_VER/ -mpiinc=/usr/include/mpi/ -mpilib=/usr/lib/"
            ;;
        *)
            echo
            echo "download_and_compile: __internal_error__"
            exit 1
            ;;
    esac
    echo
    echo "********************************** ${1} ********************************"

    # enter pkg
    cd $MECHSYS_ROOT/pkg

    # erase existent directory ?
    if [ -d "$MECHSYS_ROOT/pkg/$PKG_DIR" ]; then
        if [ "$RECOMPILE" -eq 1 ]; then
            echo "    Recompiling $PKG"
            rm -rf $MECHSYS_ROOT/pkg/$PKG_DIR
        else
            echo "    Using existing $PKG_DIR"
            return
        fi
    else
        echo "    New compilation of $PKG"
    fi

    # download and uncompressing package
    if [ -e $PKG ]; then
        echo "    Using local <$PKG>"
    else
        echo "    Downloading package"
        if [ -z "$LOCATION" ]; then
            echo "    Please, download <$PKG> first"
            return
        else
            wget $LOCATION
        fi
    fi
    echo "        . . . uncompressing . . ."
    tar xzf $PKG

    # enter pkg dir
    cd $PKG_DIR

    # configure
    if [ "$DO_CONF" -eq 1 ]; then
        echo "        . . . configuring . . ."
        ./configure $CONF_PRMS 2> /dev/null
    fi

    # compilation
    if [ "$DO_MAKE" -eq 1 ]; then
        echo "        . . . compiling . . ."
        make > /dev/null 2> /dev/null
    fi

    # local installation
    if [ "$DO_LINST" -eq 1 ]; then
        echo "        . . . local installation . . ."
        make install > /dev/null 2> /dev/null
    fi
}

download_and_compile pdt
download_and_compile tau

echo
echo "Finished ###################################################################"
echo
