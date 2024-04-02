#!/bin/bash

set -e

sudo apt-get install g++ wget cmake patch gfortran mercurial subversion libgsl0-dev libcgal-dev libvtk5-dev \
                     libblitz0-dev libparmetis-dev cmake-curses-gui libhdf5-serial-dev libsuitesparse-dev libboost-python-dev \
                     python-tk python-numpy python-scipy python-matplotlib openmpi-bin libopenmpi-dev
