#!/bin/bash

set -e

#sudo apt-get install software-properties-common
#sudo add-apt-repository ppa:george-edison55/cmake-3.x
#sudo apt-get update

sudo apt-get install \
    wget patch \
    g++ gfortran make cmake-curses-gui \
    libgsl0-dev \
    libboost-python-dev \
    libxml2-dev \
    cmake opencl-headers ocl-icd-opencl-dev \
    nvidia-cuda-dev \
    nvidia-cuda-toolkit \
    libblas-dev liblapack-dev 
    #libhdf5-serial-dev
    #python-tk python-numpy python-scipy python-matplotlib \
    #libvtk6-dev \

#sudo apt-get upgrade

# note: libxml2-dev is for igraph
#       libmumps-dev will install libopenmpi-dev

# subversion 
# libblitz0-dev
# libgtk2.0-dev libfltk1.3-dev libhdf5-serial-dev libxml2-dev libvtk5-dev \
# mencoder pkg-config
