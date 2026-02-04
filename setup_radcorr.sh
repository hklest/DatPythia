#!/bin/bash

# Function to add to LD_LIBRARY_PATH if the path is not already present
add_to_ld_library_path() {
  if [[ ":$LD_LIBRARY_PATH:" != *":$1:"* ]]; then
    export LD_LIBRARY_PATH=$1:$LD_LIBRARY_PATH
  fi
}

# Function to add to PATH if the path is not already present
add_to_path() {
  if [[ ":$PATH:" != *":$1:"* ]]; then
    export PATH=$1:$PATH
  fi
}

# Set paths for LHAPDF
add_to_ld_library_path /cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/lhapdf-5.9.1/lib
add_to_path /work/clas12/klest/DatPyth/bin
export LHAPATH=/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/lhapdf-5.9.1/share/lhapdf/PDFsets

# Confirm setup
echo "LD_LIBRARY_PATH:"
echo $LD_LIBRARY_PATH
echo "PATH:"
echo $PATH
echo "PYTHONPATH:"
echo $PYTHONPATH
