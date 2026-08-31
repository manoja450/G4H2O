#!/bin/bash

# --- 1. Environment Setup ---

export PREFIX="/raid1/general"

# Only source these if they exist
if [ -f "$PREFIX/ROOT_install/bin/thisroot.sh" ]; then
    source "$PREFIX/ROOT_install/bin/thisroot.sh"
fi

if [ -f "$PREFIX/Geant4_11_2_1/bin/geant4.sh" ]; then
    source "$PREFIX/Geant4_11_2_1/bin/geant4.sh"
fi

# Your CRY paths
export CRYHOME="/home/manoja450/cry_v1.7"
export CRYDATA="/home/manoja450/cry_v1.7/data"
export LD_LIBRARY_PATH="${CRYHOME}/lib:${LD_LIBRARY_PATH}"

# Project directory = current directory
export G4D2O_DIR="$(pwd)"
export DUMMY_LIB_DIR="${G4D2O_DIR}/dummy_libs"
export LD_LIBRARY_PATH="${DUMMY_LIB_DIR}:${LD_LIBRARY_PATH}"

# Only add OpenSSL path if it exists
if [ -d "/home/manoja450/openssl/lib" ]; then
    export LD_LIBRARY_PATH="/home/manoja450/openssl/lib:${LD_LIBRARY_PATH}"
fi

# --- 2. Check for beamOn.dat ---
if [ ! -f "beamOn.dat" ]; then
    echo "Error: beamOn.dat not found in current directory!"
    echo "Please copy beamOn.dat.bak to beamOn.dat and configure it."
    exit 1
fi

# --- 3. Run the Simulation ---
echo "Running G4d2o locally..."
echo "Output will be saved to output.log"

./G4d2o

echo "Simulation complete. Check output.log."
