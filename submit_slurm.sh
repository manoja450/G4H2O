#!/bin/bash
#SBATCH -J G4h2o
#SBATCH -p shortjobs
#SBATCH -t 8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH -o /home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/NEXTmodify/G4H2O/SLURMOUT/slurm-%j.out
#SBATCH -e /home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/NEXTmodify/G4H2O/SLURMOUT/slurm-%j.out

TIMESTAMP=$(date +%Y%m%d_%H%M%S)

BASE_DIR="/home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/NEXTmodify/G4H2O"
OUTPUT_DIR="${BASE_DIR}/OUTPUT/run_${TIMESTAMP}"
SLURM_DIR="${BASE_DIR}/SLURMOUT"

mkdir -p ${OUTPUT_DIR}
mkdir -p ${SLURM_DIR}

cd ${BASE_DIR}

export LD_LIBRARY_PATH=""
export ROOTSYS="/usr"
export LD_LIBRARY_PATH="/usr/lib64:${LD_LIBRARY_PATH}"

export G4INSTALL="/home/manoja450/geant4-install"
source ${G4INSTALL}/bin/geant4.sh
export LD_LIBRARY_PATH="${G4INSTALL}/lib64:${LD_LIBRARY_PATH}"

export CRYHOME="/home/manoja450/cry_v1.7"
export LD_LIBRARY_PATH="${CRYHOME}/lib:${LD_LIBRARY_PATH}"

export G4NEUTRONHPDATA="/home/manoja450/geant4-data/G4NDL4.7"
export G4LEDATA="/home/manoja450/geant4-data/G4EMLOW8.6.1"
export G4LEVELGAMMADATA="/home/manoja450/geant4-data/PhotonEvaporation5.7"
export G4RADIOACTIVEDATA="/home/manoja450/geant4-data/RadioactiveDecay5.6"
export G4PARTICLEXSDATA="/home/manoja450/geant4-data/G4PARTICLEXS4.1"
export G4PIIDATA="/home/manoja450/geant4-data/G4PII1.3"
export G4REALSURFACEDATA="/home/manoja450/geant4-data/RealSurface2.2"
export G4SAIDXSDATA="/home/manoja450/geant4-data/G4SAIDDATA2.0"
export G4ABLADATA="/home/manoja450/geant4-data/G4ABLA3.3"
export G4INCLDATA="/home/manoja450/geant4-data/G4INCL1.2"

echo "========================================="
echo "Simulation started at $(date)"
echo "Output directory: ${OUTPUT_DIR}"
echo "========================================="

echo "/run/beamOn" > ${OUTPUT_DIR}/input.txt

./build/G4d2o  

echo "========================================="
echo "Simulation finished at $(date)"
echo "========================================="

grep "Number of events processed" ${OUTPUT_DIR}/simulation.log
