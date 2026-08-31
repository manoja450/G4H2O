#!/bin/bash
#SBATCH -J G4dataAnalysis
#SBATCH -p shortjobs
#SBATCH -t 1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH -o /home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/NEXTmodify/G4d2o_DATA_DRIVEN_COPY/mac/SLURMOUT/angular_analysis-%j.out
#SBATCH -e /home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/NEXTmodify/G4d2o_DATA_DRIVEN_COPY/mac/SLURMOUT/angular_analysis-%j.out

# Get timestamp for this run
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Base directory (mac folder)
BASE_DIR="/home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/NEXTmodify/G4d2o_DATA_DRIVEN_COPY/mac"
SLURM_DIR="${BASE_DIR}/SLURMOUT"

# Create SLURM output directory if it doesn't exist
mkdir -p ${SLURM_DIR}

# Change to base directory
cd ${BASE_DIR}

echo "========================================="
echo "Angular Analysis started at $(date)"
echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "Working directory: $(pwd)"
echo "========================================="

# Set up Python environment
# If using conda (uncomment and modify as needed):
# source /home/manoja450/miniconda3/etc/profile.d/conda.sh
# conda activate your_env_name

# If using virtualenv (uncomment and modify as needed):
# source /path/to/your/venv/bin/activate

# Print Python version and packages for debugging
echo "Python version:"
python3 --version

echo -e "\nInstalled packages:"
python3 -c "import numpy, uproot, matplotlib, sys; print(f'Python: {sys.version}'); print(f'Numpy: {numpy.__version__}'); print(f'Uproot: {uproot.__version__}'); print(f'Matplotlib: {matplotlib.__version__}')" 2>/dev/null || echo "Some packages not found"

echo "========================================="
echo "Running angular_analysis.py..."
echo "========================================="

# Run the Python script
#python3 tighttoleranceAngular.py
python3 angular_analysis.py
#python3 angular_analysisV2.py
#python3 two_models_comparison.py
# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "========================================="
    echo "✅ Angular analysis completed successfully at $(date)"
    echo "========================================="
    
    # Find the latest output directory
    LATEST_OUTPUT=$(ls -td ${BASE_DIR}/*_angular_* 2>/dev/null | head -1)
    
    if [ -n "${LATEST_OUTPUT}" ]; then
        echo "Output files saved to: ${LATEST_OUTPUT}"
        echo -e "\nGenerated files:"
        ls -la ${LATEST_OUTPUT}/ 2>/dev/null || echo "No files found"
        echo -e "\nTotal size of output:"
        du -sh ${LATEST_OUTPUT} 2>/dev/null
    else
        echo "Warning: Could not find output directory"
    fi
else
    echo "========================================="
    echo "❌ Angular analysis FAILED at $(date)"
    echo "========================================="
    echo "Check the error log for details"
    exit 1
fi

echo "========================================="
echo "SLURM job completed at $(date)"
echo "========================================="
