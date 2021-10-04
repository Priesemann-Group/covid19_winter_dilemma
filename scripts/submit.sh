#!/bin/bash [could also be /bin/tcsh]
#$ -S /bin/bash
#$ -N Winterdilemma_Sweeps
#$ -pe mvapich2-sam 32
#$ -cwd
#$ -o ./log/outputs
#$ -e ./log/errors
#$ -t 1:2532:1

# avoid multithreading in numpy
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

# >>>  conda initialize >>>
. $HOME/anaconda3/etc/profile.d/conda.sh
conda activate soccer
# # >>>  conda initialize >>>

 cd $HOME/Repositories/covid19_winter_dilemma/scripts/

python -u ./sweep_r.py -i $SGE_TASK_ID
