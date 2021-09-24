#!/bin/bash [could also be /bin/tcsh]
#$ -S /bin/bash
#$ -N Winterdilemma_Sweeps
#$ -pe mvapich2-sam 1
#$ -cwd
#$ -t 1:100:1

# avoid multithreading in numpy
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1



python -u ./sweep.py -i $SGE_TASK_ID
