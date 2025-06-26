#! /bin/sh
#SBATCH -p cp4
#SBATCH -N 1
#SBATCH -n 1
module load openblas-0.3.20-gcc-8.5.0-ntzahrh 
BLAS_INCLUDE_PATH=/fs2/software/spack/optv018/linux-rhel8-icelake/gcc-8.5.0/openblas-0.3.20-ntzahrh/include
BLAS_LIBRARY_PATH=/fs2/software/spack/optv018/linux-rhel8-icelake/gcc-8.5.0/openblas-0.3.20-ntzahrh/lib
g++  test_day9.cpp -o test -L $BLAS_LIBRARY_PATH -I $BLAS_INCLUDE_PATH -lopenblas 
srun ./test