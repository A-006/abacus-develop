#!/bin/bash
> abacus.txt
# rm -rf build/
#source /opt/intel/oneapi/setvars.sh intel64 --force
cmake -B build -DCMAKE_CXX_COMPILER=/usr/bin/g++  \
    -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx \
    -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.4/bin/nvcc \
    -DUSE_CUDA=OFF \
    -DENABLE_LIBXC=ON \
   #  -DENABLE_MPI=OFF \
   #  -DENABLE_LCAO=OFF \
   #  -DUSE_OPENMP=OFF
   #  -DENABLE_LIBRI=ON 
    # -DMPI_CXX_LIBRARIES="-L/opt/intel/oneapi/mpi/2021.4.0/lib -lmpi -lpmpi" \
    
cmake --build build -j50
# cmake --install build
# cd /home/ubuntu/desktop/github/abacus/A-006/abacus-develop/tests/integrate/212_NO_wfc_ienvelope
cd ./examples/scf/lcao_Cu
# cd /home/ubuntu/desktop/github/abacus/A-006/abacus-develop/test
# cd ./tests/performance/P101_si32_lcao
# cd /home/ubuntu/desktop/github/abacus/A-006/abacus-develop/tests/integrate/207_NO_OK
OMP_NUM_THREADS=8 mpirun -n 4 /home/ubuntu/desktop/github/abacus/A-006/abacus-develop/build/abacus >> \
   /home/ubuntu/desktop/github/abacus/A-006/abacus-develop/abacus.txt
# OMP_NUM_THREADS=8 mpirun -n 4 /home/ubuntu/desktop/github/abacus/deepmodeling/abacus-develop/build/abacus >> \
#    /home/ubuntu/desktop/github/abacus/A-006/abacus-develop/abacus.txt
# nsys profile --force-overwrite true --trace=cuda -o report /home/ubuntu/desktop/github/cpu_latest/develop/abacus-develop/build/abacus
# OMP_NUM_THREADS=2 mpirun -n 1 gdb /home/ubuntu/desktop/github/abacus/A-006/abacus-develop/build/abacus >> \

# dvlcal  TOTAL-PRESSURE: -3819.245352 KPAR
