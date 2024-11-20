
#include "fftw3.h"
#if defined(__FFTW3_MPI) && defined(__MPI)
#include <fftw3-mpi.h>
//#include "fftw3-mpi_mkl.h"
#endif

#if defined(__CUDA) || defined(__UT_USE_CUDA)
#include "cufft.h"
#include "cuda_runtime.h"
#endif

#if defined(__ROCM) || defined(__UT_USE_ROCM)
#include <hipfft/hipfft.h>
#include <hip/hip_runtime.h>
#endif


#include "module_psi/psi.h"



