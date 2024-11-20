#ifndef FFT_CUDA_FUNC_H
#define FFT_CUDA_FUNC_H
#include "cufft.h"
#include "cuda_runtime.h"

static const char* _cufftGetErrorString(cufftResult_t error)
{
    switch (error)
    {
    case CUFFT_SUCCESS:
        return "CUFFT_SUCCESS";
    case CUFFT_INVALID_PLAN:
        return "CUFFT_INVALID_PLAN";
    case CUFFT_ALLOC_FAILED:
        return "CUFFT_ALLOC_FAILED";
    case CUFFT_INVALID_TYPE:
        return "CUFFT_INVALID_TYPE";
    case CUFFT_INVALID_VALUE:
        return "CUFFT_INVALID_VALUE";
    case CUFFT_INTERNAL_ERROR:
        return "CUFFT_INTERNAL_ERROR";
    case CUFFT_EXEC_FAILED:
        return "CUFFT_EXEC_FAILED";
    case CUFFT_SETUP_FAILED:
        return "CUFFT_SETUP_FAILED";
    case CUFFT_INVALID_SIZE:
        return "CUFFT_INVALID_SIZE";
    case CUFFT_UNALIGNED_DATA:
        return "CUFFT_UNALIGNED_DATA";
    case CUFFT_INCOMPLETE_PARAMETER_LIST:
        return "CUFFT_INCOMPLETE_PARAMETER_LIST";
    case CUFFT_INVALID_DEVICE:
        return "CUFFT_INVALID_DEVICE";
    case CUFFT_PARSE_ERROR:
        return "CUFFT_PARSE_ERROR";
    case CUFFT_NO_WORKSPACE:
        return "CUFFT_NO_WORKSPACE";
    case CUFFT_NOT_IMPLEMENTED:
        return "CUFFT_NOT_IMPLEMENTED";
    case CUFFT_LICENSE_ERROR:
        return "CUFFT_LICENSE_ERROR";
    case CUFFT_NOT_SUPPORTED:
        return "CUFFT_NOT_SUPPORTED";
    }
    return "<unknown>";
}

#define CHECK_CUFFT(func)                                                                                              \
    {                                                                                                                  \
        cufftResult_t status = (func);                                                                                 \
        if (status != CUFFT_SUCCESS)                                                                                   \
        {                                                                                                              \
            printf("In File %s : CUFFT API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,               \
                   _cufftGetErrorString(status), status);                                                              \
        }                                                                                                              \
    }
#endif // FFT_CUDA_FUNC_H
