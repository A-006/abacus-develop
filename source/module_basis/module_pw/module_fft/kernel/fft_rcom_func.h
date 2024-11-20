#ifndef FFT_ROCM_FUNC_H
#define FFT_ROCM_FUNC_H
#include <hipfft/hipfft.h>
#include <hip/hip_runtime.h>
static const char* _hipfftGetErrorString(hipfftResult_t error)
{
    switch (error)
    {
    case HIPFFT_SUCCESS:
        return "HIPFFT_SUCCESS";
    case HIPFFT_INVALID_PLAN:
        return "HIPFFT_INVALID_PLAN";
    case HIPFFT_ALLOC_FAILED:
        return "HIPFFT_ALLOC_FAILED";
    case HIPFFT_INVALID_TYPE:
        return "HIPFFT_INVALID_TYPE";
    case HIPFFT_INVALID_VALUE:
        return "HIPFFT_INVALID_VALUE";
    case HIPFFT_INTERNAL_ERROR:
        return "HIPFFT_INTERNAL_ERROR";
    case HIPFFT_EXEC_FAILED:
        return "HIPFFT_EXEC_FAILED";
    case HIPFFT_SETUP_FAILED:
        return "HIPFFT_SETUP_FAILED";
    case HIPFFT_INVALID_SIZE:
        return "HIPFFT_INVALID_SIZE";
    case HIPFFT_UNALIGNED_DATA:
        return "HIPFFT_UNALIGNED_DATA";
    case HIPFFT_INCOMPLETE_PARAMETER_LIST:
        return "HIPFFT_INCOMPLETE_PARAMETER_LIST";
    case HIPFFT_INVALID_DEVICE:
        return "HIPFFT_INVALID_DEVICE";
    case HIPFFT_PARSE_ERROR:
        return "HIPFFT_PARSE_ERROR";
    case HIPFFT_NO_WORKSPACE:
        return "HIPFFT_NO_WORKSPACE";
    case HIPFFT_NOT_IMPLEMENTED:
        return "HIPFFT_NOT_IMPLEMENTED";
    case HIPFFT_NOT_SUPPORTED:
        return "HIPFFT_NOT_SUPPORTED";
    }
    return "<unknown>";
}
#define CHECK_CUFFT(func)                                                                                              \
    {                                                                                                                  \
        hipfftResult_t status = (func);                                                                                \
        if (status != HIPFFT_SUCCESS)                                                                                  \
        {                                                                                                              \
            printf("In File %s : HIPFFT API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,              \
                   _hipfftGetErrorString(status), status);                                                             \
        }                                                                                                              \
    }
#endif // FFT_ROCM_FUNC_H