#include "fft_base.h"
namespace ModulePW
{
template <typename FPTYPE>
FFT_BASE<FPTYPE>::FFT_BASE()
{
}
template <typename FPTYPE>
FFT_BASE<FPTYPE>::~FFT_BASE()
{
}

template FFT_BASE<float>::FFT_BASE();
template FFT_BASE<double>::FFT_BASE();
template FFT_BASE<float>::~FFT_BASE();
template FFT_BASE<double>::~FFT_BASE();
}