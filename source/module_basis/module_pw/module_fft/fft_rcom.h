#include "fft_base.h"
#include <hipfft/hipfft.h>
#include <hip/hip_runtime.h>
#ifndef FFT_ROCM_H
#define FFT_ROCM_H
namespace ModulePW
{
template <typename FPTYPE>
class FFT_ROCM : public FFT_BASE<FPTYPE>
{
    public:
        FFT_ROCM();
        ~FFT_ROCM(); 
        //init fftw_plans
        void setupFFT() override; 

        void clear() override;

        void cleanFFT() override;

        std::complex<FPTYPE>* get_auxr_3d_data() const override;
        
        void fft3D_forward(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;

        void fft3D_backward(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;
    private:
        hipfftHandle c_handle = {};
        hipfftHandle z_handle = {};
        mutable std::complex<float>* c_auxr_3d = nullptr;  // fft space
        mutable std::complex<double>* z_auxr_3d = nullptr; // fft space

};
}// namespace ModulePW
#endif
