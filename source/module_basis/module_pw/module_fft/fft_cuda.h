#include "fft_base.h"
#include "cufft.h"
#include "cuda_runtime.h"

#ifndef FFT_CUDA_H
#define FFT_CUDA_H
namespace ModulePW
{
template <typename FPTYPE>
class FFT_CUDA : public FFT_BASE<FPTYPE>
{
    public:
        FFT_CUDA();
        ~FFT_CUDA(); 
        //init fftw_plans
	    void setupFFT() override; 

        void clear() override;

        void cleanFFT() override;

        std::complex<FPTYPE>* get_auxr_3d_data() const override;
        
        void fft3D_forward(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;

        void fft3D_backward(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;
    private:
        cufftHandle c_handle = {};
        cufftHandle z_handle = {};
       
        std::complex<float>* c_auxr_3d = nullptr;  // fft space
        std::complex<double>* z_auxr_3d = nullptr; // fft space

};
} // namespace ModulePW
#endif