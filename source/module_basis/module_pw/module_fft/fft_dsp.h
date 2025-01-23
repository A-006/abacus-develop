#ifndef FFT_CUDA_H
#define FFT_CUDA_H

#include "fft_base.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "hthread_host.h"
#include "mtfft.h"
#include "fftw3.h"

namespace ModulePW
{
template <typename FPTYPE>
class FFT_DSP : public FFT_BASE<FPTYPE>
{
    public:
        FFT_DSP(){};
        ~FFT_DSP(){}; 
        
	    void setupFFT() override; 

        void clear() override;

        void cleanFFT() override;

        /** 
        * @brief Initialize the fft parameters
        * @param nx_in  number of grid points in x direction
        * @param ny_in  number of grid points in y direction
        * @param nz_in  number of grid points in z direction
        * 
        */
        virtual __attribute__((weak)) 
        void initfft(int nx_in, 
                     int ny_in, 
                     int nz_in) override;
        
        /**
         * @brief Get the real space data
         * @return real space data
         */
        virtual __attribute__((weak)) 
        std::complex<FPTYPE>* get_auxr_3d_data() const override;
        
        /**
         * @brief Forward FFT in 3D
         * @param in  input data, complex FPTYPE
         * @param out  output data, complex FPTYPE
         * 
         * This function performs the forward FFT in 3D.
         */
         virtual __attribute__((weak)) 
        void fft3D_forward(std::complex<FPTYPE>* in, 
                           std::complex<FPTYPE>* out) const override;
        /**
         * @brief Backward FFT in 3D
         * @param in  input data, complex FPTYPE
         * @param out  output data, complex FPTYPE
         * 
         * This function performs the backward FFT in 3D.
         */
         virtual __attribute__((weak)) 
        void fft3D_backward(std::complex<FPTYPE>* in, 
                            std::complex<FPTYPE>* out) const override;
    public:
        INT cluster_id=0;
        INT   b_id;
        INT thread_id_for=0;
        PLAN* ptr_plan_forward=nullptr;
        PLAN* ptr_plan_backward=nullptr;
        mutable unsigned long args_for[2];
        mutable unsigned long args_back[2];
        mutable E *   forward_in;
        std::complex<float>* c_auxr_3d = nullptr;  // fft space
        std::complex<double>* z_auxr_3d = nullptr; // fft space

};
void test_fft_dsp();
} // namespace ModulePW
#endif