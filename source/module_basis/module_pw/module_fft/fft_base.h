#include <complex>
#include <string>
#include "fftw3.h"
#ifndef FFT_BASE_H
#define FFT_BASE_H
namespace ModulePW
{
template <typename FPTYPE>
class FFT_BASE
{
public:

	FFT_BASE();
	virtual  ~FFT_BASE(); 
	
	// init parameters of fft
	virtual __attribute__((weak))
    void initfft(int nx_in, int ny_in, int nz_in, int lixy_in, int rixy_in, int ns_in, int nplane_in, 
				 int nproc_in, bool gamma_only_in, bool xprime_in = true, bool mpifft_in = false);

	//init fftw_plans
	virtual void setupFFT()=0; 

	//destroy fftw_plans
	virtual void cleanFFT()=0;
    //clear fftw_data
    virtual void clear()=0;
    
    // access the real space data
    virtual __attribute__((weak)) FPTYPE* get_rspace_data() const;

    virtual __attribute__((weak)) std::complex<FPTYPE>* get_auxr_data() const;

    virtual __attribute__((weak)) std::complex<FPTYPE>* get_auxg_data() const;

    virtual __attribute__((weak)) std::complex<FPTYPE>* get_auxr_3d_data() const;

    //forward fft in x-y direction
    virtual __attribute__((weak)) void fftxyfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;

    virtual __attribute__((weak)) void fftxybac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;

    virtual __attribute__((weak)) void fftzfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
    
    virtual __attribute__((weak)) void fftzbac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;

    virtual __attribute__((weak)) void fftxyr2c(FPTYPE* in, std::complex<FPTYPE>* out) const;

    virtual __attribute__((weak)) void fftxyc2r(std::complex<FPTYPE>* in, FPTYPE* out) const;
    
    virtual __attribute__((weak)) void fft3D_forward(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;

    virtual __attribute__((weak)) void fft3D_backward(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;

protected:
	int ny=0;
    int nx=0;
    int nz=0;

};
}
#endif // FFT_BASE_H
