#include <complex>
#include <string>
#include "fftw3.h"
#ifndef FFT_BASE_H
#define FFT_BASE_H
template <typename FPTYPE>
class FFT_BASE
{
public:

	FFT_BASE();
	virtual ~FFT_BASE(); 
	
	// init parameters of fft
	virtual void initfft(int nx_in, int ny_in, int nz_in, int lixy_in, int rixy_in, int ns_in, int nplane_in, 
				 int nproc_in, bool gamma_only_in, bool xprime_in = true, bool mpifft_in = false);

    virtual void initfftmode(int fft_mode_in);
    
	//init fftw_plans
	virtual void setupFFT()=0; 

	//destroy fftw_plans
	virtual void cleanFFT()=0;

    virtual void clear()=0;
    
    virtual FPTYPE* get_rspace_data() const;

    virtual void fftxyfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;

    virtual void fftzfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
    
    virtual void fft3D_forward(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
    // template <typename Device>
    // virtual void fft3D_forward(const Device* ctx, std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;

protected:
    int initflag = 0; // 0: not initialized; 1: initialized
	int fftnx=0;
    int fftny=0;
	int fftnxy=0;
	int ny=0;
    int nx=0;
    int nz=0;
	int nxy=0;
    int nplane=0; //number of x-y planes
    bool gamma_only = false;
    int lixy=0;
    int rixy=0;// lixy: the left edge of the pw ball in the y direction; rixy: the right edge of the pw ball in the x or y direction
    bool mpifft = false; // if use mpi fft, only used when define __FFTW3_MPI
    int maxgrids = 0; // maxgrids = (nsz > nrxx) ? nsz : nrxx;
    bool xprime = true; // true: when do recip2real, x-fft will be done last and when doing real2recip, x-fft will be done first; false: y-fft
                         // For gamma_only, true: we use half x; false: we use half y
	int ns=0; //number of sticks
	int nproc=1; // number of proc.
	int fft_mode = 0; ///< fftw mode 0: estimate, 1: measure, 2: patient, 3: exhaustive 
    
public:
    void set_device(std::string device_);
    void set_precision(std::string precision_);

};
#endif // FFT_BASE_H
