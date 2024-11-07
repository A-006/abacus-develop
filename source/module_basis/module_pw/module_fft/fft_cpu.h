#include "fft_base.h"
#include "fftw3.h"

// #ifdef __ENABLE_FLOAT_FFTW

// #endif
// #endif
#ifndef FFT_CPU_H
#define FFT_CPU_H
namespace ModulePW
{
template <typename FPTYPE>
class FFT_CPU : public FFT_BASE<FPTYPE>
{
    public:
    FFT_CPU();
    FFT_CPU(const int fft_mode_in);
    ~FFT_CPU(); 

    //init fftw_plans
    // __attribute__((weak)) 
    __attribute__((weak)) void initfft(int nx_in, int ny_in, int nz_in, int lixy_in, int rixy_in, int ns_in, int nplane_in, 
				            int nproc_in, bool gamma_only_in, bool xprime_in = true, bool mpifft_in = false) override;
	__attribute__((weak)) void setupFFT() override; 

	// void initplan(const unsigned int& flag = 0);
    __attribute__((weak)) void cleanFFT() override;

    __attribute__((weak)) void clear() override;

    __attribute__((weak)) FPTYPE* get_rspace_data() const override;

    __attribute__((weak)) std::complex<FPTYPE>* get_auxr_data() const;

    __attribute__((weak)) std::complex<FPTYPE>* get_auxg_data() const;

    __attribute__((weak)) void fftxyfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;

    __attribute__((weak)) void fftxybac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;

    __attribute__((weak)) void fftzfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;

    __attribute__((weak)) void fftzbac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;

    __attribute__((weak)) void fftxyr2c(FPTYPE* in, std::complex<FPTYPE>* out) const override;

    __attribute__((weak)) void fftxyc2r(std::complex<FPTYPE>* in, FPTYPE* out) const override;
    private:
        void clearfft(fftw_plan& plan);
        void clearfft(fftwf_plan& plan);

        fftw_plan planzfor  = NULL;//create a special pointer pointing to the fftw_plan class as a plan for performing FFT
        fftw_plan planzbac  = NULL;
        fftw_plan planxfor1 = NULL;
        fftw_plan planxbac1 = NULL;
        fftw_plan planxfor2 = NULL;
        fftw_plan planxbac2 = NULL;
        fftw_plan planyfor  = NULL;
        fftw_plan planybac  = NULL;
        fftw_plan planxr2c  = NULL;
        fftw_plan planxc2r  = NULL;
        fftw_plan planyr2c  = NULL;
        fftw_plan planyc2r  = NULL;

        fftwf_plan planfzfor = NULL;
        fftwf_plan planfzbac = NULL;
        fftwf_plan planfxfor1= NULL;
        fftwf_plan planfxbac1= NULL;
        fftwf_plan planfxfor2= NULL;
        fftwf_plan planfxbac2= NULL;
        fftwf_plan planfyfor = NULL;
        fftwf_plan planfybac = NULL;
        fftwf_plan planfxr2c = NULL;
        fftwf_plan planfxc2r = NULL;
        fftwf_plan planfyr2c = NULL;
        fftwf_plan planfyc2r = NULL;
        
        std::complex<float>*c_auxg = nullptr;
        std::complex<float>*c_auxr = nullptr;  // fft space
        std::complex<double>*z_auxg = nullptr;
        std::complex<double>*z_auxr = nullptr; // fft space

        float* s_rspace = nullptr;  // real number space for r, [nplane * nx *ny]
        double* d_rspace = nullptr; // real number space for r, [nplane * nx *ny]

        int initflag = 0; // 0: not initialized; 1: initialized
        int fftnx=0;
        int fftny=0;
        int fftnxy=0;
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
};
}
#endif // FFT_CPU_H