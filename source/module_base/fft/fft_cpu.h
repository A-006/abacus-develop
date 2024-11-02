#include "fft_base.h"
#include "fftw3.h"

// #ifdef __ENABLE_FLOAT_FFTW

// #endif
// #endif
#ifndef FFT_CPU_H
#define FFT_CPU_H

template <typename FPTYPE>
class FFT_CPU : public FFT_BASE<FPTYPE>
{
    public:
    FFT_CPU();
    ~FFT_CPU(); 

    void initfftmode(int fft_mode_in) override;

    //init fftw_plans
	void setupFFT() override; 

	// void initplan(const unsigned int& flag = 0);
    void cleanFFT() override;

    void clear() override;

    FPTYPE* get_rspace_data() const override;

    std::complex<FPTYPE>* get_auxr_data() const;

    std::complex<FPTYPE>* get_auxg_data() const;

    void fftxyfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;

    void fftxybac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;

    void fftzfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;

    void fftzbac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;

    void fftxyr2c(FPTYPE* in, std::complex<FPTYPE>* out) const override;

    void fftxyc2r(std::complex<FPTYPE>* in, FPTYPE* out) const override;
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
};
#endif // FFT_CPU_H