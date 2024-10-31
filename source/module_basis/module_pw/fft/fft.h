#include "fft_base.h"
#ifndef FFT_H
#define FFT_H
class FFT
{
    public:
        FFT();
        FFT(std::string device_in,std::string precision_in);
        ~FFT();
        // void clear();

        void initfft(int nx_in, int ny_in, int nz_in, int lixy_in, int rixy_in, int ns_in, int nplane_in, 
                     int nproc_in, bool gamma_only_in, bool xprime_in = true, bool mpifft_in = false);
        
        void setupFFT();

        void clearFFT();
        
        void clear();
        
        template <typename FPTYPE>
        FPTYPE* get_rspace_data() const;
        template <typename FPTYPE>
        std::complex<FPTYPE>* get_auxr_data() const;
        template <typename FPTYPE>
        std::complex<FPTYPE>* get_auxg_data() const;
        template <typename FPTYPE>
        std::complex<FPTYPE>* get_auxr_3d_data() const;
        
        template <typename FPTYPE>
        void fftzfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;

        template <typename FPTYPE>
        void fftxyfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
        // template <typename FPTYPE>
        // void fftzbac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
        // template <typename FPTYPE>
        // void fftxyfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
        // template <typename FPTYPE>
        // void fftxybac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
        // template <typename FPTYPE>
        // void fftxyr2c(FPTYPE* in, std::complex<FPTYPE>* out) const;
        // template <typename FPTYPE>
        // void fftxyc2r(std::complex<FPTYPE>* in, FPTYPE* out) const;

        template <typename FPTYPE, typename Device>
        void fft3D_forward(const Device* ctx, std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
        // template <typename FPTYPE, typename Device>
        // void fft3D_backward(const Device* ctx, std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;

        void set_device(std::string device_in); 

        void set_precision(std::string precision_in);

    private:
        char fftfalg=0;
        FFT_BASE<float>* fft_float=nullptr;
        FFT_BASE<double>* fft_double=nullptr;
        
        std::string device = "cpu";
        std::string precision = "mixing";
};   

#endif // FFT_H