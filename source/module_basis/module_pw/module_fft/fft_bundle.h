#include "fft_base.h"
#include <memory>
// #include "module_psi/psi.h"
#ifndef FFT_TEMP_H
#define FFT_TEMP_H
namespace ModulePW
{
class FFT_Bundle
{
    public:
        FFT_Bundle();
        FFT_Bundle(std::string device_in,std::string precision_in);
        ~FFT_Bundle();

        void setfft(std::string device_in,std::string precision_in);
        void initfft(int nx_in, int ny_in, int nz_in, int lixy_in, int rixy_in, int ns_in, int nplane_in, 
                     int nproc_in, bool gamma_only_in, bool xprime_in = true, bool mpifft_in = false);
        
        void initfftmode(int fft_mode_in);

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
        template <typename FPTYPE>
        void fftzbac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
        template <typename FPTYPE>
        void fftxybac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
        template <typename FPTYPE>
        void fftxyr2c(FPTYPE* in, std::complex<FPTYPE>* out) const;
        template <typename FPTYPE>
        void fftxyc2r(std::complex<FPTYPE>* in, FPTYPE* out) const;

        template <typename FPTYPE, typename Device>
        void fft3D_forward(const Device* ctx, std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
        template <typename FPTYPE, typename Device>
        void fft3D_backward(const Device* ctx, std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;

        void set_device(std::string device_in); 

        void set_precision(std::string precision_in);

    private:
        int fft_mode = 0; ///< fftw mode 0: estimate, 1: measure, 2: patient, 3: exhaustive
        bool float_flag=false;
        bool float_define=false;
        bool double_flag=false;
        // FFT_BASE<float>* fft_float=nullptr; // Remove the qualified name and use a raw pointer instead
        std::shared_ptr<FFT_BASE<float>> fft_float=nullptr;
        std::shared_ptr<FFT_BASE<double>> fft_double=nullptr;
        
        std::string device = "cpu";
        std::string precision = "double";
};   
} // namespace ModulePW
#endif // FFT_H

