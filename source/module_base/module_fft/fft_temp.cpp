#include <cassert>
#include "fft_temp.h"
#include "fft_cpu.h"
#include "module_base/module_device/device.h"
// #if defined(__CUDA)
// #include "fft_cuda.h"
// #endif
// #if defined(__ROCM)
// #include "fft_rcom.h"
// #endif

template<typename FFT_BASE, typename... Args>
std::unique_ptr<FFT_BASE> make_unique(Args &&... args)
{
    return std::unique_ptr<FFT_BASE>(new FFT_BASE(std::forward<Args>(args)...));
}
// #include "fft_gpu.h"
FFT_TEMP::FFT_TEMP()
{
}
FFT_TEMP::FFT_TEMP(std::string device_in,std::string precision_in)
{
    assert(device_in=="cpu" || device_in=="gpu");
    assert(precision_in=="single" || precision_in=="double" || precision_in=="mixing");
    this->device = device_in;
    this->precision = precision_in;
    if (device=="cpu")
    {
        fft_float = make_unique<FFT_CPU<float>>();
        fft_double = make_unique<FFT_CPU<double>>();
    }
    // else if (device=="gpu")
    // {        
    //     #if defined(__ROCM)
    //         fft_float = new FFT_RCOM<float>();
    //         fft_double = new FFT_RCOM<double>();
    //     #elif defined(__CUDA)
    //         fft_float = new FFT_CUDA<float>();
    //         fft_double = new FFT_CUDA<double>();
    //     #endif
    // }
}

FFT_TEMP::~FFT_TEMP()
{
}

void FFT_TEMP::set_device(std::string device_in)
{
    this->device = device_in;
}

void FFT_TEMP::set_precision(std::string precision_in)
{
    this->precision = precision_in;
}
void FFT_TEMP::setfft(std::string device_in,std::string precision_in)
{
    assert(device_in=="cpu" || device_in=="gpu");
    assert(precision_in=="single" || precision_in=="double" || precision_in=="mixing");
    this->device = device_in;
    this->precision = precision_in;

}
void FFT_TEMP::initfft(int nx_in, int ny_in, int nz_in, int lixy_in, int rixy_in, int ns_in, int nplane_in, 
                     int nproc_in, bool gamma_only_in, bool xprime_in , bool mpifft_in)
{
    if (device=="cpu")
    {
        fft_float = make_unique<FFT_CPU<float>>();
        fft_double = make_unique<FFT_CPU<double>>();
        // fft_double = new FFT_CPU<double>();
        #ifndef __ENABLE_FLOAT_FFTW
            float_define = false;
        #endif
    }
    if (device=="gpu")
    {
        // #if defined(__ROCM)
        //     fft_float = new FFT_RCOM<float>();
        //     fft_double = new FFT_RCOM<double>();
        // #elif defined(__CUDA)
        //     fft_float = make_unique<FFT_CUDA<float>>();
        //     fft_double = make_unique<FFT_CUDA<double>>();
        // #endif
    }
    if (this->precision=="single")
    {
        #ifndef __ENABLE_FLOAT_FFTW
        float_define = false;
        #endif
        float_flag = float_flag;
        double_flag = true;
        
        
    }
    else if (this->precision=="double")
    {
        double_flag = true;
    }

    if (float_flag && float_define)
    {
        fft_float->initfftmode(this->fft_mode);
        fft_float->initfft(nx_in,ny_in,nz_in,lixy_in,rixy_in,ns_in,nplane_in,nproc_in,gamma_only_in,xprime_in,mpifft_in);
    }
    if (double_flag)
    {
        fft_double->initfftmode(this->fft_mode);
        fft_double->initfft(nx_in,ny_in,nz_in,lixy_in,rixy_in,ns_in,nplane_in,nproc_in,gamma_only_in,xprime_in,mpifft_in);
    }
}
void FFT_TEMP::initfftmode(int fft_mode_in)
{
    this->fft_mode = fft_mode_in;
}

void FFT_TEMP::setupFFT()
{
    if (double_flag)
    {
        fft_double->setupFFT();
    }
    if (float_flag)
    {
        fft_float->setupFFT();
    }
}

void FFT_TEMP::clearFFT()
{
    if (double_flag)
    {
        fft_double->cleanFFT();
    }
    if (float_flag)
    {
        fft_float->cleanFFT();
    }
}
void FFT_TEMP::clear()
{
    this->clearFFT();
    if (float_flag)
    {
        fft_float->clear();
    }
    if (double_flag)
    {
        fft_double->clear();
    }
    // if (fft_float!=nullptr)
    // {
    //     delete fft_float;
    //     fft_float=nullptr;
    //     float_flag = false;
    // }
    // if (fft_double!=nullptr)
    // {
    //     delete fft_double;
    //     fft_double=nullptr;
    //     double_flag = false;
    // }
}
// access the real space data
template <>
float* FFT_TEMP::get_rspace_data() const
{
    return fft_float->get_rspace_data();
}

template <>
double* FFT_TEMP::get_rspace_data() const
{
    return fft_double->get_rspace_data();
}
template <>
std::complex<float>* FFT_TEMP::get_auxr_data() const
{
    return fft_float->get_auxr_data();
}
template <>
std::complex<double>* FFT_TEMP::get_auxr_data() const
{
    return fft_double->get_auxr_data();
}
template <>
std::complex<float>* FFT_TEMP::get_auxg_data() const
{
    return fft_float->get_auxg_data();
}
template <>
std::complex<double>* FFT_TEMP::get_auxg_data() const
{
    return fft_double->get_auxg_data();
}
template <>
std::complex<float>* FFT_TEMP::get_auxr_3d_data() const
{
    return fft_float->get_auxr_3d_data();
}
template <>
std::complex<double>* FFT_TEMP::get_auxr_3d_data() const
{
    return fft_double->get_auxr_3d_data();
}
template <>
void FFT_TEMP::fftxyfor(std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fftxyfor(in,out);
}

template <>
void FFT_TEMP::fftxyfor(std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fftxyfor(in,out);
}

template <>
void FFT_TEMP::fftzfor(std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fftzfor(in,out);
}
template <>
void FFT_TEMP::fftzfor(std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fftzfor(in,out);
}

template <>
void FFT_TEMP::fftxybac(std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fftxybac(in,out);
}
template <>
void FFT_TEMP::fftxybac(std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fftxybac(in,out);
}

template <>
void FFT_TEMP::fftzbac(std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fftzbac(in,out);
}
template <>
void FFT_TEMP::fftzbac(std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fftzbac(in,out);
}
template <>
void FFT_TEMP::fftxyr2c(float* in, std::complex<float>* out) const
{
    fft_float->fftxyr2c(in,out);
}
template <>
void FFT_TEMP::fftxyr2c(double* in, std::complex<double>* out) const
{
    fft_double->fftxyr2c(in,out);
}

template <>
void FFT_TEMP::fftxyc2r(std::complex<float>* in, float* out) const
{
    fft_float->fftxyc2r(in,out);
}
template <>
void FFT_TEMP::fftxyc2r(std::complex<double>* in, double* out) const
{
    fft_double->fftxyc2r(in,out);
}

template <>
void  FFT_TEMP::fft3D_forward(const base_device::DEVICE_GPU* ctx, std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fft3D_forward(in, out);
}

template <>
void  FFT_TEMP::fft3D_forward(const base_device::DEVICE_GPU* ctx, std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fft3D_forward(in, out);
}
template <>
void  FFT_TEMP::fft3D_backward(const base_device::DEVICE_GPU* ctx, std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fft3D_backward(in, out);
}
template <>
void  FFT_TEMP::fft3D_backward(const base_device::DEVICE_GPU* ctx, std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fft3D_backward(in, out);
}