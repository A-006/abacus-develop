#include <cassert>
#include "fft.h"
#include "fft_cpu.h"
#if defined(__CUDA)
#include "fft_cuda.h"
#endif
#if defined(__ROCM)
#include "fft_rcom.h"
#endif
#include "module_base/module_device/device.h"
// #include "fft_gpu.h"
FFT1::FFT1()
{
    fft_float = nullptr;
    fft_double = nullptr;
}
FFT1::FFT1(std::string device_in,std::string precision_in)
{
    assert(device_in=="cpu" || device_in=="gpu");
    assert(precision_in=="single" || precision_in=="double" || precision_in=="mixing");
    this->device = device_in;
    this->precision = precision_in;
    if (device=="cpu")
    {
        fft_float = new FFT_CPU<float>();
        fft_double = new FFT_CPU<double>();
    }
    else if (device=="gpu")
    {        
        #if defined(__ROCM)
            fft_float = new FFT_RCOM<float>();
            fft_double = new FFT_RCOM<double>();
        #elif defined(__CUDA)
            fft_float = new FFT_CUDA<float>();
            fft_double = new FFT_CUDA<double>();
        #endif
    }
}

FFT1::~FFT1()
{
    if (fft_float!=nullptr)
    {
        delete fft_float;
        fft_float=nullptr;
    }
    if (fft_double!=nullptr)
    {
        delete fft_double;
        fft_double=nullptr;
    }
}

void FFT1::set_device(std::string device_in)
{
    this->device = device_in;
}

void FFT1::set_precision(std::string precision_in)
{
    this->precision = precision_in;
}
void FFT1::setfft(std::string device_in,std::string precision_in)
{
    assert(device_in=="cpu" || device_in=="gpu");
    assert(precision_in=="single" || precision_in=="double" || precision_in=="mixing");
    this->device = device_in;
    this->precision = precision_in;
    if (device=="cpu")
    {
        fft_float = new FFT_CPU<float>();
        fft_double = new FFT_CPU<double>();
    }
    else if (device=="gpu")
    {      
        #if defined(__ROCM)
            fft_float = new FFT_RCOM<float>();
            fft_double = new FFT_RCOM<double>();
        #elif defined(__CUDA)
            fft_float = new FFT_CUDA<float>();
            fft_double = new FFT_CUDA<double>();
        #endif
    }
}
void FFT1::initfft(int nx_in, int ny_in, int nz_in, int lixy_in, int rixy_in, int ns_in, int nplane_in, 
                     int nproc_in, bool gamma_only_in, bool xprime_in , bool mpifft_in)
{
    if (this->precision=="single")
    {
        float_flag = 1;
    }
    else if (this->precision=="double")
    {
        double_flag = 1;
    }
    else if (this->precision=="mixing")
    {
        float_flag = 1;
        double_flag = 1;
    }
    if (float_flag)
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
void FFT1::initfftmode(int fft_mode_in)
{
    this->fft_mode = fft_mode_in;
}

void FFT1::setupFFT()
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

void FFT1::clearFFT()
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
void FFT1::clear()
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
}
// access the real space data
template <>
float* FFT1::get_rspace_data() const
{
    return fft_float->get_rspace_data();
}

template <>
double* FFT1::get_rspace_data() const
{
    return fft_double->get_rspace_data();
}
template <>
std::complex<float>* FFT1::get_auxr_data() const
{
    return fft_float->get_auxr_data();
}
template <>
std::complex<double>* FFT1::get_auxr_data() const
{
    return fft_double->get_auxr_data();
}
template <>
std::complex<float>* FFT1::get_auxg_data() const
{
    return fft_float->get_auxg_data();
}
template <>
std::complex<double>* FFT1::get_auxg_data() const
{
    return fft_double->get_auxg_data();
}
template <>
std::complex<float>* FFT1::get_auxr_3d_data() const
{
    return fft_float->get_auxr_3d_data();
}
template <>
std::complex<double>* FFT1::get_auxr_3d_data() const
{
    return fft_double->get_auxr_3d_data();
}
template <>
void FFT1::fftxyfor(std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fftxyfor(in,out);
}

template <>
void FFT1::fftxyfor(std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fftxyfor(in,out);
}

template <>
void FFT1::fftzfor(std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fftzfor(in,out);
}
template <>
void FFT1::fftzfor(std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fftzfor(in,out);
}

template <>
void FFT1::fftxybac(std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fftxybac(in,out);
}
template <>
void FFT1::fftxybac(std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fftxybac(in,out);
}

template <>
void FFT1::fftzbac(std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fftzbac(in,out);
}
template <>
void FFT1::fftzbac(std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fftzbac(in,out);
}
template <>
void FFT1::fftxyr2c(float* in, std::complex<float>* out) const
{
    fft_float->fftxyr2c(in,out);
}
template <>
void FFT1::fftxyr2c(double* in, std::complex<double>* out) const
{
    fft_double->fftxyr2c(in,out);
}

template <>
void FFT1::fftxyc2r(std::complex<float>* in, float* out) const
{
    fft_float->fftxyc2r(in,out);
}
template <>
void FFT1::fftxyc2r(std::complex<double>* in, double* out) const
{
    fft_double->fftxyc2r(in,out);
}

template <>
void  FFT1::fft3D_forward(const base_device::DEVICE_GPU* ctx, std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fft3D_forward(in, out);
}

template <>
void  FFT1::fft3D_forward(const base_device::DEVICE_GPU* ctx, std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fft3D_forward(in, out);
}
template <>
void  FFT1::fft3D_backward(const base_device::DEVICE_GPU* ctx, std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fft3D_backward(in, out);
}
template <>
void  FFT1::fft3D_backward(const base_device::DEVICE_GPU* ctx, std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fft3D_backward(in, out);
}