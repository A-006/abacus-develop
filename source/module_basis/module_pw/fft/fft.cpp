#include <cassert>
#include "fft.h"
#include "fft_cpu.h"
// #include "fft_gpu.h"

FFT::FFT(std::string device_in,std::string precision_in)
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
        
        // fft_float = new FFT_GPU<float>();
        // fft_double = new FFT_GPU<double>();
    }
}

FFT::~FFT()
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

void FFT::set_device(std::string device_in)
{
    this->device = device_in;
}

void FFT::set_precision(std::string precision_in)
{
    this->precision = precision_in;
}

void FFT::initfft(int nx_in, int ny_in, int nz_in, int lixy_in, int rixy_in, int ns_in, int nplane_in, 
                     int nproc_in, bool gamma_only_in, bool xprime_in , bool mpifft_in)
{

    fft_float->initfft(nx_in,ny_in,nz_in,lixy_in,rixy_in,ns_in,nplane_in,nproc_in,gamma_only_in,xprime_in,mpifft_in);
    fft_float->initfftmode(this->fft_mode);
    fft_double->initfft(nx_in,ny_in,nz_in,lixy_in,rixy_in,ns_in,nplane_in,nproc_in,gamma_only_in,xprime_in,mpifft_in);
    fft_double->initfftmode(this->fft_mode);
}
void FFT::initfftmode(int fft_mode_in)
{
    this->fft_mode = fft_mode_in;
}

void FFT::setupFFT()
{
    if (precision=="single")
    {
        fft_float->setupFFT();
    }
    else if (precision=="double")
    {
        fft_double->setupFFT();
    }
    else if (precision=="mixing")
    {
        fft_float->setupFFT();
        fft_double->setupFFT();
    }
}

void FFT::clearFFT()
{
    if (precision=="single")
    {
        fft_float->clear();
    }
    else if (precision=="double")
    {
        fft_double->clear();
    }
    else if (precision=="mixing")
    {
        fft_float->clear();
        fft_double->clear();
    }
}

template <>
float* FFT::get_rspace_data() const
{
    return fft_float->get_rspace_data();
}

template <>
double* FFT::get_rspace_data() const
{
    return fft_double->get_rspace_data();
}

template <>
void FFT::fftxyfor(std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fftxyfor(in,out);
}

template <>
void FFT::fftxyfor(std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fftxyfor(in,out);
}

template <>
void FFT::fftzfor(std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fftzfor(in,out);
}

template <>
void FFT::fftzfor(std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fftzfor(in,out);
}

template <>
void  FFT::fft3D_forward(const char* ctx, std::complex<float>* in, std::complex<float>* out) const
{
    fft_float->fft3D_forward(in, out);
}

template <>
void  FFT::fft3D_forward(const char* ctx, std::complex<double>* in, std::complex<double>* out) const
{
    fft_double->fft3D_forward(in, out);
}