#include "fft_base.h"
template <typename FPTYPE>
FFT_BASE<FPTYPE>::FFT_BASE()
{
}
template <typename FPTYPE>
FFT_BASE<FPTYPE>::~FFT_BASE()
{

}
template <typename FPTYPE>
void FFT_BASE<FPTYPE>::initfft(int nx_in, int ny_in, int nz_in, int lixy_in, int rixy_in, int ns_in, int nplane_in, 
				 int nproc_in, bool gamma_only_in, bool xprime_in, bool mpifft_in)
{
    this->gamma_only = gamma_only_in;
    this->xprime = xprime_in;
    this->fftnx = this->nx = nx_in;
    this->fftny = this->ny = ny_in;
    if (this->gamma_only)
    {
        if (xprime)
            this->fftnx = int(nx / 2) + 1;
        else
            this->fftny = int(ny / 2) + 1;
    }
    this->nz = nz_in;
    this->ns = ns_in;
    this->lixy = lixy_in;
    this->rixy = rixy_in;
    this->nplane = nplane_in;
    this->nproc = nproc_in;
    this->mpifft = mpifft_in;
    this->nxy = this->nx * this->ny;
    this->fftnxy = this->fftnx * this->fftny;
    const int nrxx = this->nxy * this->nplane;
    const int nsz = this->nz * this->ns;
    this->maxgrids = (nsz > nrxx) ? nsz : nrxx;
}

template FFT_BASE<float>::FFT_BASE();
template FFT_BASE<double>::FFT_BASE();
template FFT_BASE<float>::~FFT_BASE();
template FFT_BASE<double>::~FFT_BASE();