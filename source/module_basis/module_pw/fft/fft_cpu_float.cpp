#include "fft_cpu.h"
#include "fftw3.h"

template <>
FFT_CPU<float>::FFT_CPU()
{
    
}
template <>
FFT_CPU<float>::~FFT_CPU()
{

}
template <>
void FFT_CPU<float>::initfftmode(int fft_mode_in)
{
    this->fft_mode = fft_mode_in;
}
template <>
void FFT_CPU<float>::setupFFT()
{
    unsigned int flag = FFTW_ESTIMATE;
    switch (this->fft_mode)
    {
    case 0:
        flag = FFTW_ESTIMATE;
        break;
    case 1:
        flag = FFTW_MEASURE;
        break;
    case 2:
        flag = FFTW_PATIENT;
        break;
    case 3:
        flag = FFTW_EXHAUSTIVE;
        break;
    default:
        break;
    }
    if (!this->mpifft)
    {
        c_auxg = (std::complex<float>*)fftwf_malloc(sizeof(fftwf_complex) * this->maxgrids); 
        c_auxr = (std::complex<float>*)fftwf_malloc(sizeof(fftwf_complex) * maxgrids);
        s_rspace = (float*)c_auxg;
        //---------------------------------------------------------
        //                              1 D
        //---------------------------------------------------------

        //               fftw_plan_many_dft(int rank,          const int *n,       int howmany,
        //					                fftw_complex *in,  const int *inembed, int istride, int idist,
        //					                fftw_complex *out, const int *onembed, int ostride, int odist, int sign, unsigned
        //flags);

        this->planfzfor = fftwf_plan_many_dft(1, &this->nz, this->ns, (fftwf_complex*)c_auxg, &this->nz, 1, this->nz,
                                            (fftwf_complex*)c_auxg, &this->nz, 1, this->nz, FFTW_FORWARD, flag);

        this->planfzbac = fftwf_plan_many_dft(1, &this->nz, this->ns, (fftwf_complex*)c_auxg, &this->nz, 1, this->nz,
                                            (fftwf_complex*)c_auxg, &this->nz, 1, this->nz, FFTW_BACKWARD, flag);
        //---------------------------------------------------------
        //                              2 D
        //---------------------------------------------------------

        int* embed = nullptr;
        int npy = this->nplane * this->ny;
        if (this->xprime)
        {
            this->planfyfor = fftwf_plan_many_dft(1, &this->ny, this->nplane, (fftwf_complex*)c_auxr, embed, nplane, 1,
                                                (fftwf_complex*)c_auxr, embed, nplane, 1, FFTW_FORWARD, flag);
            this->planfybac = fftwf_plan_many_dft(1, &this->ny, this->nplane, (fftwf_complex*)c_auxr, embed, nplane, 1,
                                                (fftwf_complex*)c_auxr, embed, nplane, 1, FFTW_BACKWARD, flag);
            if (this->gamma_only)
            {
                this->planfxr2c = fftwf_plan_many_dft_r2c(1, &this->nx, npy, s_rspace, embed, npy, 1,
                                                        (fftwf_complex*)c_auxr, embed, npy, 1, flag);
                this->planfxc2r = fftwf_plan_many_dft_c2r(1, &this->nx, npy, (fftwf_complex*)c_auxr, embed, npy, 1,
                                                        s_rspace, embed, npy, 1, flag);
            }
            else
            {
                this->planfxfor1 = fftwf_plan_many_dft(1, &this->nx, npy, (fftwf_complex*)c_auxr, embed, npy, 1,
                                                    (fftwf_complex*)c_auxr, embed, npy, 1, FFTW_FORWARD, flag);
                this->planfxbac1 = fftwf_plan_many_dft(1, &this->nx, npy, (fftwf_complex*)c_auxr, embed, npy, 1,
                                                    (fftwf_complex*)c_auxr, embed, npy, 1, FFTW_BACKWARD, flag);
            }
        }
        else
        {
            this->planfxfor1 = fftwf_plan_many_dft(1, &this->nx, this->nplane * (lixy + 1), (fftwf_complex*)c_auxr, embed,
                                                npy, 1, (fftwf_complex*)c_auxr, embed, npy, 1, FFTW_FORWARD, flag);
            this->planfxbac1 = fftwf_plan_many_dft(1, &this->nx, this->nplane * (lixy + 1), (fftwf_complex*)c_auxr, embed,
                                                npy, 1, (fftwf_complex*)c_auxr, embed, npy, 1, FFTW_BACKWARD, flag);
            if (this->gamma_only)
            {
                this->planfyr2c = fftwf_plan_many_dft_r2c(1, &this->ny, this->nplane, s_rspace, embed, this->nplane, 1,
                                                        (fftwf_complex*)c_auxr, embed, this->nplane, 1, flag);
                this->planfyc2r = fftwf_plan_many_dft_c2r(1, &this->ny, this->nplane, (fftwf_complex*)c_auxr, embed,
                                                        this->nplane, 1, s_rspace, embed, this->nplane, 1, flag);
            }
            else
            {
                this->planfxfor2
                    = fftwf_plan_many_dft(1, &this->nx, this->nplane * (this->ny - rixy), (fftwf_complex*)c_auxr, embed,
                                        npy, 1, (fftwf_complex*)c_auxr, embed, npy, 1, FFTW_FORWARD, flag);
                this->planfxbac2
                    = fftwf_plan_many_dft(1, &this->nx, this->nplane * (this->ny - rixy), (fftwf_complex*)c_auxr, embed,
                                        npy, 1, (fftwf_complex*)c_auxr, embed, npy, 1, FFTW_BACKWARD, flag);
                this->planfyfor
                    = fftwf_plan_many_dft(1, &this->ny, this->nplane, (fftwf_complex*)c_auxr, embed, this->nplane, 1,
                                        (fftwf_complex*)c_auxr, embed, this->nplane, 1, FFTW_FORWARD, flag);
                this->planfybac
                    = fftwf_plan_many_dft(1, &this->ny, this->nplane, (fftwf_complex*)c_auxr, embed, this->nplane, 1,
                                        (fftwf_complex*)c_auxr, embed, this->nplane, 1, FFTW_BACKWARD, flag);
            }
        }
        }
    #if defined(__FFTW3_MPI) && defined(__MPI)
        else
        {
            // this->initplan_mpi();
            // if (this->precision == "single") {
            //     this->initplanf_mpi();
            // }
        }
    #endif
    return;
}

template <>
void FFT_CPU<float>::fftxyfor(std::complex<float>* in, std::complex<float>* out) const
{
    int npy = this->nplane * this->ny;
    if (this->xprime)
    {
        fftwf_execute_dft(this->planfxfor1, (fftwf_complex*)in, (fftwf_complex*)out);

        for (int i = 0; i < this->lixy + 1; ++i)
        {
            fftwf_execute_dft(this->planfyfor, (fftwf_complex*)&in[i * npy], (fftwf_complex*)&out[i * npy]);
        }
        for (int i = rixy; i < this->nx; ++i)
        {
            fftwf_execute_dft(this->planfyfor, (fftwf_complex*)&in[i * npy], (fftwf_complex*)&out[i * npy]);
        }
    }
    else
    {
        for (int i = 0; i < this->nx; ++i)
        {
            fftwf_execute_dft(this->planfyfor, (fftwf_complex*)&in[i * npy], (fftwf_complex*)&out[i * npy]);
        }

        fftwf_execute_dft(this->planfxfor1, (fftwf_complex*)in, (fftwf_complex*)out);
        fftwf_execute_dft(this->planfxfor2, (fftwf_complex*)&in[rixy * nplane], (fftwf_complex*)&out[rixy * nplane]);
    }
}

template <>
void FFT_CPU<float>::fftzfor(std::complex<float>* in, std::complex<float>* out) const
{
    fftwf_execute_dft(this->planfzfor, (fftwf_complex*)in, (fftwf_complex*)out);
}

template <>
void FFT_CPU<float>::clearfft(fftwf_plan& plan)
{
    if (plan)
    {
        fftwf_destroy_plan(plan);
        plan = NULL;
    }
}

template <>
void FFT_CPU<float>::cleanFFT()
{
    printf("in the float cleanFFT\n");
    clearfft(planfzfor);
    clearfft(planfzbac);
    clearfft(planfxfor1);
    clearfft(planfxbac1);
    clearfft(planfxfor2);
    clearfft(planfxbac2);
    clearfft(planfyfor);
    clearfft(planfybac);
    clearfft(planfxr2c);
    clearfft(planfxc2r);
    clearfft(planfyr2c);
    clearfft(planfyc2r);
}

template <>
void FFT_CPU<float>::clear()
{
    this->cleanFFT();
    if (c_auxg != nullptr)
    {
        fftw_free(c_auxg);
        c_auxg = nullptr;
    }
    if (c_auxr != nullptr)
    {
        fftwf_free(c_auxr);
        c_auxr = nullptr;
    }
    d_rspace = nullptr;
}

template FFT_CPU<float>::FFT_CPU();