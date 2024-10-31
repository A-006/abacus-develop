#include "fft_cpu.h"
#include "fftw3.h"

template <>
FFT_CPU<double>::FFT_CPU()
{
    
}
template <>
FFT_CPU<double>::~FFT_CPU()
{

}
template <>
void FFT_CPU<double>::initfftmode(int fft_mode_in)
{
    this->fft_mode = fft_mode_in;
}
template <>
void FFT_CPU<double>::setupFFT()
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
        z_auxg = (std::complex<double>*)fftw_malloc(sizeof(fftw_complex) * this->maxgrids);
        z_auxr = (std::complex<double>*)fftw_malloc(sizeof(fftw_complex) * this->maxgrids);
        d_rspace = (double*)z_auxg;
        this->planzfor = fftw_plan_many_dft(1, &this->nz, this->ns, (fftw_complex*)z_auxg, &this->nz, 1, this->nz,
                                        (fftw_complex*)z_auxg, &this->nz, 1, this->nz, FFTW_FORWARD, flag);

        this->planzbac = fftw_plan_many_dft(1, &this->nz, this->ns, (fftw_complex*)z_auxg, &this->nz, 1, this->nz,
                                            (fftw_complex*)z_auxg, &this->nz, 1, this->nz, FFTW_BACKWARD, flag);

        //---------------------------------------------------------
        //                              2 D - XY
        //---------------------------------------------------------
        // 1D+1D is much faster than 2D FFT!
        // in-place fft is better for c2c and out-of-place fft is better for c2r
        int* embed = nullptr;
        int npy = this->nplane * this->ny;
        if (this->xprime)
        {
            this->planyfor = fftw_plan_many_dft(1, &this->ny, this->nplane, (fftw_complex*)z_auxr, embed,this->nplane, 1,
                                                (fftw_complex*)z_auxr, embed,this->nplane, 1, FFTW_FORWARD, flag);
            this->planybac = fftw_plan_many_dft(1, &this->ny, this->nplane, (fftw_complex*)z_auxr, embed,this->nplane, 1,
                                                (fftw_complex*)z_auxr, embed,this->nplane, 1, FFTW_BACKWARD, flag);
            if (this->gamma_only)
            {
                this->planxr2c = fftw_plan_many_dft_r2c(1, &this->nx, npy, d_rspace, embed, npy, 1, (fftw_complex*)z_auxr,
                                                        embed, npy, 1, flag);
                this->planxc2r = fftw_plan_many_dft_c2r(1, &this->nx, npy, (fftw_complex*)z_auxr, embed, npy, 1, d_rspace,
                                                        embed, npy, 1, flag);
            }
            else
            {
                this->planxfor1 = fftw_plan_many_dft(1, &this->nx, npy, (fftw_complex*)z_auxr, embed, npy, 1,
                                                    (fftw_complex*)z_auxr, embed, npy, 1, FFTW_FORWARD, flag);
                this->planxbac1 = fftw_plan_many_dft(1, &this->nx, npy, (fftw_complex*)z_auxr, embed, npy, 1,
                                                    (fftw_complex*)z_auxr, embed, npy, 1, FFTW_BACKWARD, flag);
            }
        }
        else
        {
            this->planxfor1 = fftw_plan_many_dft(1, &this->nx, this->nplane * (this->lixy + 1), (fftw_complex*)z_auxr, embed, npy,
                                                1, (fftw_complex*)z_auxr, embed, npy, 1, FFTW_FORWARD, flag);
            this->planxbac1 = fftw_plan_many_dft(1, &this->nx, this->nplane * (this->lixy + 1), (fftw_complex*)z_auxr, embed, npy,
                                                1, (fftw_complex*)z_auxr, embed, npy, 1, FFTW_BACKWARD, flag);
            if (this->gamma_only)
            {
                this->planyr2c = fftw_plan_many_dft_r2c(1, &this->ny, this->nplane, d_rspace, embed, this->nplane, 1,
                                                        (fftw_complex*)z_auxr, embed, this->nplane, 1, flag);
                this->planyc2r = fftw_plan_many_dft_c2r(1, &this->ny, this->nplane, (fftw_complex*)z_auxr, embed,
                                                        this->nplane, 1, d_rspace, embed, this->nplane, 1, flag);
            }
            else
            {

                this->planxfor2 = fftw_plan_many_dft(1, &this->nx, this->nplane * (this->ny - this->rixy), (fftw_complex*)z_auxr, embed,
                                                    npy, 1, (fftw_complex*)z_auxr, embed, npy, 1, FFTW_FORWARD, flag);
                this->planxbac2 = fftw_plan_many_dft(1, &this->nx, this->nplane * (this->ny - this->rixy), (fftw_complex*)z_auxr, embed,
                                                    npy, 1, (fftw_complex*)z_auxr, embed, npy, 1, FFTW_BACKWARD, flag);
                this->planyfor = fftw_plan_many_dft(1, &this->ny, this->nplane, (fftw_complex*)z_auxr, embed, this->nplane,
                                                    1, (fftw_complex*)z_auxr, embed, this->nplane, 1, FFTW_FORWARD, flag);
                this->planybac = fftw_plan_many_dft(1, &this->ny, this->nplane, (fftw_complex*)z_auxr, embed, this->nplane,
                                                    1, (fftw_complex*)z_auxr, embed, this->nplane, 1, FFTW_BACKWARD, flag);
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
void FFT_CPU<double>::clearfft(fftw_plan& plan)
{
    if (plan)
    {
        fftw_destroy_plan(plan);
        plan = NULL;
    }
}

template <>
void FFT_CPU<double>::cleanFFT()
{
    printf("in the double cleanFFT\n");
    clearfft(planzfor);
    clearfft(planzbac);
    clearfft(planxfor1);
    clearfft(planxbac1);
    clearfft(planxfor2);
    clearfft(planxbac2);
    clearfft(planyfor);
    clearfft(planybac);
    clearfft(planxr2c);
    clearfft(planxc2r);
    clearfft(planyr2c);
    clearfft(planyc2r);
}


template <>
void FFT_CPU<double>::clear()
{
    this->cleanFFT();
    if (z_auxg != nullptr)
    {
        fftw_free(z_auxg);
        z_auxg = nullptr;
    }
    if (z_auxr != nullptr)
    {
        fftw_free(z_auxr);
        z_auxr = nullptr;
    }
    d_rspace = nullptr;
}

template <> 
float* FFT_CPU<float>::get_rspace_data() const
{
    return s_rspace;
}

template <>
double* FFT_CPU<double>::get_rspace_data() const
{
    return d_rspace;
}

template <>
void FFT_CPU<double>::fftxyfor(std::complex<double>* in, std::complex<double>* out) const
{
    int npy = this->nplane * this->ny;
    if (this->xprime)
    {
        fftw_execute_dft(this->planxfor1, (fftw_complex*)in, (fftw_complex*)out);
        printf("the first element of out is %f\n",out[0].real());
        for (int i = 0; i < this->lixy + 1; ++i)
        {
            fftw_execute_dft(this->planyfor, (fftw_complex*)&in[i * npy], (fftw_complex*)&out[i * npy]);
        }
        for (int i = rixy; i < this->nx; ++i)
        {
            fftw_execute_dft(this->planyfor, (fftw_complex*)&in[i * npy], (fftw_complex*)&out[i * npy]);
        }
    }
    else
    {
        for (int i = 0; i < this->nx; ++i)
        {
            fftw_execute_dft(this->planyfor, (fftw_complex*)&in[i * npy], (fftw_complex*)&out[i * npy]);
        }

        fftw_execute_dft(this->planxfor1, (fftw_complex*)in, (fftw_complex*)out);
        fftw_execute_dft(this->planxfor2, (fftw_complex*)&in[rixy * nplane], (fftw_complex*)&out[rixy * nplane]);
    }
}

template <>
void FFT_CPU<double>::fftzfor(std::complex<double>* in, std::complex<double>* out) const
{
    fftw_execute_dft(this->planzfor, (fftw_complex*)in, (fftw_complex*)out);
}
template FFT_CPU<double>::FFT_CPU();
