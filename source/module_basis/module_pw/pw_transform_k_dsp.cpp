#include "module_base/timer.h"
#include "module_basis/module_pw/kernels/pw_op.h"
#include "pw_basis_k.h"

#include <cassert>
#include <complex>
#include <string>
namespace ModulePW
{
    template <typename FPTYPE>
    void PW_Basis_K::real2recip_3d(const std::complex<FPTYPE>* in,
                            std::complex<FPTYPE>* out,
                            const int ik,
                            const bool add,
                            const FPTYPE factor) const
    {
        ModuleBase::timer::tick(this->classname,"real2recip_3d");
        const base_device::DEVICE_CPU* ctx;
        const base_device::DEVICE_GPU* gpux;
        assert(this->gamma_only == false);
        auto* auxr = this->fft_bundle.get_auxr_3d_data<double>();
        
        const int startig = ik * this->npwk_max;
        const int npw_k   = this->npwk[ik];
        memcpy(auxr,in,this->nrxx*2*8);
        this->fft_bundle.fft3D_forward(gpux, 
                                       auxr,
                                       auxr);
        set_real_to_recip_output_op<double, base_device::DEVICE_CPU>()(ctx,
                                                                npw_k,
                                                                this->nxyz,
                                                                add,
                                                                factor,
                                                                this->ig2ixyz_k_cpu + startig,
                                                                this->fft_bundle.get_auxr_3d_data<double>(),
                                                                out);
        ModuleBase::timer::tick(this->classname,"real2recip_3d");
    }

    template <typename FPTYPE>
    void PW_Basis_K::recip2real_3d(const std::complex<FPTYPE>* in,
                            std::complex<FPTYPE>* out,
                            const int ik,
                            const bool add,
                            const FPTYPE factor) const
    {
        ModuleBase::timer::tick(this->classname,"recip2real_3d");

        assert(this->gamma_only == false);
        const base_device::DEVICE_CPU* ctx;
        const base_device::DEVICE_GPU* gpux;
        auto* auxr = this->fft_bundle.get_auxr_3d_data<double>();
        memset(auxr,0,this->nrxx*2*8);
        const int startig = ik * this->npwk_max;
        const int npw_k   = this->npwk[ik];

        set_3d_fft_box_op<double, base_device::DEVICE_CPU>()(ctx,
                                                        npw_k,
                                                        this->ig2ixyz_k_cpu + startig,
                                                        in,
                                                        auxr);
        this->fft_bundle.fft3D_backward(gpux,auxr,auxr);
        set_recip_to_real_output_op<double, base_device::DEVICE_CPU>()(ctx,
                                                                this->nrxx,
                                                                add,
                                                                factor,
                                                                auxr,
                                                                out);
        ModuleBase::timer::tick(this->classname,"recip2real_3d");
    }

template void PW_Basis_K::real2recip_3d<double>(const std::complex<double>* in,
                                             std::complex<double>* out,
                                             const int ik,
                                             const bool add,
                                             const double factor) const; // in:(nplane,nx*ny)  ; out(nz, ns)
template void PW_Basis_K::recip2real_3d<double>(const std::complex<double>* in,
                                             std::complex<double>* out,
                                             const int ik,
                                             const bool add,
                                             const double factor) const; // in:(nz, ns)  ; out(nplane,nx*ny)
}