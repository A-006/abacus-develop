#include "gint.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

void Gint::cpu_vlocal_interface(Gint_inout* inout)
{
    ModuleBase::TITLE("Gint_interface", "cal_gint_vlocal");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * ucell.nwmax;
    const int lgd = this->gridt->lgd;
    const int nnrg = this->gridt->nnrg;
    const int ncyz = this->ny * this->nplane; 
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;
    double* pvpR_thread = nullptr;
    hamilt::HContainer<double>* hRGint_thread = nullptr; 

    if (!pvpR_alloc_flag)
    {
        ModuleBase::WARNING_QUIT("Gint_interface::cal_gint", "pvpR has not been allocated yet!");
    }
    else
    {
        ModuleBase::GlobalFunc::ZEROS(this->pvpR_reduced[inout->ispin], nnrg);
    }

    if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
        pvpR_thread = new double[nnrg]();
        ModuleBase::GlobalFunc::ZEROS(pvpR_thread, nnrg);
    }

    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        hRGint_thread = new hamilt::HContainer<double>(*this->hRGint);
    }

#ifdef _OPENMP
#pragma omp for
#endif
    for (int grid_index = 0; grid_index < this->nbxx; grid_index++)
    {
        const int na_grid = this->gridt->how_many_atoms[grid_index];
        if (na_grid == 0)
        {
            continue;
        }
        double* vldr3 = Gint_Tools::get_vldr3(inout->vl,
                                                this->bxyz,
                                                this->bx,
                                                this->by,
                                                this->bz,
                                                this->nplane,
                                                this->gridt->start_ind[grid_index],
                                                ncyz,
                                                dv);
#ifdef _OPENMP
        if ((GlobalV::GAMMA_ONLY_LOCAL && lgd > 0) || !GlobalV::GAMMA_ONLY_LOCAL)
        {
            this->gint_kernel_vlocal(na_grid,
                                        grid_index,
                                        delta_r,
                                        vldr3,
                                        LD_pool,
                                        pvpR_thread,
                                        ucell,
                                        hRGint_thread);
        }
#else
        if (GlobalV::GAMMA_ONLY_LOCAL && lgd > 0)
        {
            this->gint_kernel_vlocal(na_grid, grid_index, delta_r, vldr3, LD_pool, ucell, nullptr);
        }
        if (!GlobalV::GAMMA_ONLY_LOCAL)
        {
            this->gint_kernel_vlocal(na_grid,
                                        grid_index,
                                        delta_r,
                                        vldr3,
                                        LD_pool,
                                        ucell,
                                        this->pvpR_reduced[inout->ispin]);
        }
#endif
    delete[] vldr3;
    }
    if (GlobalV::GAMMA_ONLY_LOCAL && lgd > 0)
    {
#ifdef _OPENMP
#pragma omp critical(gint_gamma)
#endif
        {
            BlasConnector::axpy(this->hRGint->get_nnr(),
                                1.0,
                                hRGint_thread->get_wrapper(),
                                1,
                                this->hRGint->get_wrapper(),
                                1);
        }
        delete hRGint_thread;
    }
    if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
#ifdef _OPENMP
#pragma omp critical(gint_k)
#endif
        {
            BlasConnector::axpy(nnrg, 1.0, pvpR_thread, 1, pvpR_reduced[inout->ispin], 1);
        }
        delete[] pvpR_thread;   
    }
    ModuleBase::TITLE("Gint_interface", "cal_gint_vlocal");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal");
}

void Gint::cpu_dvlocal_interface(Gint_inout* inout)
{
    ModuleBase::TITLE("Gint_interface", "cal_gint_dvlocal");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_dvlocal");
    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        ModuleBase::WARNING_QUIT("Gint_interface::cal_gint", "dvlocal only for k point!");
    }
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * ucell.nwmax;
    const int lgd = this->gridt->lgd;
    const int nnrg = this->gridt->nnrg;
    const int ncyz = this->ny * this->nplane; 
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;
    

    double* pvdpRx_thread = nullptr;
    double* pvdpRy_thread = nullptr;
    double* pvdpRz_thread = nullptr;
    hamilt::HContainer<double>* hRGint_thread = nullptr; 

    pvdpRx_thread = new double[nnrg];   
    ModuleBase::GlobalFunc::ZEROS(pvdpRx_thread, nnrg);
    pvdpRy_thread = new double[nnrg];
    ModuleBase::GlobalFunc::ZEROS(pvdpRy_thread, nnrg);
    pvdpRz_thread = new double[nnrg];
    ModuleBase::GlobalFunc::ZEROS(pvdpRz_thread, nnrg);
    hRGint_thread = new hamilt::HContainer<double>(*this->hRGint);

    ModuleBase::GlobalFunc::ZEROS(this->pvdpRx_reduced[inout->ispin], nnrg);
    ModuleBase::GlobalFunc::ZEROS(this->pvdpRy_reduced[inout->ispin], nnrg);
    ModuleBase::GlobalFunc::ZEROS(this->pvdpRz_reduced[inout->ispin], nnrg);
#ifdef _OPENMP
#pragma omp for
#endif
    for (int grid_index = 0; grid_index < this->nbxx; grid_index++)
    {
        const int na_grid = this->gridt->how_many_atoms[grid_index];
        if (na_grid == 0)
        {
            continue;
        }
        double* vldr3 = Gint_Tools::get_vldr3(inout->vl,
                                                this->bxyz,
                                                this->bx,
                                                this->by,
                                                this->bz,
                                                this->nplane,
                                                this->gridt->start_ind[grid_index],
                                                ncyz,
                                                dv);
#ifdef _OPENMP
        this->gint_kernel_dvlocal(na_grid,
                                    grid_index,
                                    delta_r,
                                    vldr3,
                                    LD_pool,
                                    pvdpRx_thread,
                                    pvdpRy_thread,
                                    pvdpRz_thread,
                                    ucell);
#else
        this->gint_kernel_dvlocal(na_grid,
                                    grid_index,
                                    delta_r,
                                    vldr3,
                                    LD_pool,
                                    this->pvdpRx_reduced[inout->ispin],
                                    this->pvdpRy_reduced[inout->ispin],
                                    this->pvdpRz_reduced[inout->ispin],
                                    ucell);
#endif
        delete[] vldr3;
    }
    delete[] pvdpRx_thread;
    delete[] pvdpRy_thread;
    delete[] pvdpRz_thread;
    delete hRGint_thread;
    ModuleBase::TITLE("Gint_interface", "cal_gint_dvlocal");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_dvlocal");
}

void Gint::cpu_vlocal_meta_interface(Gint_inout* inout)
{
    ModuleBase::TITLE("Gint_interface", "cal_gint_vlocal_meta");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_vlocal_meta");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * ucell.nwmax;
    const int lgd = this->gridt->lgd;
    const int nnrg = this->gridt->nnrg;
    const int ncyz = this->ny * this->nplane; 
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;
    double* pvpR_thread = nullptr;
    hamilt::HContainer<double>* hRGint_thread = nullptr; 
    if (!pvpR_alloc_flag)
    {
        ModuleBase::WARNING_QUIT("Gint_interface::cal_gint", "pvpR has not been allocated yet!");
    }
    else
    {
        ModuleBase::GlobalFunc::ZEROS(this->pvpR_reduced[inout->ispin], nnrg);
    }

    if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
        pvpR_thread = new double[nnrg]();
        ModuleBase::GlobalFunc::ZEROS(pvpR_thread, nnrg);
    }

    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        hRGint_thread = new hamilt::HContainer<double>(*this->hRGint);
    }

    for (int grid_index = 0; grid_index < this->nbxx; grid_index++)
    {
        const int na_grid = this->gridt->how_many_atoms[grid_index];
        if (na_grid == 0)
        {
            continue;
        }
        double* vldr3 = Gint_Tools::get_vldr3(inout->vl,
                                                this->bxyz,
                                                this->bx,
                                                this->by,
                                                this->bz,
                                                this->nplane,
                                                this->gridt->start_ind[grid_index],
                                                ncyz,
                                                dv);
        double* vkdr3 = Gint_Tools::get_vldr3(inout->vofk,
                                                this->bxyz,
                                                this->bx,
                                                this->by,
                                                this->bz,
                                                this->nplane,
                                                this->gridt->start_ind[grid_index],
                                                ncyz,
                                                dv);
#ifdef _OPENMP
        if ((GlobalV::GAMMA_ONLY_LOCAL && lgd > 0) || !GlobalV::GAMMA_ONLY_LOCAL)
        {
            this->gint_kernel_vlocal_meta(na_grid,
                                            grid_index,
                                            delta_r,
                                            vldr3,
                                            vkdr3,
                                            LD_pool,
                                            pvpR_thread,
                                            ucell,
                                            hRGint_thread);
        }
#else
        if (GlobalV::GAMMA_ONLY_LOCAL && lgd > 0)
        {
            this->gint_kernel_vlocal_meta(na_grid,
                                            grid_index,
                                            delta_r,
                                            vldr3,
                                            vkdr3,
                                            LD_pool,
                                            ucell,
                                            nullptr);
        }
        if (!GlobalV::GAMMA_ONLY_LOCAL)
        {
            this->gint_kernel_vlocal_meta(na_grid,
                                            grid_index,
                                            delta_r,
                                            vldr3,
                                            vkdr3,
                                            LD_pool,
                                            ucell,
                                            this->pvpR_reduced[inout->ispin]);
        }
#endif
        delete[] vldr3;
        delete[] vkdr3;                                     
    }
if (GlobalV::GAMMA_ONLY_LOCAL)
    {
#ifdef _OPENMP
#pragma omp critical(gint_gamma)
#endif
        {
            BlasConnector::axpy(this->hRGint->get_nnr(),
                                1.0,
                                hRGint_thread->get_wrapper(),
                                1,
                                this->hRGint->get_wrapper(),
                                1);
        }
        delete hRGint_thread;
    }
    if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
#ifdef _OPENMP
#pragma omp critical(gint_k)
#endif
        {
            BlasConnector::axpy(nnrg, 1.0, pvpR_thread, 1, pvpR_reduced[inout->ispin], 1);
        }
        delete[] pvpR_thread;   
    }
}