#include "gint.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

void Gint::gint_kernel_rho(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_rho");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_rho");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * ucell.nwmax;
    const int lgd = this->gridt->lgd;
    const int ncyz = this->ny * this->nplane;
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;

    int block_iw[max_size];
    ModuleBase::GlobalFunc::ZEROS(block_iw, max_size);
    int block_index[max_size+1];
    ModuleBase::GlobalFunc::ZEROS(block_index, max_size+1);
    int block_size[max_size];
    ModuleBase::GlobalFunc::ZEROS(block_size, max_size);
    int vindex[bxyz];
    ModuleBase::GlobalFunc::ZEROS(vindex, bxyz);
#ifdef _OPENMP
#pragma omp parallel private(block_iw, block_index, block_size,vindex)
{
#pragma omp for
#endif
    for (int grid_index = 0; grid_index < this->nbxx; grid_index++) {
        const int na_grid = this->gridt->how_many_atoms[grid_index];
        if (na_grid == 0) {
            continue;
        }
        Gint_Tools::get_vindex_rho(this->bxyz,
                                    this->bx,
                                    this->by,
                                    this->bz,
                                    this->nplane,
                                    this->gridt->start_ind[grid_index],
                                    ncyz,
                                    vindex);
         // prepare block information
        // int *block_iw, *block_index, *block_size;
        ModuleBase::Array_Pool<bool> cal_flag(this->bxyz,max_size);
        Gint_Tools::get_block_info_vlocal(*this->gridt,
                                this->bxyz,
                                na_grid,
                                grid_index,
                                block_iw,
                                block_index,
                                block_size,
                                cal_flag.get_ptr_2D());

    // evaluate psi on grids
        ModuleBase::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
        Gint_Tools::cal_psir_ylm(*this->gridt,
                                this->bxyz,
                                na_grid,
                                grid_index,
                                delta_r,
                                block_index,
                                block_size,
                                cal_flag.get_ptr_2D(),
                                psir_ylm.get_ptr_2D());

        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            ModuleBase::Array_Pool<double> psir_DM(this->bxyz, LD_pool);
            ModuleBase::GlobalFunc::ZEROS(psir_DM.get_ptr_1D(), this->bxyz * LD_pool);
            if (GlobalV::GAMMA_ONLY_LOCAL)
            {
                Gint_Tools::mult_psi_DM_new(*this->gridt,
                                            this->bxyz,
                                            grid_index,
                                            na_grid,
                                            LD_pool,
                                            block_iw,
                                            block_size,
                                            block_index,
                                            cal_flag.get_ptr_2D(),
                                            psir_ylm.get_ptr_2D(),
                                            psir_DM.get_ptr_2D(),
                                            this->DMRGint[is],
                                            inout->if_symm);
            }
            else
            {
                // calculating g_mu(r) = sum_nu rho_mu,nu psi_nu(r)
                Gint_Tools::mult_psi_DMR(*this->gridt,
                                        this->bxyz,
                                        LD_pool,
                                        grid_index,
                                        na_grid,
                                        block_index,
                                        block_size,
                                        cal_flag.get_ptr_2D(),
                                        psir_ylm.get_ptr_2D(),
                                        psir_DM.get_ptr_2D(),
                                        this->DMRGint[is],
                                        inout->if_symm);
            }

            // do sum_mu g_mu(r)psi_mu(r) to get electron density on grid
            this->cal_meshball_rho(na_grid, block_index, vindex, psir_ylm.get_ptr_2D(), psir_DM.get_ptr_2D(), inout->rho[is]);
        }
        }
#ifdef _OPENMP
}
#endif
    ModuleBase::TITLE("Gint_interface", "cal_gint_rho");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_rho");
}

void Gint::gint_kernel_tau(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint_tau");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_tau");
    const UnitCell& ucell = *this->ucell;
    const int max_size = this->gridt->max_atom;
    const int LD_pool = max_size * ucell.nwmax;
    const int lgd = this->gridt->lgd;
    const int ncyz = this->ny * this->nplane;
    const double dv = ucell.omega / this->ncxyz;
    const double delta_r = this->gridt->dr_uniform;

    int block_iw[max_size];
    ModuleBase::GlobalFunc::ZEROS(block_iw, max_size);
    int block_index[max_size+1];
    ModuleBase::GlobalFunc::ZEROS(block_index, max_size+1);
    int block_size[max_size];
    ModuleBase::GlobalFunc::ZEROS(block_size, max_size);
    int vindex[bxyz];
    ModuleBase::GlobalFunc::ZEROS(vindex, bxyz);
#ifdef _OPENMP
#pragma omp parallel private(block_iw, block_index, block_size,vindex)
{
#pragma omp for
#endif
    for (int grid_index = 0; grid_index < this->nbxx; grid_index++) {
        const int na_grid = this->gridt->how_many_atoms[grid_index];
        if (na_grid == 0) {
            continue;
        }
        // int* vindex = Gint_Tools::get_vindex(ncyz, ibx, jby, kbz);
        Gint_Tools::get_vindex_rho(this->bxyz,
                                this->bx,
                                this->by,
                                this->bz,
                                this->nplane,
                                this->gridt->start_ind[grid_index],
                                ncyz,
                                vindex);
        //prepare block information
        ModuleBase::Array_Pool<bool> cal_flag(this->bxyz,max_size);
        Gint_Tools::get_block_info_vlocal(*this->gridt, this->bxyz, na_grid, grid_index, 
                                            block_iw, block_index, block_size, cal_flag.get_ptr_2D());

    //evaluate psi and dpsi on grids
        ModuleBase::Array_Pool<double> psir_ylm(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_x(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_y(this->bxyz, LD_pool);
        ModuleBase::Array_Pool<double> dpsir_ylm_z(this->bxyz, LD_pool);

        Gint_Tools::cal_dpsir_ylm(*this->gridt, 
            this->bxyz, na_grid, grid_index, delta_r,
            block_index, block_size, 
            cal_flag.get_ptr_2D(),
            psir_ylm.get_ptr_2D(),
            dpsir_ylm_x.get_ptr_2D(),
            dpsir_ylm_y.get_ptr_2D(),
            dpsir_ylm_z.get_ptr_2D());

        for(int is=0; is<GlobalV::NSPIN; ++is)
        {
            ModuleBase::Array_Pool<double> dpsix_DM(this->bxyz, LD_pool);
            ModuleBase::Array_Pool<double> dpsiy_DM(this->bxyz, LD_pool);
            ModuleBase::Array_Pool<double> dpsiz_DM(this->bxyz, LD_pool);
            ModuleBase::GlobalFunc::ZEROS(dpsix_DM.get_ptr_1D(), this->bxyz*LD_pool);
            ModuleBase::GlobalFunc::ZEROS(dpsiy_DM.get_ptr_1D(), this->bxyz*LD_pool);
            ModuleBase::GlobalFunc::ZEROS(dpsiz_DM.get_ptr_1D(), this->bxyz*LD_pool);

            //calculating g_i,mu(r) = sum_nu rho_mu,nu d/dx_i psi_nu(r), x_i=x,y,z
            if(GlobalV::GAMMA_ONLY_LOCAL)
            {
                Gint_Tools::mult_psi_DM_new(
                    *this->gridt,this->bxyz, grid_index, na_grid, LD_pool,
                    block_iw, block_size,
                    block_index, cal_flag.get_ptr_2D(),
                    dpsir_ylm_x.get_ptr_2D(),
                    dpsix_DM.get_ptr_2D(),
                    this->DMRGint[is], 1);
                Gint_Tools::mult_psi_DM_new(
                    *this->gridt, this->bxyz, grid_index, na_grid, LD_pool,
                    block_iw, block_size,
                    block_index, cal_flag.get_ptr_2D(),
                    dpsir_ylm_y.get_ptr_2D(),
                    dpsiy_DM.get_ptr_2D(),
                    this->DMRGint[is], 1);	
                Gint_Tools::mult_psi_DM_new(
                    *this->gridt, this->bxyz, grid_index, na_grid, LD_pool,
                    block_iw, block_size,
                    block_index, cal_flag.get_ptr_2D(),
                    dpsir_ylm_z.get_ptr_2D(),
                    dpsiz_DM.get_ptr_2D(),
                    this->DMRGint[is], 1);
            }
            else
            {
                Gint_Tools::mult_psi_DMR(
                    *this->gridt, this->bxyz,
                    LD_pool,
                    grid_index, na_grid,
                    block_index, block_size,
                    cal_flag.get_ptr_2D(), 
                    dpsir_ylm_x.get_ptr_2D(),
                    dpsix_DM.get_ptr_2D(),
                    this->DMRGint[is],
                    1);
                Gint_Tools::mult_psi_DMR(
                    *this->gridt, this->bxyz,
                    LD_pool,
                    grid_index, na_grid,
                    block_index, block_size,
                    cal_flag.get_ptr_2D(),
                    dpsir_ylm_y.get_ptr_2D(),
                    dpsiy_DM.get_ptr_2D(),
                    this->DMRGint[is],
                    1);
                Gint_Tools::mult_psi_DMR(
                    *this->gridt, this->bxyz,
                    LD_pool,
                    grid_index, na_grid,
                    block_index, block_size,
                    cal_flag.get_ptr_2D(), 
                    dpsir_ylm_z.get_ptr_2D(),
                    dpsiz_DM.get_ptr_2D(),
                    this->DMRGint[is],
                    1);
            }

        //do sum_i,mu g_i,mu(r) * d/dx_i psi_mu(r) to get kinetic energy density on grid
            if(inout->job==Gint_Tools::job_type::tau)
            {
                this->cal_meshball_tau(
                    na_grid, block_index,
                    vindex,
                    dpsir_ylm_x.get_ptr_2D(), dpsir_ylm_y.get_ptr_2D(), dpsir_ylm_z.get_ptr_2D(),
                    dpsix_DM.get_ptr_2D(), dpsiy_DM.get_ptr_2D(), dpsiz_DM.get_ptr_2D(),
                    inout->rho[is]);
            }
        }
    }
#ifdef _OPENMP
}
#endif
    ModuleBase::TITLE("Gint_interface", "cal_gint_tau");
    ModuleBase::timer::tick("Gint_interface", "cal_gint_tau");
}
