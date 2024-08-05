//=========================================================
// REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#ifndef GINT_TOOLS_H
#define GINT_TOOLS_H
#include "../grid_technique.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_base/array_pool.h"

#include <cstdlib>

namespace Gint_Tools
{
enum class job_type
{
    vlocal,
    rho,
    force,
    tau,
    vlocal_meta,
    force_meta,
    dvlocal
};
// Hamiltonian, electron density, force, kinetic energy density, Hamiltonian for mGGA

// if exponent is an integer between 0 and 5 (the most common cases in gint),
// pow_int is much faster than std::pow
inline double pow_int(const double base, const int exp)
{
    switch (exp)
    {
    case 0:
        return 1.0;
    case 1:
        return base;
    case 2:
        return base * base;
    case 3:
        return base * base * base;
    case 4:
        return base * base * base * base;
    case 5:
        return base * base * base * base * base;
    default:
        double result = std::pow(base, exp);
        return result;
    }
}
// vindex[pw.bxyz]

/**
 * @brief Get the vindex form the grid index
 * @param bxyz number of big grids
 * @param bx number of big grids in x direction
 * @param by number of big grids in y direction
 * @param bz number of big grids in z direction
 * @param nplane Currently using Z-axis 1D division, 
 * recording the number of the Z-axis process
 * (nbz in the current process).
 * @param start_ind start index of the grid in the 1D FFT grid
 * @param ncyz number of grids in yz plane
 * @param vindex the index of the grid 
*/
void get_vindex(const int bxyz, const int bx, const int by,
                    const int bz, const int nplane, 
                    const int start_ind,const int ncyz,int* vindex);

/**
 * @brief Get the vldr3 form the grid index
 * @param vldr3 the local potential multiplied by the grid volume
 * @param vlocal the local potential
 * @param bxyz number of grids
 * @param bx number of grids in x direction
 * @param by number of grids in y direction
 * @param bz number of grids in z direction
 * @param nplane Currently using Z-axis 1D division, 
 * recording the number of the Z-axis process
 * (nbz in the current process).
 * @param start_ind start index of the grid in the 1D FFT grid
 * @param ncyz number of grids in yz plane
 * @param dv the volume of the grid
*/
void get_gint_vldr3(double* vldr3,
                    const double* const vlocal,
                    const int bxyz,
                    const int bx,
                    const int by,
                    const int bz,
                    const int nplane,
                    const int start_ind,
                    const int ncyz,
                    const double dv);

/**
 * @brief Get the information of a big grid index
 * @param gt the grid technique, which contains the tools of the grid intergration
 * @param bxyz number of grids
 * @param na_grid number of atoms on this grid
 * @param grid_index 1d index of FFT index (i,j,k)
 * @param block_iw track the atom orbitals in all atoms
 * @param block_index count total number of atomis orbitals
 * @param block_size count the number of atomis orbitals in each atom
 * @param cal_flag whether the atom-grid distance is larger than cutoff
*/                    
void get_block_info(const Grid_Technique& gt, 
                    const int bxyz, 
                    const int na_grid, 
                    const int grid_index,
                    int* block_iw, 
                    int* block_index, 
                    int* block_size, 
                    bool** cal_flag);

void init_orb(double& dr_uniform,
              std::vector<double>& rcuts,
              UnitCell& ucell,
              std::vector<std::vector<double>>& psi_u,
              std::vector<std::vector<double>>& dpsi_u,
              std::vector<std::vector<double>>& d2psi_u);
/**
 * @brief Get the psi and dpsi from the it type atom
 * @param gt Grid_Technique
 * @param it index of the it type atom
 * @param atom the atom type contianing the nw, iw2_new
 * @param it_psi_uniform psi of the it type atom
 * @param it_dpsi_uniform dpsi of the it type atom
*/
void get_psi_dpsi(const Grid_Technique& gt,
                  int it,
                  Atom* atom,
                  std::vector<const double*>& it_psi_uniform,
                  std::vector<const double*>& it_dpsi_uniform);
/**
 * @brief Obtain the distance between the grid points and the atoms,
 *  and the type of atoms.
 * @param gt Grid_Technique.
 * @param bcell_start start index of the big cell.
 * @param it index of the atom type.
 * @param mt the distance between the big grid and the atom.
*/
void get_grid_bigcell_distance(const Grid_Technique& gt,
                                const int bcell_start,
								int& it,
								double* mt);
/**
 * @brief Calculate the distance between the meshcell and the atoms.
 * 
 * @param distance the distance between the meshcell and the atoms.
 * @param dr The three-dimensional distance from the starting position
 *  of the small grid point to the atom.
 * @param mt the distance between the big grid and the atom.
 * @param meshcell_pos The distance between the starting positions 
 * of the small grid point and the large grid point.
*/
void cal_grid_atom_distance(double &distance,
                            double* dr,
                            const double* mt,
                            const double* meshcell_pos);

/**
 * @brief Calculate the spherical harmonic functions Ylm.
 * 
 * @param distance the distance between the meshcell and the atoms.
 * @param delta_r the interval of atom segmentation.
 * @param atom the atom type contianing the nw, iw2_new
 * @param ylma the spherical harmonic functions Ylm.
 * @param it_psi_uniform psi of the it type atom
 * @param it_dpsi_uniform dpsi of the it type atom
*/
void spl_intrp(const double distance,
							const double delta_r,
							Atom*& atom,
							std::vector<double>& ylma,
							std::vector<const double*>& it_psi_uniform,
							std::vector<const double*>& it_dpsi_uniform,
							double *p);
/**
 * @brief Calculate the gradient of the 
 * spherical harmonic functions Ylm.
 * 
 * @param distance the distance between the meshcell and the atoms.
 * @param dr The three-dimensional distance from the starting position
 * of the small grid point to the atom.
 * @param delta_r the interval of atom segmentation.
 * @param atom the atom type contianing the nw, iw2_new
 * @param rly the spherical harmonic functions Ylm.
 * @param grly the gradient of the spherical harmonic functions Ylm.
 * @param it_psi_uniform psi of the it type atom
 * @param it_dpsi_uniform dpsi of the it type atom
 * @param p_psi psi of between the meshcell and the it type atom
 * @param p_dpsi_x gradient psi_x of the grid point
 * @param p_dpsi_y gradient psi_y of the grid point
 * @param p_dpsi_z gradient psi_z of the grid point
*/
void dpsi_spl_intrp(const double distance,
                    const double* dr,
                    const double delta_r,
                    Atom*& atom,
                    double* rly,
                    double** grly,
                    std::vector<const double*>& it_psi_uniform,
                    std::vector<const double*>& it_dpsi_uniform,
                    double *p_psi,
                    double *p_dpsi_x,
                    double *p_dpsi_y,
                    double *p_dpsi_z);

/**
 * @brief Calculate the spherical harmonic functions Ylm and its derivatives.
 * 
 * @param distance the distance between the meshcell and the atoms.
 * @param dr The three-dimensional distance from the starting position
 * of the small grid point to the atom.
 * @param delta_r the interval of atom segmentation.
 * @param atom the atom type contianing the nw, iw2_new
 * @param rly the spherical harmonic functions Ylm.
 * @param grly the gradient of the spherical harmonic functions Ylm.
 * @param it_psi_uniform psi of the it type atom.
 * @param it_dpsi_uniform dpsi of the it type atom.
 * @param psi psi and three gradient psi of 
 * between the meshcell and the it type atom.
*/
void dpsi_spl_intrp(const double distance1,
								const double* dr1,
								const double delta_r,
								const int i,
								Atom*& atom,
                                double* rly,
								double** grly,
                                std::vector<const double*>& it_psi_uniform,
                                std::vector<const double*>& it_dpsi_uniform,
                                double ***dpsi);
// psir_ylm[pw.bxyz][LD_pool]
void cal_psir_ylm(const Grid_Technique& gt,
                  const int bxyz,
                  const int na_grid,            // number of atoms on this grid
                  const int grid_index,         // 1d index of FFT index (i,j,k)
                  const double delta_r,         // delta_r of the uniform FFT grid
                  const int* const block_index, // count total number of atomis orbitals
                  const int* const block_size,
                  const bool* const* const cal_flag,
                  double* const* const psir_ylm); // whether the atom-grid distance is larger than cutoff

// psir_ylm and dpsir_ylm, both[pw.bxyz][LD_pool]
void cal_dpsir_ylm(
    const Grid_Technique& gt,
    const int bxyz,
    const int na_grid,                 // number of atoms on this grid
    const int grid_index,              // 1d index of FFT index (i,j,k)
    const double delta_r,              // delta_r of the uniform FFT grid
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const int* const block_size,       // block_size[na_grid],	number of columns of a band
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    double* const* const psir_ylm,
    double* const* const dpsir_ylm_x,
    double* const* const dpsir_ylm_y,
    double* const* const dpsir_ylm_z);

// dpsir_ylm * (r-R), R is the atomic position
void cal_dpsirr_ylm(
    const Grid_Technique& gt, const int bxyz,
    const int na_grid,                 // number of atoms on this grid
    const int grid_index,              // 1d index of FFT index (i,j,k)
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const int* const block_size,       // block_size[na_grid],	number of columns of a band
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    double* const* const dpsir_ylm_x, double* const* const dpsir_ylm_y, double* const* const dpsir_ylm_z,
    double* const* const dpsir_ylm);

void cal_ddpsir_ylm(
    const Grid_Technique& gt,
    const int bxyz,
    const int na_grid,                 // number of atoms on this grid
    const int grid_index,              // 1d index of FFT index (i,j,k)
    const double delta_r,              // delta_r of the uniform FFT grid
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const int* const block_size,       // block_size[na_grid],	number of columns of a band
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    double* const* const ddpsir_ylm_xx,
    double* const* const ddpsir_ylm_xy,
    double* const* const ddpsir_ylm_xz,
    double* const* const ddpsir_ylm_yy,
    double* const* const ddpsir_ylm_yz,
    double* const* const ddpsir_ylm_zz);

// psir_ylm * vldr3
ModuleBase::Array_Pool<double> get_psir_vlbr3(
    const int bxyz,
    const int na_grid, // how many atoms on this (i,j,k) grid
    const int LD_pool,
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    const double* const vldr3,         // vldr3[bxyz]
    const double* const* const psir_ylm); // psir_ylm[bxyz][LD_pool]

// sum_nu rho_mu,nu psi_nu, for gamma point
void mult_psi_DM(
    const Grid_Technique& gt,
    const int bxyz,
    const int na_grid, // how many atoms on this (i,j,k) grid
    const int LD_pool,
    const int* const block_iw,         // block_iw[na_grid],	index of wave functions for each block
    const int* const block_size,       // block_size[na_grid],	number of columns of a band
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    const double* const* const psi,    // psir_vlbr3[bxyz][LD_pool]
    double** psi_DM,
    const double* const* const DM,
    const bool if_symm);

// sum_nu,R rho_mu,nu(R) psi_nu, for multi-k
void mult_psi_DMR(const Grid_Technique& gt,
                  const int bxyz,
                  const int LD_pool,
                  const int& grid_index,
                  const int& na_grid,
                  const int* const block_index,
                  const int* const block_size,
                  bool** cal_flag,
                  double** psi,
                  double** psi_DMR,
                  const hamilt::HContainer<double>* DM,
                  const bool if_symm);

// sum_nu rho_mu,nu psi_nu, for gamma point
void mult_psi_DM_new(
    const Grid_Technique& gt,
    const int bxyz,
    const int& grid_index,
    const int na_grid, // how many atoms on this (i,j,k) grid
    const int LD_pool,
    const int* const block_iw,         // block_iw[na_grid],	index of wave functions for each block
    const int* const block_size,       // block_size[na_grid],	number of columns of a band
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    const double* const* const psi,    // psir_vlbr3[bxyz][LD_pool]
    double** psi_DM,
    const hamilt::HContainer<double>* DM,
    const bool if_symm);

} // namespace Gint_Tools
#endif