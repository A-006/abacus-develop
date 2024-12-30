#ifndef UPDATE_CELL_H
#define UPDATE_CELL_H

#include "unitcell_data.h"
#include "unitcell.h"  

/*
this file is used to update the cell,contains the following functions:
1. remake_cell: for constrained vc-relaxation where type of lattice 
is fixed, adjust the lattice vectors
2. setup_cell_after_vc: setup cell after vc-relaxation
the functions are defined in the namespace UnitCell,
Accually, the functions are focused on the cell-relax part functions
of the UnitCell class.
3. periodic_boundary_adjustment: adjust the boundary of the cell
*/
namespace unitcell
{
    // for constrained vc-relaxation where type of lattice
    // is fixed, adjust the lattice vectors
    void remake_cell(Lattice& lat);

    void setup_cell_after_vc(UnitCell& ucell, std::ofstream& log);

    void periodic_boundary_adjustment(Atom* atoms,
                                      const ModuleBase::Matrix3& latvec,
                                      const int ntype);

    /** 
    * update the position and tau of the atoms
    * 
    * @param lat: the lattice of the atoms [in]
    * @param pos: the position of the atoms [in]
    * @param ntype: the number of types of the atoms [in]
    * @param nat: the number of atoms [in]
    * @param atoms: the atoms to be updated [out]
    */
    void update_pos_tau(const Lattice& lat,
                        const double* pos,
                        const int ntype,
                        const int nat,
                        Atom* atoms);
}
//
#endif // UPDATE_CELL_H