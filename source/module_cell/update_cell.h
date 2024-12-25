#ifndef UPDATE_CELL_H
#define UPDATE_CELL_H

#include "unitcell_data.h"
#include "module_cell/unitcell.h"  
    // for constrained vc-relaxation where type of lattice
    // is fixed, adjust the lattice vectors
    void remake_cell(Lattice& lat);

    void setup_cell_after_vc(UnitCell& ucell, std::ofstream& log);
#endif