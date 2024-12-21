#ifndef UPDATE_CELL_H
#define UPDATE_CELL_H

#include "unitcell_data.h"
    // for constrained vc-relaxation where type of lattice
    // is fixed, adjust the lattice vectors
    void remake_cell(Lattice& lat);

#endif