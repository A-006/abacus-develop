#ifndef READ_STRU_H
#define READ_STRU_H

#include "atom_spec.h"
#include "module_cell/unitcell.h"
namespace unitcell
{
    bool check_tau(const Atom* atoms,
                   const int ntype,
                   const int lat0);
                   
    int read_atom_species(std::ifstream& ifa,
                          std::ofstream& ofs_running,
                          UnitCell& ucell); // read in the atom information for each type of atom

}
#endif // READ_STRU_H