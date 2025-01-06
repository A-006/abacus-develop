#ifndef READ_STRU_H
#define READ_STRU_H

#include "atom_spec.h"
namespace unitcell
{
    bool check_tau(const Atom* atoms,
                   const int ntype,
                   const int lat0);
}
#endif // READ_STRU_H