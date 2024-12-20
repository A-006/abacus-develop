#include "bcast_cell.h"
#include "mpi.h"
void bcast_atoms_tau(Atom* atoms) {
    int ntype = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < ntype; i++) {
        atoms[i].bcast_atom(); // bcast tau array
    }
}