#ifndef ELECSTATE_GETTERS_H
#define ELECSTATE_GETTERS_H

#include <string>

// Description: Getters for elecstate module
namespace elecstate
{

/// @brief get the value of GlobalC::ucell.omega
double get_ucell_omega();
/// @brief get the value of GlobalC::ucell.tpiba
double get_ucell_tpiba();
/// @brief get the value of XC_Functional::func_type
int get_xc_func_type();
/// @brief get the type of KS_SOLVER
std::string get_ks_solver_type();

} // namespace elecstate

#endif // ELECSTATE_GETTERS_H
