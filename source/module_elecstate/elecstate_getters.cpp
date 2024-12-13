#include "module_elecstate/elecstate_getters.h"

#include "module_cell/unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_parameter/parameter.h"
#include "module_hamilt_general/module_xc/xc_functional.h"

namespace elecstate
{

double get_ucell_omega()
{
    return GlobalC::ucell.omega;
}

double get_ucell_tpiba()
{
    return GlobalC::ucell.tpiba;
}

int get_xc_func_type()
{
    return XC_Functional::get_func_type();
}

std::string get_input_vdw_method()
{
    return PARAM.inp.vdw_method;
}

std::string get_ks_solver_type()
{
    return PARAM.inp.ks_solver;
}

} // namespace elecstate
