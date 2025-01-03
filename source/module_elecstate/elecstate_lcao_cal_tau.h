#ifndef ELECSTATE_LCAO_CAL_TAU_H
#define ELECSTATE_LCAO_CAL_TAU_H

namespace elecstate
{

    void cal_tau_k(Gint_k& gint_k, 
                   Charge& charge);

    void cal_tau_gamma(Gint_Gamma& gint_gamma,
                       Charge& charge);

}
#endif