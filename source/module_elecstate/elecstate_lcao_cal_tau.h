#ifndef ELECSTATE_LCAO_CAL_TAU_H
#define ELECSTATE_LCAO_CAL_TAU_H

namespace elecstate
{

    void lcao_cal_tau_k(Gint_k* gint_k, 
                        Charge* charge);

    void lcao_cal_tau_gamma(Gint_Gamma* gint_gamma,
                            Charge* charge);

    template <typename T>
    void lcao_cal_tau(Gint_Gamma* gint_gamma, 
                      Gint_k* gint_k, 
                      Charge* charge);

}
#endif