#ifndef FOKKERPLANCK_HPP_INCLUDED
#define FOKKERPLANCK_HPP_INCLUDED

#include "dynarray.hpp"
#include "dynmatrix.hpp"

#include <vector>

template<class T>
class FokkerPlanck
{
private:

    std::vector<T> mExcitationEnergies;
    std::vector< dynarray<T> > mExcitationRates;
    T mIonizationEnergy;
    dynarray<T> mIonizationRate;


    dynarray<T> mF;
    T mTe;
    T mNe;
    T mRateMultiplier;

    dynarray<T> m_tmp1_v_c;
    dynarray<T> m_tmp2_v_c;
    dynarray<T> m_tmp3_v_c;
    dynarray<T> m_tmp4_v_c;

    dynarray<T> m_tmp1_v_d;
    dynarray<T> m_tmp2_v_d;
    dynarray<T> m_tmp3_v_d;
    dynarray<T> m_tmp4_v_d;



public:

    FokkerPlanck(
        const std::vector<T>& excitation_energies,
        const std::vector< dynarray<T> >& excitation_rates,
        const T ionization_energy,
        const dynarray<T> & ionization_rate,
        unsigned int num_velocities,
        T max_velocity,
        T Te,
        T ne,
        T rate_multiplier);

    dynarray<T>& f();

    T step(
        T dt,
        T dt_rate_recalc,
        T min_equib_rate,
        unsigned int max_iterations,
        T& Te_eff);

    void calc_D(const dynarray<T> & F, dynarray<T> & out, dynarray<T> & tmp);
    void calc_R(const dynarray<T> & F, const dynarray<T> & F_derivative, dynarray<T> & out, dynarray<T> & tmp);

    void calc_K_exc(dynarray<T> & out);

    void calc_K_ion(dynarray<T> & out);

    T collision_factor();

};

#endif // FOKKERPLANCK_HPP_INCLUDED
