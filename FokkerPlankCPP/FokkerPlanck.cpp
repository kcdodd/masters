#include "FokkerPlanck.hpp"
#include "Constants.hpp"
#include <cmath>
#include <iostream>

#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

using namespace std;

template<class T>
FokkerPlanck<T>::FokkerPlanck(
        const std::vector<T>& excitation_energies,
        const std::vector< dynarray<T> >& excitation_rates,
        const T ionization_energy,
        const dynarray<T> & ionization_rate,
        unsigned int num_velocities,
        T max_velocity,
        T Te,
        T ne,
        T rate_multiplier)
    :
        mExcitationEnergies(excitation_energies),
        mIonizationEnergy(ionization_energy),

        mF(num_velocities, 0, max_velocity),
        mTe(Te),
        mNe(ne),
        mRateMultiplier(rate_multiplier),

        m_tmp1_v_c(num_velocities, 0, max_velocity),
        m_tmp2_v_c(num_velocities, 0, max_velocity),
        m_tmp3_v_c(num_velocities, 0, max_velocity),
        m_tmp4_v_c(num_velocities, 0, max_velocity),
        m_tmp1_v_d(num_velocities, 0, max_velocity),
        m_tmp2_v_d(num_velocities, 0, max_velocity),
        m_tmp3_v_d(num_velocities, 0, max_velocity),
        m_tmp4_v_d(num_velocities, 0, max_velocity)


    {
        T energy_factor = 0.5*ELECTRON_MASS/ELEMENTAL_CHARGE;

        mExcitationRates.resize(excitation_rates.size());

        // put excitation rates into velocity space
        for (unsigned int r = 0; r < mExcitationRates.size(); ++r)
        {
            dynarray<T> & rates_v = mExcitationRates[r];
            const dynarray<T> & rates_e = excitation_rates[r];

            rates_v.resize(mF.length());
            rates_v.set_xminmax(mF.xmin(), mF.xmax());


            for(unsigned int i = 0; i < rates_v.length(); ++i)
            {
                T energy = energy_factor*rates_v.x(i)*rates_v.x(i);
                rates_v(i) = rates_e.interp(energy);
            }

        }


        mIonizationRate.resize(mF.length());
        mIonizationRate.set_xminmax(mF.xmin(), mF.xmax());

        for(unsigned int i = 0; i < mIonizationRate.length(); ++i)
        {
            T energy = energy_factor*mIonizationRate.x(i)*mIonizationRate.x(i);
            mIonizationRate(i) = ionization_rate.interp(energy);
        }

        // inital F is maxwellian
        for (unsigned int i = 0; i < mF.length(); ++i)
        {
            mF(i) = exp(-0.5*ELECTRON_MASS*pow(mF.x(i),2.0)/(ELEMENTAL_CHARGE*mTe));
        }

        dynarray<T> F_tmp(mF.length(), mF.xmin(), mF.xmax());

        // normalize

        for(unsigned int i = 0; i < mF.length(); ++i)
        {
            F_tmp(i) = 4.0*PI_CONST*mF(i)*pow(mF.x(i),2.0);

        }

        T A = F_tmp.integral_total();


        for(unsigned int i = 0; i < mF.length(); ++i)
        {
            mF(i) /= A;
        }
    }

template<class T>
dynarray<T>& FokkerPlanck<T>::f()
    {
        return mF;
    }

template<class T>
T FokkerPlanck<T>::step(
        T dt,
        T dt_rate_recalc,
        T min_equib_rate,
        unsigned int max_iterations,
        T& Te_eff)
    {

        T D0 = 0; // heating diffusion term
        T dfdt = 0;
        T equib_rate = 0;
        unsigned int iterations = 0;

        unsigned int num_steps_between_recalc = ceil(dt_rate_recalc/dt);


        dynarray<T> F_derivative(mF.length(), mF.xmin(), mF.xmax());

        dynarray<T> F_2nd_derivative(mF.length(), mF.xmin(), mF.xmax());
        dynarray<T> D(mF.length(), mF.xmin(), mF.xmax());
        dynarray<T> R(mF.length(), mF.xmin(), mF.xmax());
        dynarray<T> K_exc(mF.length(), mF.xmin(), mF.xmax());
        dynarray<T> K_ion(mF.length(), mF.xmin(), mF.xmax());
        dynarray<T> D_derivative(mF.length(), mF.xmin(), mF.xmax());
        dynarray<T> R_derivative(mF.length(), mF.xmin(), mF.xmax());


        dynarray<T> F_tmp(mF.length(), mF.xmin(), mF.xmax());

        dynarray<T> tmp1_v_a(mF.length(), mF.xmin(), mF.xmax());

        dynarray<T> tmp1_v_b(mF.length(), mF.xmin(), mF.xmax());


        dynarray<T> M1_tmp(mF.length(), mF.xmin(), mF.xmax());
        dynarray<T> M2_tmp(mF.length(), mF.xmin(), mF.xmax());

        dynarray<T> K_power_tmp(mF.length(), mF.xmin(), mF.xmax());
        dynarray<T> D0_power_tmp(mF.length(), mF.xmin(), mF.xmax());

        // adjustable to help match target Te
        T heating_factor = 1;


        do{


            mF.derivative(F_derivative);
            mF.derivative_2(F_2nd_derivative);

            if (iterations%num_steps_between_recalc == 0 || iterations == 0)
            {

                calc_D(mF, D, tmp1_v_a);
                calc_R(mF, F_derivative, R, tmp1_v_b);
                calc_K_exc(K_exc);
                calc_K_ion(K_ion);

/*
                boost::thread_group rate_threads;

                rate_threads.create_thread(boost::bind(&FokkerPlanck<T>::calc_D, this, boost::ref(mF), boost::ref(D), boost::ref(tmp1_v_a)));
                rate_threads.create_thread(boost::bind(&FokkerPlanck<T>::calc_R, this, boost::ref(mF), boost::ref(F_derivative), boost::ref(R), boost::ref(tmp1_v_b)));
                rate_threads.create_thread(boost::bind(&FokkerPlanck<T>::calc_K_exc, this, boost::ref(K_exc)));
                rate_threads.create_thread(boost::bind(&FokkerPlanck<T>::calc_K_ion, this, boost::ref(K_ion)));

                rate_threads.join_all();
*/
                D.derivative(D_derivative);
                R.derivative(R_derivative);

            }


            for(unsigned int i = 1; i < mF.length(); ++i)
            {

                M1_tmp(i) = 2.0*F_derivative(i)/mF.x(i);
                M2_tmp(i) = 2.0*R(i)/mF.x(i);

            }

            M1_tmp(0) = M1_tmp(1); // NaN b/c of 1/v, but does converge so use nearest point
            M2_tmp(0) = M2_tmp(1);

            // determin D0 to balance power loss from K
            for (unsigned int i = 0; i < mF.length(); ++i)
            {

                K_power_tmp(i) = (K_exc(i) + K_ion(i))*pow(mF.x(i),4.0);

                D0_power_tmp(i) = (M1_tmp(i) + F_2nd_derivative(i))*pow(mF.x(i),4.0);
            }

            D0 = -K_power_tmp.integral_total()/D0_power_tmp.integral_total();


            //cerr << K_power_tmp.integral_total() << endl;

            // calculate effective temperature
            for(unsigned int i = 0; i < mF.length(); ++i)
            {
                tmp1_v_a(i) = 4.0*PI_CONST*mF(i)*pow(mF.x(i),4.0);

            }

            Te_eff = (1.0/3.0)*ELECTRON_MASS*tmp1_v_a.integral_total()/ELEMENTAL_CHARGE;

            //cerr << heating_factor << endl;

            heating_factor = 1 + (mTe-Te_eff)/mTe;

            //cerr << D0 << endl;

            // redo the M1

            for(unsigned int i = 1; i < mF.length(); ++i)
            {

                M1_tmp(i) = 2.0*F_derivative(i)*(D(i) + heating_factor*D0)/mF.x(i);

            }

            M1_tmp(0) = M1_tmp(1);


            // calc next time step
            for(unsigned int i = 0; i < mF.length(); ++i)
            {

                dfdt = M1_tmp(i) + (D_derivative(i) - R(i))*F_derivative(i) + F_2nd_derivative(i)*(D(i) + heating_factor*D0) - (M2_tmp(i) + R_derivative(i))*mF(i) + K_exc(i) + K_ion(i);

                //std::cerr << dfdt << std::endl;

                F_tmp(i) = max(0.0, mF(i) + dt*dfdt);
            }

            // re-normalize

            for(unsigned int i = 0; i < mF.length(); ++i)
            {
                tmp1_v_a(i) = 4.0*PI_CONST*F_tmp(i)*pow(F_tmp.x(i),2.0);

            }

            T A = tmp1_v_a.integral_total();


            for(unsigned int i = 0; i < mF.length(); ++i)
            {
                F_tmp(i) = F_tmp(i)/A;
            }

            // calculate equib rate
            for(unsigned int i = 0; i < mF.length(); ++i)
            {
                tmp1_v_a(i) = 4.0*PI_CONST*abs(F_tmp(i)-mF(i))*pow(mF.x(i),2.0);

                // can now copy over new value of F.
                mF(i) = F_tmp(i);

            }

            T delta = tmp1_v_a.integral_total();

            equib_rate = delta/dt;


            ++iterations;

            if (iterations%10000 == 0)
            {

                //std::cerr << "Iteration " << iterations << ", Te_eff = " << Te_eff << ", equib_rate = " << equib_rate << " (" << delta << ")" << std::endl;
            }

        }while(equib_rate > min_equib_rate && iterations < max_iterations);

        return equib_rate;


    }

template<class T>
void FokkerPlanck<T>::calc_D(
    const dynarray<T> & F,
    dynarray<T> & out,
    dynarray<T> & tmp)
    {
        if (out.length() != F.length())
        {
            out.resize(F.length());
        }

        out.set_xminmax(F.xmin(), F.xmax());

        T factor = collision_factor();


        for (unsigned int i = 0; i < F.length(); ++i)
        {

            for (unsigned int j = 0; j < i; ++j)
            {
                // velocities below current velocity

                tmp(j) = F(j)*F.x(j)*(4.0/3.0)*pow(F.x(j)/F.x(i), 3.0);

            }

            for (unsigned int j = i; j < F.length(); ++j)
            {
                // velocities above current velocity

                tmp(j) = F(j)*F.x(j)*(4.0/3.0);

            }

            out(i) = factor*2.0*PI_CONST*tmp.integral_total();

        }
    }

template<class T>
void FokkerPlanck<T>::calc_R(
    const dynarray<T> & F,
    const dynarray<T> & F_derivative,
    dynarray<T> & out,
    dynarray<T> & tmp)
    {
        if (out.length() != F.length())
        {
            out.resize(F.length());
        }

        out.set_xminmax(F.xmin(), F.xmax());

        T factor = collision_factor();

        //cerr << factor << endl;

        for (unsigned int i = 0; i < F.length(); ++i)
        {

            for (unsigned int j = 0; j < i; ++j)
            {
                // velocities below current velocity

                tmp(j) = F_derivative(j)*(4.0/3.0)*pow(F.x(j)/F.x(i), 3.0);

            }

            for (unsigned int j = i; j < F.length(); ++j)
            {
                // velocities above current velocity

                tmp(j) = F_derivative(j)*(4.0/3.0);
            }

            out(i) = factor*2.0*PI_CONST*F.x(i)*tmp.integral_total();

        }
    }

template<class T>
void FokkerPlanck<T>::calc_K_exc(dynarray<T> & out)
    {


        dynarray<T>& loss_v = m_tmp1_v_c;
        dynarray<T>& gain_v = m_tmp2_v_c;

        T energy_factor = 0.5*ELECTRON_MASS/ELEMENTAL_CHARGE;


        for (unsigned int i = 0; i < mF.length(); ++i)
        {
            loss_v(i) = 0;
            gain_v(i) = 0;

            T v2 = pow(mF.x(i), 2.0);

            for (unsigned int r = 0; r < mExcitationRates.size(); ++r)
            {
                const dynarray<T> & rate_function_r = mExcitationRates[r];


                if (energy_factor*v2 > mExcitationEnergies[r])
                {

                    loss_v(i) += mF(i)*rate_function_r(i);

                }

                T s_v2 = v2 + mExcitationEnergies[r]/energy_factor;
                T source_velocity = sqrt(s_v2);

                if (source_velocity <= mF.xmax())
                {
                    // ????
                    gain_v(i) += (source_velocity/mF.x(i))*mF.interp(source_velocity)*rate_function_r.interp(source_velocity);
                }


            }

        }



        gain_v(0) = gain_v(1);


        // make sure total rate will be exactly zero

        dynarray<T>& gain_integral = m_tmp3_v_c;
        dynarray<T>& loss_integral = m_tmp4_v_c;

        for(unsigned int i = 0; i < mF.length(); ++i)
        {
            gain_integral(i) = gain_v(i)*4*PI_CONST*mF.x(i)*mF.x(i);
            loss_integral(i) = loss_v(i)*4*PI_CONST*mF.x(i)*mF.x(i);
        }

        T G_O_v = gain_integral.integral_total();
        T L_O_v = loss_integral.integral_total();



        T gain_norm_factor = L_O_v/G_O_v;

        //cerr << gain_norm_factor << endl;

        for(unsigned int i = 0; i < mF.length(); ++i)
        {
            out(i) = mRateMultiplier*(gain_norm_factor*gain_v(i) - loss_v(i));

        }


    }


template<class T>
void FokkerPlanck<T>::calc_K_ion(dynarray<T> & out)
    {

        dynarray<T>& loss_v = m_tmp1_v_d;
        dynarray<T>& gain_v = m_tmp2_v_d;
        dynarray<T>& gain_density = m_tmp3_v_d;
        dynarray<T>& gain_integral = m_tmp4_v_d;

        T energy_factor = 0.5*ELECTRON_MASS/ELEMENTAL_CHARGE;

        // construct gain density function
        for(unsigned int i = 0; i < mIonizationRate.length(); ++i)
        {
            T v2 = pow(mF.x(i),2.0);
            T energy = energy_factor*v2;

            if (energy > mIonizationEnergy)
            {

                gain_density(i) = v2*mF(i)*mIonizationRate(i)/(energy - mIonizationEnergy);
            }else{
                gain_density(i) = 0;
            }

            //cerr << gain_density(i);
        }

        gain_density.integral(gain_integral);

        // use limits of integral to calculate gain
        for(unsigned int i = 1; i < mF.length(); ++i)
        {
            T v_min = sqrt(pow(mF.x(i),2.0) + mIonizationEnergy/energy_factor);

            unsigned int i_min = min(mF.length()-1, (unsigned int)ceil((v_min - mF.xmin())/mF.dx()));

            // factor of 2 comes from the assumption that secondary electrons are distributed the same as final energy of incident electrons
            gain_v(i) = 2*(gain_integral(mF.length()-1) - gain_integral(i_min))*ELECTRON_MASS/(ELEMENTAL_CHARGE*mF.x(i));
        }

        gain_v(0) = gain_v(1);

        // calculate total loss rate
        for(unsigned int i = 0; i < mF.length(); ++i)
        {

            loss_v(i) = mF(i)*mIonizationRate(i)*4*PI_CONST*mF.x(i)*mF.x(i);
        }

        // calulate total ionization rate
        T L_0 = loss_v.integral_total();

        // total loss rate needs to be twice the ionization rate, to get equilibrium solution

        for(unsigned int i = 0; i < mF.length(); ++i)
        {

            loss_v(i) = mF(i)*(L_0 + mIonizationRate(i));
        }


        // make sure total rate will be exactly zero

        dynarray<T>& loss_integral = m_tmp3_v_d;

        for(unsigned int i = 0; i < mF.length(); ++i)
        {
            gain_integral(i) = gain_v(i)*4*PI_CONST*mF.x(i)*mF.x(i);
            loss_integral(i) = loss_v(i)*4*PI_CONST*mF.x(i)*mF.x(i);
        }

        T G_O_v = gain_integral.integral_total();
        T L_O_v = loss_integral.integral_total();


        if (L_O_v > 0)
        {
            T gain_norm_factor = L_O_v/G_O_v;

            //cerr << gain_norm_factor << endl;


            for(unsigned int i = 0; i < mF.length(); ++i)
            {
                out(i) = mRateMultiplier*(gain_norm_factor*gain_v(i) - loss_v(i));

            }
        }else{

            for(unsigned int i = 0; i < mF.length(); ++i)
            {
                out(i) = 0;

            }
        }


    }

template<class T>
T FokkerPlanck<T>::collision_factor()
    {
        // from NRL plasma formulary for e-e collisions
        T coulomb_logarithm = 23.5 - 0.5*log(mNe/1E6) + (5.0/4.0)*log(mTe) - sqrt(1E-5 + pow(log(mTe)-2.0, 2)/16.0);

        //cerr << coulomb_logarithm << endl;


        return mNe*pow(ELEMENTAL_CHARGE,4.0)*coulomb_logarithm/(8.0*PI_CONST*pow(ELECTRON_MASS*ELECTRIC_CONST, 2.0));
    }
