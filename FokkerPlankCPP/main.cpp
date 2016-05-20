#include <iostream>

#include "dynarray.hpp"
#include "dynmatrix.hpp"

#include "LoadInelasticRates.hpp"

#include "FokkerPlanck.hpp"
#include "FokkerPlanck.cpp"

#include "SaveF.hpp"

#include <cmath>
#include "Constants.hpp"

#include <sstream>

#include <stdio.h>

#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

#include <time.h>

using namespace std;

template<class T>
class FP_Case{
private:

        const std::vector<T>& m_excitation_energies;
        const std::vector< dynarray<T> >& m_excitation_rates;
        T m_ionization_energy;
        const dynarray<T> & m_ionization_rate;
        unsigned int m_num_velocities;

public:

    FP_Case(
        const std::vector<T>& excitation_energies,
        const std::vector< dynarray<T> >& excitation_rates,
        T ionization_energy,
        const dynarray<T> & ionization_rate,
        unsigned int num_velocities)
        :
        m_excitation_energies(excitation_energies),
        m_excitation_rates(excitation_rates),
        m_ionization_energy(ionization_energy),
        m_ionization_rate(ionization_rate),
        m_num_velocities(num_velocities)

        {}



    void run(
        T Te,
        T ne,
        T rate_multiplier,
        T dt,
        string outfile,
        T vth_mult)
        {


                double Te_eff;
                double vth = sqrt(2*ELEMENTAL_CHARGE*Te/ELECTRON_MASS);

                FokkerPlanck<double> solver(m_excitation_energies, m_excitation_rates, m_ionization_energy, m_ionization_rate, m_num_velocities, vth_mult*vth, Te, ne, rate_multiplier);

                double equib_rate = solver.step(dt, 1000*dt, 50.0, 100E6, Te_eff);


                SaveF<double> saver;

                saver.save(outfile, solver.f(), Te, ne, rate_multiplier*1E17, equib_rate, Te_eff);



            cout << "Completed Te_ff = " << Te_eff << ", equib_rate = " << equib_rate << ": " << outfile << endl;
        }
};

int main()
{

    LoadInelasticRates<double> loader("ar_inelastic_rates_1E17_640eV.txt");
    std::vector<double> excitation_energies;
    std::vector< dynarray<double> > excitation_rates;
    double ionization_energy;
    dynarray<double> ionization_rate;


    loader.load(excitation_energies, excitation_rates, ionization_energy, ionization_rate);

    //cout << excitation_rates[0].dx() << endl;

    unsigned int num_Te = 6;
    double Te[] = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0};
    double dt[] = {5E-14, 1E-13, 1E-13, 1E-11, 1E-11, 1E-11};
    double vth_mult[] = {7, 5, 5, 5, 5, 5};


    unsigned int num_ne = 8;
    double ne[] = {1E16, 2E16, 5E16, 1E17, 2E17, 5E17, 1E18, 2E18};

    unsigned int num_n0 = 7;
    double n0[] = {1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0};


    FP_Case<double> cases(excitation_energies, excitation_rates, ionization_energy, ionization_rate, 501);

    cout << "Beginning Calculations" << endl;

    boost::thread_group case_threads;

    for (unsigned int Te_i = 0; Te_i < num_Te; ++Te_i)
    {
        for (unsigned int ne_i = 0; ne_i < num_ne; ++ne_i)
        {
            for (unsigned int n0_i = 0; n0_i < num_n0; ++n0_i)
            {
                stringstream outfile;

                outfile << "output7/ar_evdf_Te" << Te_i << "_ne" << ne_i << "_n0" << n0_i << ".txt";

                case_threads.create_thread(boost::bind(&FP_Case<double>::run, &cases, Te[Te_i], ne[ne_i], n0[n0_i], dt[Te_i], outfile.str(), vth_mult[Te_i]));

                //cases.run(Te[Te_i], ne[ne_i], n0[n0_i], dt[Te_i], outfile.str());

            } // n0

        } // ne

    } // Te

    case_threads.join_all();


    cout << "Fokker-Planck Done." << endl;

    return 0;
}
