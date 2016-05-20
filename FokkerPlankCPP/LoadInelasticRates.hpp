#ifndef LOADINELASTICRATES_HPP_INCLUDED
#define LOADINELASTICRATES_HPP_INCLUDED

#include "dynarray.hpp"
#include <vector>

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

template<class T>
class LoadInelasticRates
{
private:

    string mFilename;

public:

    LoadInelasticRates(const string& filename)
    :
        mFilename(filename)
    {

    }

    void load(
        std::vector<T>& excitation_energies,
        std::vector< dynarray<T> >& excitation_rates,
        T& ionization_energy,
        dynarray<T> & ionization_rate)
    {
        ifstream inputfile(mFilename.c_str());

        unsigned int num_energies;
        unsigned int num_rates;


        if (inputfile.is_open())
        {
            double tmp;
            inputfile >> tmp;
            num_energies = (unsigned int) tmp;
            inputfile >> tmp;
            num_rates = (unsigned int) tmp;

            T emin;
            T emax;


            excitation_energies.resize(num_rates-1);
            excitation_rates.resize(num_rates-1);

            inputfile >> emin;

            for(unsigned int i = 1; i < num_energies-1; ++i)
            {

                inputfile >> tmp;

            }

            inputfile >> emax;


            for(unsigned int i = 0; i < num_rates-1; ++i)
            {

                inputfile >> excitation_energies[i];

            }

            // last one is ionization energy
            inputfile >> ionization_energy;



            for( unsigned int r = 0; r < num_rates-1; ++r)
            {
                dynarray<T> &rates_i = excitation_rates[r];
                rates_i.resize(num_energies);
                rates_i.set_xminmax(emin, emax);

                for(unsigned int i = 0; i < num_energies; ++i)
                {
                    inputfile >> rates_i(i);
                }
            }

            ionization_rate.resize(num_energies);
            ionization_rate.set_xminmax(emin, emax);

            for(unsigned int i = 0; i < num_energies; ++i)
            {
                inputfile >> ionization_rate(i);
            }

            inputfile.close();
        }

    }

};

#endif // LOADINELASTICRATES_HPP_INCLUDED
