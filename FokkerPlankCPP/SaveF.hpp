#ifndef SAVEF_HPP_INCLUDED
#define SAVEF_HPP_INCLUDED

#include "dynarray.hpp"
#include <vector>

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

template<class T>
class SaveF
{

    public:

    SaveF() {}
    ~SaveF(){}

    void save(const string& filename, const dynarray<T> & F, T Te, T ne, T n0, T equib_rate, T Te_eff)
    {
        ofstream outfile(filename.c_str());

        if (outfile.is_open())
        {

            outfile << Te << endl;
            outfile << ne << endl;
            outfile << n0 << endl;
            outfile << equib_rate << endl;
            outfile << Te_eff << endl;

            // save dv
            outfile << F.dx() << endl;

            // save dist. function
            for (unsigned int i = 0; i < F.length(); ++i)
            {
                outfile << F(i) << endl;
            }
        }
    }
};

#endif // SAVEF_HPP_INCLUDED
