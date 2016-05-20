#ifndef LOADFUNCTION1D_HPP_INCLUDED
#define LOADFUNCTION1D_HPP_INCLUDED

#include <string>
#include <iostream>
#include <fstream>


using namespace std;

template<class T>
class LoadFunction1D
{
private:

    string mFilename;

public:

    LoadFunction1D(const string& filename)
    :
        mFilename(filename)
    {

    }

    void load(
        dynarry<T>& function)
    {
        ifstream inputfile(mFilename.c_str());


        if (inputfile.is_open())
        {
            unsigned int input_num = 0;

            inputfile >> input_num;

            if (input_numr != numr)
            {
                cerr << "Input error: numr not correct in input file." << endl;
            }

            double tmp;

            // r
            for(unsigned int nr = 0; nr < numr; ++nr)
            {
                inputfile >> tmp;

                //cerr << tmp << endl;
            }

            // Te
            for(unsigned int nr = 0; nr < numr; ++nr)
            {
                inputfile >> tmp;
            }

            // ne
            for(unsigned int nr = 0; nr < numr; ++nr)
            {
                inputfile >> tmp;
            }

            for(unsigned int nr = 0; nr < numr; ++nr)
            {

                inputfile >> Li[nr];

                //cerr << Li[nr] << endl;

            }

            inputfile.close();
        }

    }

};

#endif // LOADATOMICDATA_HPP_INCLUDED
