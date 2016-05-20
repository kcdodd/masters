#ifndef EXPORTSOLUTION_HPP_INCLUDED
#define EXPORTSOLUTION_HPP_INCLUDED

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
class ExportSolution
{
public:
    typedef RadialSolver<numr,numz,numv> SolverT;

private:
    SolverT& mSolver;

public:
    ExportSolution(SolverT& solver)
    :
        mSolver(solver)
    {}

    void save(const std::string& filename)
    {
        ofstream outfile(filename.c_str());

        if (outfile.is_open())
        {
            // numr
            outfile << numr << endl;
            // numz
            outfile << 2*numz-1 << endl;
            // numr_bd
            outfile << numr_bd << endl;
            // numz_bd
            outfile << numz_bd << endl;
            // numv
            outfile << numv << endl;
            // avg n
            outfile << mSolver.avg_n() << endl;
            // avg_ion_rate
            outfile << mSolver.avg_ion_rate() << endl;



                // export outer boundary
                for (unsigned int nz = 0; nz < numz_bd; ++nz)
                {
                    outfile << mSolver.boundaryN(nz, 0) << "\t";
                }

                outfile << endl;

                // inner boundary
                for (unsigned int nz = 0; nz < numz_bd; ++nz)
                {
                    outfile << mSolver.boundaryN(nz, 1) << "\t";
                }

                outfile << endl;

                // top boundary
                for(unsigned int nr = 0; nr < numr_bd; ++nr)
                {
                    outfile << mSolver.boundaryN(nr, 2) << "\t";
                }

                outfile << endl;

                // bottom boundary
                for(unsigned int nr = 0; nr < numr_bd; ++nr)
                {
                    outfile << mSolver.boundaryN(nr, 3) << "\t";
                }

                outfile << endl;

                for(unsigned int nz = 0; nz < 2*numz-1; ++nz)

                {
                    for(unsigned int nr = 0; nr < numr; ++nr)
                    {

                        outfile << mSolver.solutionN(nr, nz) << "\t";

                    }//r

                    outfile << endl;

                }//z

                outfile << endl;


                for(unsigned int nz = 0; nz < 2*numz-1; ++nz)

                {
                    for(unsigned int nr = 0; nr < numr; ++nr)
                    {

                        outfile << mSolver.solutionVr(nr, nz) << "\t";

                    }//r

                    outfile << endl;

                }//z

                outfile << endl;



                for(unsigned int nz = 0; nz < 2*numz-1; ++nz)

                {
                    for(unsigned int nr = 0; nr < numr; ++nr)
                    {

                        outfile << mSolver.solutionVz(nr, nz) << "\t";

                    }//r

                    outfile << endl;

                }//z

                outfile << endl;


            outfile.close();
        }
    }
};

#endif // EXPORTSOLUTION_HPP_INCLUDED
