#include <iostream>

#include "RadialSolver.hpp"
#include "RadialSolver.cpp"

#include <string>
#include <iostream>
#include <fstream>

#include "ExportSolution.hpp"
#include "LoadAtomicData.hpp"

using namespace std;

int main()
{
    const unsigned int numr = 101;
    const unsigned int numz = 101;
    const unsigned int numr_bd = numr-1;
    const unsigned int numz_bd = 2*(numz-1);
    const unsigned int numv = 50;

    unsigned int num_files = 5;
    string inputfiles[] = {"ion_ArI_8p3eV.txt", "ion_ArI_6p5eV.txt", "ion_ArI_7p5eV.txt", "ion_ArI_10eV.txt", "ion_ArI_12p5eV.txt"};
    string outputfiles[] = {"sol_ArI_8p3eV_86.txt", "sol_ArI_6p5eV_86.txt", "sol_ArI_7p5eV_86.txt", "sol_ArI_10eV_86.txt", "sol_ArI_12p5eV_86.txt"};

    double zmax[] = {0.86, 0.86, 0.86, 0.86, 0.86};

    //double dev = -1;
    //double ion_precision = 0.01;
    double bd_precision = 1E-6;


    for(unsigned int i = 0; i < num_files; ++i)
    {
        cerr << "Calculating solution for input file: " << inputfiles[i] << "." << endl;

        LoadAtomicData<numr> loader_atomic_data(inputfiles[i]);

        double Li[numr];

        loader_atomic_data.load(Li);

        RadialSolver<numr,numz,numv> solver_high(0.6, 1.6, zmax[i], 40*1.67e-27, 300, Li, "M.txt", "N.txt", 0.0, 0.75);

        cerr << "Calculating Boundary." << endl;

        solver_high.solve_bd(bd_precision, 1000, true);

        cerr << "Calculating Solution..." << endl;

        solver_high.solve(true);

        cerr << "Saving solution to: " << outputfiles[i] << "." << endl;

        ExportSolution<numr,numz,numr_bd,numz_bd,numv> saver(solver_high);

        saver.save(outputfiles[i]);

    }

    return 0;
}
