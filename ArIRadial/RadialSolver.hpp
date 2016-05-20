/*


    Copyright (C) 2014 K. Carter Dodd (carter.dodd@gmail.com)


    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of
    the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

*/

#ifndef RADIALSOLVER_HPP_INCLUDED
#define RADIALSOLVER_HPP_INCLUDED

#include "VectorN.hpp"
#include "dynarray.hpp"

/**

	Solves neutral flow in cylindrical geometry due to ionization.

	numr : number of radial points
	numz : number of vertical points
	numv : number of points in velocity space
*/
template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
class RadialSolver
{

private:
    static const unsigned int numr_bd = numr-1;
    static const unsigned int numz_bd = 2*(numz-1);
    static const unsigned int numz_tot = 2*numz-1;
    static const unsigned int num_bd = 2*(numr_bd + numz_bd);
    static const unsigned int total_numv = 2*numv-1;

    static const unsigned int bdOUTER = 0;
    static const unsigned int bdINNER = 1;
    static const unsigned int bdTOP = 2;
    static const unsigned int bdBOTTOM = 3;

    double mMinr; ///< Inner radiusŔ
    double mMaxr; ///< Outer radius
    double mMaxz; ///< 1/2 height of chamber. Center is z=0.

    double mMass; ///< mass of single neutral atom

    double mT0; ///< temperature of neutrals
    double mVth;

    double mLi[numr]; ///< ionization loss rate in 1/s

    dynarray<double> mM;
    dynarray<double> mN;

    double mUpDownSymmetry;
    double mLeftRightSymmetry;

    double mSolutionN[numr][numz_tot];
    double mSolutionVr[numr][numz_tot];
    double mSolutionVz[numr][numz_tot];


    // outer=0, inner=1, top=2, bottom=3
    double mSolN_BD_top[numr_bd];
    double mSolN_BD_bottom[numr_bd];
    double mSolN_BD_outer[numz_bd];
    double mSolN_BD_inner[numz_bd];

    double mSolN_BD_top_tmp[numr_bd];
    double mSolN_BD_bottom_tmp[numr_bd];
    double mSolN_BD_outer_tmp[numz_bd];
    double mSolN_BD_inner_tmp[numz_bd];

    bool mBD_matrix_compiled;
    double mBD_matrix[num_bd][num_bd];

    double mDr;
    double mDz;

    double mDtheta;
    double mDphi;

public:

    /**
        Constructor for solver.

        @param minr Inner radius of chamber
        @param maxr Outer radius of chamber
        @param maxz 1/2 of chamber height: z=0 is in center of chamber
        @param mass Neutral atom mass
        @param T0 Thermal temperature of neutrals
        @param Li Ionization loss rate of neutral atoms in units of 1/s as a function of radius, starting at inner radius
                and ending at outer radius.
        @param M_file path to file containing M integral
        @param N_file path to file containing N integral
        @param up_down_symmetry 1.0=perfectly up/down symmetric, 0.0 = completely asymmetric
        @param left_right_symmetry only affected if up_down_symmetry < 1.0. left_right_symmetry = 0.5 = left/right is perfectly symmetric, 0.0 it's left heavy, 1.0 it's right heavy
    */
    RadialSolver(
        double minr,
        double maxr,
        double maxz,
        double mass,
        double T0,
        const double (&Li)[numr],
        const std::string& M_file,
        const std::string& N_file,
        double up_down_symmetry,
        double left_right_symmetry);

    double solutionN(unsigned int nr, unsigned int nz);

    double solutionVr(unsigned int nr, unsigned int nz);

    double solutionVz(unsigned int nr, unsigned int nz);

    double boundaryN(unsigned int n, unsigned int bd);

    /**

    */
    double avg_ion_rate();

    double avg_n();


    void solve_bd(
        double precision,
        unsigned int max_iterations,
        bool multithread);

    /**
         Solves the neutral profile on the boundary

	 bd: outer=0, inner=1, top=2, bottom=3
    */
    void solve_bd(
        unsigned int n,
        unsigned int bd);

    /** 
	Computes coefficients for boundary solution

	bd: outer=0, inner=1, top=2, bottom=3
    */
    void compile_bd(
        unsigned int n,
        unsigned int bd);

    /**
	Solves neutral flow in volume of cylinder
    */
    void solve(bool multithread);

    /**
	Solves neutral flow at particular point
    */
    void solve(
        unsigned int nr,
        unsigned int nz);

    /**
	Integrates velocity distribution
    */
    void integrate(
        double r,
        double z,
        double theta,
        double phi,
        double& nM,
        double& nN,
        unsigned int& bdi);

    /**
	computes index for boundary element
    */
    unsigned int bd_index(
        unsigned int n,
        unsigned int bd);

    double neTe_profile_function(double r, double rp, double gamma);
};

#endif // ARIRADIALSOLVER_HPP_INCLUDED
