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

#include "RadialSolver.hpp"
#include "Constants.hpp"

#include <cmath>
#include <algorithm>

#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

#include <iostream>

#include "dynarray.hpp"

#include <time.h>

using namespace std;

//----------------------------------------------------------------------------------------------------
// RadialSolver()
//----------------------------------------------------------------------------------------------------

template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
RadialSolver<
    numr,
    numz,
    numv>::RadialSolver(
        double minr,
        double maxr,
        double maxz,
        double mass,
        double T0,
        const double (&Li)[numr],
        const std::string& M_file,
        const std::string& N_file,
        double up_down_symmetry,
        double left_right_symmetry)
    :
        mMinr(minr),
        mMaxr(maxr),
        mMaxz(maxz),
        mMass(mass),
        mT0(T0),
        mVth(sqrt(2*BOLTZMANN_CONST*mT0/mMass)),
        mBD_matrix_compiled(false),
        mUpDownSymmetry(up_down_symmetry),
        mLeftRightSymmetry(left_right_symmetry)
    {
        for(unsigned int nr = 0; nr < numr_bd; ++nr)
        {

            mSolN_BD_top[nr] = 0.5;
            mSolN_BD_top_tmp[nr] = 0.5;

            mSolN_BD_bottom[nr] = 0.5;
            mSolN_BD_bottom_tmp[nr] = 0.5;
        }

        for (unsigned int nz = 0; nz < numz_bd; ++nz)
        {
            mSolN_BD_outer[nz]=0.5;
            mSolN_BD_outer_tmp[nz]=0.5;

            mSolN_BD_inner[nz]=0.5;
            mSolN_BD_inner_tmp[nz]=0.5;
        }

        for(unsigned int nr = 0; nr < numr; ++nr)
        {

            mLi[nr] = Li[nr];

            for (unsigned int nz = 0; nz < numz_tot; ++nz)
            {

                    mSolutionN[nr][nz] = 1;
                    mSolutionVr[nr][nz] = 0;
                    mSolutionVz[nr][nz] = 0;

            }
        }

        for (unsigned int i=0; i < num_bd; ++i)
        {
            for (unsigned int j=0; j < num_bd; ++j)
            {
                mBD_matrix[i][j] = 0;
            }
        }

        mM.load_text_file(M_file);
        mN.load_text_file(N_file);

        mDtheta = PI_CONST/(total_numv-1);
        mDphi = PI_CONST/(total_numv-1);

        mDr = (mMaxr - mMinr)/(numr-1);

        mDz = mMaxz/(numz-1);


    }

//----------------------------------------------------------------------------------------------------
// solutionN()
//----------------------------------------------------------------------------------------------------
template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
double RadialSolver<
    numr,
    numz,
    numv>::solutionN(unsigned int nr, unsigned int nz)
    {
        return mSolutionN[nr][nz];
    }

//----------------------------------------------------------------------------------------------------
// solutionVr()
//----------------------------------------------------------------------------------------------------
template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
double RadialSolver<
    numr,
    numz,
    numv>::solutionVr(unsigned int nr, unsigned int nz)
    {
        return mSolutionVr[nr][nz];
    }

//----------------------------------------------------------------------------------------------------
// solutionVz()
//----------------------------------------------------------------------------------------------------
template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
double RadialSolver<
    numr,
    numz,
    numv>::solutionVz(unsigned int nr, unsigned int nz)
    {
        return mSolutionVz[nr][nz];
    }

//----------------------------------------------------------------------------------------------------
// boundaryN()
//----------------------------------------------------------------------------------------------------

template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
double RadialSolver<
    numr,
    numz,
    numv>::boundaryN(
        unsigned int n,
        unsigned int bd)
    {
        // outer
        if (bd == bdOUTER)
        {
            return mSolN_BD_outer[n];
        }

        // inner
        if (bd == bdINNER)
        {
            return mSolN_BD_inner[n];
        }

        // top
        if (bd == bdTOP)
        {
            return mSolN_BD_top[n];
        }

        // bottom
        if (bd == bdBOTTOM)
        {
            return mSolN_BD_bottom[n];
        }

        return -1;
    }

//----------------------------------------------------------------------------------------------------
// ionLevel()
//----------------------------------------------------------------------------------------------------

template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
double RadialSolver<
    numr,
    numz,
    numv>::avg_ion_rate()
    {


        double sum=0;
        double weights=0;
        double r_weight = 0;

        for (unsigned int nr = 1; nr < numr-1; ++nr)
        {
            r_weight = nr*mDr+mMinr;

            for(unsigned int nz=1; nz < numz_tot-1; ++nz)
            {
                sum += mLi[nr]*mSolutionN[nr][nz]*r_weight;
                weights += mSolutionN[nr][nz]*r_weight;
            }
        }

        double avg = sum/weights;

        return avg;

    }


template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
double RadialSolver<
    numr,
    numz,
    numv>::avg_n()
    {


        double sum=0;
        double weights=0;
        double r_weight = 0;

        for (unsigned int nr = 1; nr < numr-1; ++nr)
        {
            r_weight = nr*mDr+mMinr;

            for(unsigned int nz=1; nz < numz_tot-1; ++nz)
            {
                sum += mSolutionN[nr][nz]*r_weight;
                weights += r_weight;
            }
        }

        double avg = sum/weights;

        return avg;

    }

//----------------------------------------------------------------------------------------------------
// solve_bd()
//----------------------------------------------------------------------------------------------------

template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
void RadialSolver<
    numr,
    numz,
    numv>::solve_bd(
        double precision,
        unsigned int max_iterations,
        bool multithread)
    {

        if (!mBD_matrix_compiled)
        {
            cerr << "First time solving boundary: compiling boundary matrix." << endl;

            // must first compile boundary matrix

            // set boundary values to unity to compile matrix so that values can be scaled to solution
            for(unsigned int nr = 0; nr < numr_bd; ++nr)
            {
                mSolN_BD_top[nr] = 1;
                mSolN_BD_bottom[nr] = 1;

            }

            for (unsigned int nz = 0; nz < numz_bd; ++nz)
            {
                mSolN_BD_outer[nz]=1;
                mSolN_BD_inner[nz]=1;
            }

            if (multithread)
            {
                // solve in parallel

                boost::thread_group calcThreads;

                for (unsigned int nz = 0; nz < numz_bd; ++nz)
                {
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numv>::compile_bd, this, nz, bdOUTER));
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numv>::compile_bd, this, nz, bdINNER));
                }

                for (unsigned int nr = 0; nr < numr_bd; ++nr)
                {
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numv>::compile_bd, this, nr, bdTOP));
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numv>::compile_bd, this, nr, bdBOTTOM));
                }

                calcThreads.join_all();

            }else{

                // solve serially
                for (unsigned int nz = 0; nz < numz_bd; ++nz)
                {
                    compile_bd(nz, bdOUTER);
                    compile_bd(nz, bdINNER);
                }

                for (unsigned int nr = 0; nr < numr_bd; ++nr)
                {
                    compile_bd(nr, bdTOP);
                    compile_bd(nr, bdBOTTOM);
                }
            }

            mBD_matrix_compiled = true;
            cerr << "Boundary matrix compiled." << endl;
        }

        unsigned int iterations = 0;

        double maxdev;

        // keep track of how long the solution takes to compute
        time_t t_begin = time(0);

        do{
            //cerr << "Solving boundary: ";

            for(unsigned int nr = 0; nr < numr_bd; ++nr)
            {
                mSolN_BD_top_tmp[nr] = mSolN_BD_top[nr];
                mSolN_BD_bottom_tmp[nr] = mSolN_BD_bottom[nr];
            }

            for (unsigned int nz = 0; nz < numz_bd; ++nz)
            {
                mSolN_BD_outer_tmp[nz]=mSolN_BD_outer[nz];
                mSolN_BD_inner_tmp[nz]=mSolN_BD_inner[nz];
            }

            // solve boundary

            if (multithread)
            {
                // solve in parallel

                boost::thread_group calcThreads;

                for (unsigned int nz = 0; nz < numz_bd; ++nz)
                {
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numv>::solve_bd, this, nz, bdOUTER));
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numv>::solve_bd, this, nz, bdINNER));
                }

                for (unsigned int nr = 0; nr < numr_bd; ++nr)
                {
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numv>::solve_bd, this, nr, bdTOP));
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numv>::solve_bd, this, nr, bdBOTTOM));
                }

                calcThreads.join_all();

            }else{

                // solve serially
                for (unsigned int nz = 0; nz < numz_bd; ++nz)
                {
                    solve_bd(nz, bdOUTER);
                    solve_bd(nz, bdINNER);
                }

                for (unsigned int nr = 0; nr < numr_bd; ++nr)
                {
                    solve_bd(nr, bdTOP);
                    solve_bd(nr, bdBOTTOM);
                }
            }

            maxdev=0;

            for(unsigned int nr = 0; nr < numr_bd; ++nr)
            {
                double dev = 2*abs(mSolN_BD_top_tmp[nr] - mSolN_BD_top[nr])/(mSolN_BD_top_tmp[nr] + mSolN_BD_top[nr]);

                maxdev = max(maxdev, dev);

                dev = 2*abs(mSolN_BD_bottom_tmp[nr] - mSolN_BD_bottom[nr])/(mSolN_BD_bottom_tmp[nr] + mSolN_BD_bottom[nr]);

                maxdev = max(maxdev, dev);
            }

            for (unsigned int nz = 0; nz < numz_bd; ++nz)
            {
                double dev = 2*abs(mSolN_BD_outer_tmp[nz] - mSolN_BD_outer[nz])/(mSolN_BD_outer_tmp[nz] + mSolN_BD_outer[nz]);

                maxdev = max(maxdev, dev);

                dev = 2*abs(mSolN_BD_inner_tmp[nz] - mSolN_BD_inner[nz])/(mSolN_BD_inner_tmp[nz] + mSolN_BD_inner[nz]);

                maxdev = max(maxdev, dev);

            }

            iterations++;

            //time_t t_cur = time(0);

            //cerr << maxdev << " deviation, iteration " << iterations << ", time = " << t_cur - t_begin << "s." << endl;

        }while(maxdev > precision && iterations < max_iterations);


        time_t t_end = time(0);

        if (iterations == max_iterations)
        {
            cerr << "Max iterations reached in attempt to solve boundary." << endl;
        }else{
            cerr << "Solution for boundary found to " << maxdev << " deviation in " << iterations << " steps and " << t_end - t_begin << "s." << endl;
        }


    }


//----------------------------------------------------------------------------------------------------
// solve_bd()
//----------------------------------------------------------------------------------------------------

template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
void RadialSolver<
    numr,
    numz,
    numv>::solve_bd(
        unsigned int n,
        unsigned int bd)
    {

        unsigned int this_bd_index = bd_index(n, bd);

        double out_flux = 0;

        for(unsigned int nr = 0; nr < numr_bd; ++nr)
        {
            out_flux += mSolN_BD_top_tmp[nr]*mBD_matrix[this_bd_index][bd_index(nr, bdTOP)];
            out_flux += mSolN_BD_bottom_tmp[nr]*mBD_matrix[this_bd_index][bd_index(nr, bdBOTTOM)];
        }

        for (unsigned int nz = 0; nz < numz_bd; ++nz)
        {
            out_flux += mSolN_BD_outer_tmp[nz]*mBD_matrix[this_bd_index][bd_index(nz, bdOUTER)];
            out_flux += mSolN_BD_inner_tmp[nz]*mBD_matrix[this_bd_index][bd_index(nz, bdINNER)];
        }




        if (bd == bdTOP)
        {
            double r = n*mDr+mMinr;

            double source_symmetric = pow(neTe_profile_function(r, 1.0, 0.4), 0.5)*neTe_profile_function(r, 1.1, 0.5)*r*2/(mMinr + mMaxr);
            double source_asymmetric = max(-sin(TWOPI_CONST*(r-mMinr)/(mMaxr-mMinr)), 0.0);


            mSolN_BD_top[n] = out_flux*sqrt(PI_CONST)/mVth + source_symmetric /* +  mUpDownSymmetry*source_symmetric + (1.0-mUpDownSymmetry)*mLeftRightSymmetry*source_asymmetric*/;

        }else if (bd == bdBOTTOM){

            double r = n*mDr+mMinr;

            double source_symmetric = pow(neTe_profile_function(r, 1.0, 0.4), 0.5)*neTe_profile_function(r, 1.1, 0.5)*r*2/(mMinr + mMaxr);
            double source_asymmetric = max(sin(TWOPI_CONST*(r-mMinr)/(mMaxr-mMinr)), 0.0);


            mSolN_BD_bottom[n] = out_flux*sqrt(PI_CONST)/mVth/* +  mUpDownSymmetry*source_symmetric + (1.0-mUpDownSymmetry)*(1.0-mLeftRightSymmetry)*source_asymmetric*/;


        }else if (bd == bdINNER){


            mSolN_BD_inner[n] = out_flux*sqrt(PI_CONST)/mVth;


        }else{


            mSolN_BD_outer[n] = out_flux*sqrt(PI_CONST)/mVth;


        }


    }

//----------------------------------------------------------------------------------------------------
// compile_bd()
//----------------------------------------------------------------------------------------------------

template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
void RadialSolver<
    numr,
    numz,
    numv>::compile_bd(
        unsigned int n,
        unsigned int bd)
    {

        double r, z;
        unsigned int min_theta, max_theta, min_phi, max_phi;
        double total_theta;
        bool need_vr_vz; // vr == true, vz == false
        double dot_product_sign;


        // outer
        if (bd == bdOUTER)
        {
            r = mMaxr;
            z = n*mDz - mMaxz + 0.5*mDz;

            need_vr_vz = true;

            min_theta = 0;
            max_theta = numv;
            total_theta = PI_CONST;

            dot_product_sign = 1;

            min_phi = 0;
            max_phi = total_numv;


        }

        // inner
        if (bd == bdINNER)
        {

            r = mMinr;
            z = n*mDz - mMaxz + 0.5*mDz;;

            need_vr_vz = true;

            min_theta = numv-1;
            max_theta = total_numv;

            total_theta = PI_CONST;

            dot_product_sign = -1;


            min_phi = 0;
            max_phi = total_numv;

        }

        // top
        if (bd == bdTOP)
        {
            r = n*mDr + mMinr + 0.5*mDr;
            z = mMaxz;

            need_vr_vz = false;

            min_phi = 0;
            max_phi = numv;

            dot_product_sign = 1;


            min_theta = 0;
            max_theta = total_numv;

            total_theta = TWOPI_CONST;


        }

        // bottom
        if (bd == bdBOTTOM)
        {
            r = n*mDr + mMinr + 0.5*mDr;
            z = -mMaxz;

            need_vr_vz = false;

            min_phi = numv-1;
            max_phi = total_numv;

            dot_product_sign = -1;


            min_theta = 0;
            max_theta = total_numv;

            total_theta = TWOPI_CONST;


        }

        double I1 = 0;
        double I2 = 0;
        double I = 0;
        unsigned int bdi;

        double int_fact;
        double int_edge_fact;

        double V_factor = 2.0*mVth/pow(PI_CONST, 3.0/2.0);


        double theta, phi;

        unsigned int this_bd_index = bd_index(n, bd);

        for(unsigned int phi_i = min_phi; phi_i < max_phi; ++phi_i)
        {
            if (phi_i == 0)
            {
                int_fact = total_theta*(1-sin(mDphi)/mDphi);

                integrate(r, z, 0, 0, I1, I2, bdi);


                if (need_vr_vz)
                {
                    I = 0;
                }else{
                    I = dot_product_sign*int_fact*V_factor*I2;
                }

                mBD_matrix[this_bd_index][bdi] += I;


            }else if (phi_i == total_numv-1){

                int_fact = total_theta*(1-sin(mDphi)/mDphi);

                integrate(r, z, 0, PI_CONST, I1, I2, bdi);

                if (need_vr_vz)
                {
                    I = 0;
                }else{
                    I = -dot_product_sign*int_fact*V_factor*I2;
                }

                mBD_matrix[this_bd_index][bdi] += I;

            }else{

                phi = phi_i*mDphi;

                if (phi_i == min_phi)
                {
                    // if it got here then min_phi != 0, so integration factor is different

                    int_fact = mDtheta*(cos(phi) - (sin(phi+mDphi)-sin(phi))/mDphi);

                }else if (phi_i == max_phi){
                    // same deal except for other limit

                    int_fact = mDtheta*((sin(phi)-sin(phi-mDphi))/mDphi - cos(phi));

                }else{

                    // not on edge of phi integration
                    int_fact = 2*mDtheta*sin(phi)*(1-cos(mDphi))/mDphi;
                }



                for (unsigned int theta_i = min_theta; theta_i < max_theta; ++theta_i)
                {
                    if (theta_i == 0)
                    {
                        integrate(r, z, 0, phi, I1, I2, bdi);

                        if (need_vr_vz)
                        {
                            I = dot_product_sign*sin(phi)*int_fact*V_factor*I2;
                        }else{
                            I = dot_product_sign*cos(phi)*int_fact*V_factor*I2;
                        }

                        mBD_matrix[this_bd_index][bdi] += I;

                    }else if (theta_i == total_numv-1){

                        integrate(r, z, PI_CONST, phi, I1, I2, bdi);

                        if (need_vr_vz)
                        {
                            I = -dot_product_sign*sin(phi)*int_fact*V_factor*I2;
                        }else{
                            I = dot_product_sign*cos(phi)*int_fact*V_factor*I2;
                        }

                        mBD_matrix[this_bd_index][bdi] += I;

                    }else{

                        if (theta_i == min_theta || min_theta == max_theta)
                        {
                            // if limits are not 0 and pi, then need to alter integration factor
                            int_edge_fact = 1/2;
                        }else{
                            int_edge_fact = 1;
                        }

                        theta = theta_i*mDtheta;

                        integrate(r, z, theta, phi, I1, I2, bdi);

                        if (need_vr_vz)
                        {
                            I = 2*dot_product_sign*int_edge_fact*sin(phi)*cos(theta)*int_fact*V_factor*I2;
                        }else{
                            I = 2*dot_product_sign*int_edge_fact*cos(phi)*int_fact*V_factor*I2;
                        }

                        mBD_matrix[this_bd_index][bdi] += I;


                    }

                }
            }

        }


    }

//----------------------------------------------------------------------------------------------------
// solve()
//----------------------------------------------------------------------------------------------------

template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
void RadialSolver<
    numr,
    numz,
    numv>::solve(bool multithread)
    {
        if (multithread)
        {
            // solve in parallel

            boost::thread_group calcThreads;

            for (unsigned int nr = 1; nr < numr-1; ++nr)
            {
                for (unsigned int nz = 1; nz < numz_tot-1; ++nz)
                {
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numv>::solve, this, nr, nz));
                }

            }

            calcThreads.join_all();

        }else{
            // solve serially

            for (unsigned int nr = 1; nr < numr-1; ++nr)
            {
                for (unsigned int nz = 1; nz < numz_tot-1; ++nz)
                {

                    solve(nr, nz);
                }

            }
        }
    }

//----------------------------------------------------------------------------------------------------
// solve()
//----------------------------------------------------------------------------------------------------

template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
void RadialSolver<
    numr,
    numz,
    numv>::solve(
        unsigned int nr,
        unsigned int nz)
    {



        double r = nr*mDr + mMinr;
        double z = nz*mDz - mMaxz;

        // reset to zero

        mSolutionN[nr][nz] = 0;
        mSolutionVr[nr][nz] = 0;
        mSolutionVz[nr][nz] = 0;

        double I1 = 0;
        double I2 = 0;
        unsigned int bdi;

        double int_fact;
        double N_factor = 2.0/pow(PI_CONST, 3.0/2.0);
        double V_factor = 2.0*mVth/pow(PI_CONST, 3.0/2.0);
        double theta, phi;

        for(unsigned int phi_i = 0; phi_i < total_numv; ++phi_i)
        {
            if (phi_i == 0)
            {
                int_fact = TWOPI_CONST*(1-sin(mDphi)/mDphi);

                integrate(r, z, 0, 0, I1, I2, bdi);

                mSolutionN[nr][nz] += int_fact*N_factor*I1;
                mSolutionVz[nr][nz] += int_fact*V_factor*I2;

            }else if (phi_i == total_numv-1){

                int_fact = TWOPI_CONST*(1-sin(mDphi)/mDphi);

                integrate(r, z, 0, PI_CONST, I1, I2, bdi);

                mSolutionN[nr][nz] += int_fact*N_factor*I1;
                mSolutionVz[nr][nz] += -int_fact*V_factor*I2;

            }else{


                phi = phi_i*mDphi;

                int_fact = 2*mDtheta*sin(phi)*(1-cos(mDphi))/mDphi;

                for (unsigned int theta_i = 0; theta_i < total_numv; ++theta_i)
                {
                    if (theta_i == 0)
                    {
                        integrate(r, z, 0, phi, I1, I2, bdi);

                        mSolutionN[nr][nz] += int_fact*N_factor*I1;
                        mSolutionVr[nr][nz] += sin(phi)*int_fact*V_factor*I2;
                        mSolutionVz[nr][nz] += cos(phi)*int_fact*V_factor*I2;

                    }else if (theta_i == total_numv-1){

                        integrate(r, z, PI_CONST, phi, I1, I2, bdi);

                        mSolutionN[nr][nz] += int_fact*N_factor*I1;
                        mSolutionVr[nr][nz] += -sin(phi)*int_fact*V_factor*I2;
                        mSolutionVz[nr][nz] += cos(phi)*int_fact*V_factor*I2;

                    }else{

                        theta = theta_i*mDtheta;

                        integrate(r, z, theta, phi, I1, I2, bdi);

                        mSolutionN[nr][nz] += 2*int_fact*N_factor*I1;
                        mSolutionVr[nr][nz] += 2*sin(phi)*cos(theta)*int_fact*V_factor*I2;
                        mSolutionVz[nr][nz] += 2*cos(phi)*int_fact*V_factor*I2;
                    }

                }
            }

        }


    }

//----------------------------------------------------------------------------------------------------
// integrate()
//----------------------------------------------------------------------------------------------------
template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
void RadialSolver<
    numr,
    numz,
    numv>::integrate(
        double r,
        double z,
        double theta,
        double phi,
        double& nM,
        double& nN,
        unsigned int& bdi)
    {
        // setup variables defining integral path


        double sintheta = sin(theta);
        double costheta = cos(theta);

        double sinphi = sin(phi);
        double cosphi = cos(phi);


        double theta_transition = asin(mMinr/r);

        double L_perp = -1;

        if (theta <= theta_transition)
        {
            //inner wall
            L_perp = r*costheta - sqrt(mMinr*mMinr-r*r*sintheta*sintheta);


        }else{
            // outer wall
            L_perp = r*costheta + sqrt(mMaxr*mMaxr-r*r*sintheta*sintheta);
        }

        double h1 = mMaxz+z;
        double h2 = mMaxz-z;

        double phi1 = atan(L_perp/h1);
        double phi2 = atan(h2/L_perp) + HALFPI_CONST;

        double L = -1;

        if (z == mMaxz && phi == HALFPI_CONST)
        {
            L = L_perp;

        }else if ((r == mMinr || r == mMaxr) && phi == 0){

            L = h1;

        }else if ((r == mMinr || r == mMaxr) && phi == PI_CONST){

            L = h2;

        }else{

            if (phi < phi1)
            {
                // floor
                L = h1/cosphi;

            }else if (phi < phi2){
                // wall
                L = L_perp/sinphi;

            }else{
                // ceiling
                L = -h2/cosphi;
            }
        }


        // get value at wall using linear interpolation
        double n0_wall = 0;

        double r_intersect = min(max(sqrt(r*r + L*L*sinphi*sinphi - 2*r*L*sinphi*costheta), mMinr), mMaxr);
        double delta_r = r - r_intersect;
        double delta_r2 = delta_r*delta_r;
        double delta_z = -L*cosphi;
        double delta_z2 = delta_z*delta_z;
        double z_intersect = min(max(z + delta_z, -mMaxz), mMaxz);


        double frac1, frac2;
        unsigned int nri, nzi;
        unsigned int nri1, nri2;

        if (phi < phi1 || phi > phi2)
        {
            // floor or ceiling

            nri = floor((r_intersect - mMinr)/mDr);

            if (phi < phi1)
            {
                // bottom

                n0_wall = mSolN_BD_bottom[nri];

                bdi = bd_index(nri, bdBOTTOM);

            }else{

                // top

                n0_wall = mSolN_BD_top[nri];

                bdi = bd_index(nri, bdTOP);
            }

        }else{
            // wall

            nzi = floor((z_intersect+mMaxz)/mDz);


            if (theta <= theta_transition)
            {
                //inner wall

                n0_wall = mSolN_BD_inner[nzi];

                bdi = bd_index(nzi, bdINNER);

            }else{

                // outer wall

                n0_wall = mSolN_BD_outer[nzi];

                bdi = bd_index(nzi, bdOUTER);

            }

        }

        if (n0_wall < 0)
        {
            cerr << "n0_wall not defined" << endl;
        }


        unsigned int N = 200;
        double dl = L/((double) N-1);
        double l = 0;
        double ri;
        double x;

        if (L == 0.0)
        {
            cerr << "L = 0" << endl;
        }

        dynarray<double> ion_rate_function(N, 0.0, L);


        for (unsigned int i=0; i < N; ++i)
        {
            // calculate element indecies
            l = i*dl;
            x = l - L;


            ri = sqrt(r*r + x*x*sinphi*sinphi + 2*r*x*sinphi*costheta);

            if (ri > mMaxr)
            {
                //cerr << "r > " << mMaxr << ": " << ri << endl;

                ri = mMaxr;
            }

            if (ri < mMinr)
            {
                //cerr << "r < " << mMinr << ": " << ri << endl;

                ri = mMinr;
            }

            nri1 = floor((ri - mMinr)/mDr);
            nri2 = nri1+1;
            frac2 = (ri - (mDr*nri1 + mMinr))/mDr;
            frac1 = 1-frac2;

            ion_rate_function(i) = frac1*mLi[nri1] + frac2*mLi[nri2];

            //ion_rate_function(i) = 0;

        }

        double vk = ion_rate_function.integral_total();

        double alpha = vk/mVth;


        if (alpha > mM.xmax())
        {
            nM = 0.0;
            nN = 0.0;
        }else{
            nM = n0_wall*mM.interp(alpha);
            nN = n0_wall*mN.interp(alpha);
        }




    }



template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
unsigned int RadialSolver<
    numr,
    numz,
    numv>::
    bd_index(
        unsigned int n,
        unsigned int bd)
    {
        if (bd == bdOUTER)
        {
            return n;
        }

        if (bd == bdINNER)
        {
            return n + numz_bd;
        }

        if (bd == bdTOP)
        {
            return n + 2*numz_bd;
        }

        if (bd == bdBOTTOM)
        {
            return n + 2*numz_bd + numr_bd;
        }

        throw -1;
    }




template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numv>
double RadialSolver<
    numr,
    numz,
    numv>::
    neTe_profile_function(
        double r,
        double rp,
        double gamma)
    {
        double xp = (rp - mMinr)/(mMaxr - mMinr);
        double x = (r - mMinr)/(mMaxr - mMinr);

        double a = (0.5 - xp)/(xp*xp - xp);
        double y = a*x*x + (1-a)*x;

        return (gamma/2)*(1-cos(PI_CONST*y)) + (1-gamma/2)*pow(sin(PI_CONST*y), 3);
    }
