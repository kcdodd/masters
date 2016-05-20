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
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
    numv>::RadialSolver(
        double minr,
        double maxr,
        double maxz,
        double mass,
        double T0,
        const double (&Li)[numr],
        const std::string& M_file,
        const std::string& N_file)
    :
        mMinr(minr),
        mMaxr(maxr),
        mMaxz(maxz),
        mMass(mass),
        mT0(T0),
        mVth(sqrt(2*BOLTZMANN_CONST*mT0/mMass)),
        mBD_matrix_compiled(false)
    {
        for(unsigned int nr = 0; nr < numr_bd; ++nr)
        {

            mSolN_BD_top[nr] = 0.5;
            mSolN_BD_top_tmp[nr] = 0.5;
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

            for (unsigned int nz = 0; nz < numz; ++nz)
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


        mDr_bd = (mMaxr - mMinr)/(numr_bd-1);
        mDz_bd = mMaxz/(numz_bd-1);

        mDr = (mMaxr - mMinr)/(numr-1);

        mDz = mMaxz/(numz-1);

        for (unsigned int nz = 0; nz < numz; ++nz)
        {
            mZ[nz] = nz*mDz;
            //cerr << mZ[nz] << endl;
        }

    }

//----------------------------------------------------------------------------------------------------
// solutionN()
//----------------------------------------------------------------------------------------------------
template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
double RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
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
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
double RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
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
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
double RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
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
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
double RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
    numv>::boundaryN(
        unsigned int n,
        unsigned int bd)
    {
        // outer
        if (bd == 0)
        {
            return mSolN_BD_outer[n];
        }

        // inner
        if (bd == 1)
        {
            return mSolN_BD_inner[n];
        }

        // top
        if (bd == 2)
        {
            return mSolN_BD_top[n];
        }

        return -1;
    }

//----------------------------------------------------------------------------------------------------
// ionLevel()
//----------------------------------------------------------------------------------------------------

template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
double RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
    numv>::avg_ion_rate()
    {


        double sum=0;
        double weights=0;
        double r_weight = 0;

        for (unsigned int nr = 1; nr < numr-1; ++nr)
        {
            r_weight = nr*mDr+mMinr;

            sum += mLi[nr]*mSolutionN[nr][0]*r_weight;
            weights += mSolutionN[nr][0]*r_weight;

            for(unsigned int nz=1; nz < numz-1; ++nz)
            {
                sum += mLi[nr]*mSolutionN[nr][nz]*2*r_weight;
                weights += mSolutionN[nr][nz]*2*r_weight;
            }
        }

        double avg = sum/weights;

        return avg;

    }


template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
double RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
    numv>::avg_n()
    {


        double sum=0;
        double weights=0;
        double r_weight = 0;

        for (unsigned int nr = 1; nr < numr-1; ++nr)
        {
            r_weight = nr*mDr+mMinr;

            sum += mSolutionN[nr][0]*r_weight;
            weights += r_weight;

            for(unsigned int nz=1; nz < numz-1; ++nz)
            {
                sum += mSolutionN[nr][nz]*2*r_weight;
                weights += 2*r_weight;
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
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
void RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
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
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numr_bd,numz_bd,numv>::compile_bd, this, nz, 0));
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numr_bd,numz_bd,numv>::compile_bd, this, nz, 1));
                }

                for (unsigned int nr = 0; nr < numr_bd; ++nr)
                {
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numr_bd,numz_bd,numv>::compile_bd, this, nr, 2));
                }

                calcThreads.join_all();

            }else{

                // solve serially
                for (unsigned int nz = 0; nz < numz_bd; ++nz)
                {
                    compile_bd(nz, 0);
                    compile_bd(nz, 1);
                }

                for (unsigned int nr = 0; nr < numr_bd; ++nr)
                {
                    compile_bd(nr, 2);
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
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numr_bd,numz_bd,numv>::solve_bd, this, nz, 0));
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numr_bd,numz_bd,numv>::solve_bd, this, nz, 1));
                }

                for (unsigned int nr = 0; nr < numr_bd; ++nr)
                {
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numr_bd,numz_bd,numv>::solve_bd, this, nr, 2));
                }

                calcThreads.join_all();

            }else{

                // solve serially
                for (unsigned int nz = 0; nz < numz_bd; ++nz)
                {
                    solve_bd(nz, 0);
                    solve_bd(nz, 1);
                }

                for (unsigned int nr = 0; nr < numr_bd; ++nr)
                {
                    solve_bd(nr, 2);
                }
            }

            maxdev=0;

            for(unsigned int nr = 0; nr < numr_bd; ++nr)
            {
                double dev = 2*abs(mSolN_BD_top_tmp[nr] - mSolN_BD_top[nr])/(mSolN_BD_top_tmp[nr] + mSolN_BD_top[nr]);

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
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
void RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
    numv>::solve_bd(
        unsigned int n,
        unsigned int bd)
    {

        unsigned int this_bd_index = bd_index(n, bd);

        double out_flux = 0;

        for(unsigned int nr = 0; nr < numr_bd; ++nr)
        {
            out_flux += mSolN_BD_top_tmp[nr]*mBD_matrix[this_bd_index][bd_index(nr, 2)];
        }

        for (unsigned int nz = 0; nz < numz_bd; ++nz)
        {
            out_flux += mSolN_BD_outer_tmp[nz]*mBD_matrix[this_bd_index][bd_index(nz, 0)];
            out_flux += mSolN_BD_inner_tmp[nz]*mBD_matrix[this_bd_index][bd_index(nz, 1)];
        }


        if (bd == 2)
        {
            double r = n*mDr_bd+mMinr;


            double source = pow(neTe_profile_function(r, 1.0, 0.4), 0.5)*neTe_profile_function(r, 1.1, 0.5)*r*2/(mMinr + mMaxr);

            if (n == 0)
            {
                mSolN_BD_top[0] = (out_flux*6*sqrt(PI_CONST)/mVth + 2*mSolN_BD_inner_tmp[numz_bd-1])/5 + source;

            }else if (n == numr_bd-1){

                mSolN_BD_top[numr_bd-1] = (out_flux*6*sqrt(PI_CONST)/mVth + 2*mSolN_BD_inner_tmp[numz_bd-1])/5 + source;
            }else{

                mSolN_BD_top[n] = out_flux*sqrt(PI_CONST)/mVth +  source;

            }
        }else if (bd == 1){

            if (n == numz_bd-1)
            {
                mSolN_BD_inner[numz_bd-1] = (out_flux*6*sqrt(PI_CONST)/mVth + 2*mSolN_BD_top_tmp[0])/5;
            }else{
                mSolN_BD_inner[n] = out_flux*sqrt(PI_CONST)/mVth;
            }

        }else{

            if (n == numz_bd-1)
            {
                mSolN_BD_outer[numz_bd-1] = (out_flux*6*sqrt(PI_CONST)/mVth + 2*mSolN_BD_top_tmp[numr_bd-1])/5;
            }else{
                mSolN_BD_outer[n] = out_flux*sqrt(PI_CONST)/mVth;
            }

        }


    }

//----------------------------------------------------------------------------------------------------
// compile_bd()
//----------------------------------------------------------------------------------------------------

template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
void RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
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
        if (bd == 0)
        {
            r = mMaxr;
            z = n*mDz_bd;

            need_vr_vz = true;

            min_theta = 0;
            max_theta = numv;
            total_theta = PI_CONST;

            dot_product_sign = 1;

            if (n == numz_bd-1)
            {
                // corner only 1/4 of distr.
                min_phi = 0;
                max_phi = numv;

            }else{
                min_phi = 0;
                max_phi = total_numv;
            }


        }

        // inner
        if (bd == 1)
        {

            r = mMinr;
            z = n*mDz_bd;

            need_vr_vz = true;

            min_theta = numv-1;
            max_theta = total_numv;

            total_theta = PI_CONST;

            dot_product_sign = -1;

            if (n == numz_bd-1)
            {
                // corner only 1/4 of distr.
                min_phi = 0;
                max_phi = numv;

            }else{
                min_phi = 0;
                max_phi = total_numv;
            }


        }

        // top
        if (bd == 2)
        {
            r = n*mDr_bd+mMinr;
            z = mMaxz;

            need_vr_vz = false;

            min_phi = 0;
            max_phi = numv;

            dot_product_sign = 1;

            if (n == 0)
            {
                // inner corner
                min_theta = numv-1;
                max_theta = total_numv;

                total_theta = PI_CONST;

            }else if (n == numr_bd-1){
                // outer corner
                min_theta = 0;
                max_theta = numv;

                total_theta = PI_CONST;

            }else{

                min_theta = 0;
                max_theta = total_numv;

                total_theta = TWOPI_CONST;

            }

        }

        double I1 = 0;
        double I2 = 0;
        double I = 0;
        bd_intersect intersect;

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

                integrate(r, z, 0, 0, I1, I2, intersect);


                if (need_vr_vz)
                {
                    I = 0;
                }else{
                    I = dot_product_sign*int_fact*V_factor*I2;
                }

                mBD_matrix[this_bd_index][bd_index(intersect.n_1, intersect.bd_1)] += intersect.fraction_1*I;
                mBD_matrix[this_bd_index][bd_index(intersect.n_2, intersect.bd_2)] += intersect.fraction_2*I;

                if (intersect.bd_3 != -1)
                {
                    mBD_matrix[this_bd_index][bd_index(intersect.n_3, intersect.bd_3)] += intersect.fraction_3*I;
                }


            }else if (phi_i == total_numv-1){

                int_fact = total_theta*(1-sin(mDphi)/mDphi);

                integrate(r, z, 0, PI_CONST, I1, I2, intersect);

                if (need_vr_vz)
                {
                    I = 0;
                }else{
                    I = -dot_product_sign*int_fact*V_factor*I2;
                }

                mBD_matrix[this_bd_index][bd_index(intersect.n_1, intersect.bd_1)] += intersect.fraction_1*I;
                mBD_matrix[this_bd_index][bd_index(intersect.n_2, intersect.bd_2)] += intersect.fraction_2*I;

                if (intersect.bd_3 != -1)
                {
                    mBD_matrix[this_bd_index][bd_index(intersect.n_3, intersect.bd_3)] += intersect.fraction_3*I;
                }

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
                        integrate(r, z, 0, phi, I1, I2, intersect);

                        if (need_vr_vz)
                        {
                            I = dot_product_sign*sin(phi)*int_fact*V_factor*I2;
                        }else{
                            I = dot_product_sign*cos(phi)*int_fact*V_factor*I2;
                        }

                        mBD_matrix[this_bd_index][bd_index(intersect.n_1, intersect.bd_1)] += intersect.fraction_1*I;
                        mBD_matrix[this_bd_index][bd_index(intersect.n_2, intersect.bd_2)] += intersect.fraction_2*I;

                        if (intersect.bd_3 != -1)
                        {
                            mBD_matrix[this_bd_index][bd_index(intersect.n_3, intersect.bd_3)] += intersect.fraction_3*I;
                        }

                    }else if (theta_i == total_numv-1){

                        integrate(r, z, PI_CONST, phi, I1, I2, intersect);

                        if (need_vr_vz)
                        {
                            I = -dot_product_sign*sin(phi)*int_fact*V_factor*I2;
                        }else{
                            I = dot_product_sign*cos(phi)*int_fact*V_factor*I2;
                        }

                        mBD_matrix[this_bd_index][bd_index(intersect.n_1, intersect.bd_1)] += intersect.fraction_1*I;
                        mBD_matrix[this_bd_index][bd_index(intersect.n_2, intersect.bd_2)] += intersect.fraction_2*I;

                        if (intersect.bd_3 != -1)
                        {
                            mBD_matrix[this_bd_index][bd_index(intersect.n_3, intersect.bd_3)] += intersect.fraction_3*I;
                        }

                    }else{

                        if (theta_i == min_theta || min_theta == max_theta)
                        {
                            // if limits are not 0 and pi, then need to alter integration factor
                            int_edge_fact = 1/2;
                        }else{
                            int_edge_fact = 1;
                        }

                        theta = theta_i*mDtheta;

                        integrate(r, z, theta, phi, I1, I2, intersect);

                        if (need_vr_vz)
                        {
                            I = 2*dot_product_sign*int_edge_fact*sin(phi)*cos(theta)*int_fact*V_factor*I2;
                        }else{
                            I = 2*dot_product_sign*int_edge_fact*cos(phi)*int_fact*V_factor*I2;
                        }

                        mBD_matrix[this_bd_index][bd_index(intersect.n_1, intersect.bd_1)] += intersect.fraction_1*I;
                        mBD_matrix[this_bd_index][bd_index(intersect.n_2, intersect.bd_2)] += intersect.fraction_2*I;

                        if (intersect.bd_3 != -1)
                        {
                            mBD_matrix[this_bd_index][bd_index(intersect.n_3, intersect.bd_3)] += intersect.fraction_3*I;
                        }
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
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
void RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
    numv>::solve(bool multithread)
    {
        if (multithread)
        {
            // solve in parallel

            boost::thread_group calcThreads;

            for (unsigned int nr = 1; nr < numr-1; ++nr)
            {
                for (unsigned int nz = 0; nz < numz-1; ++nz)
                {
                    calcThreads.create_thread(boost::bind(&RadialSolver<numr,numz,numr_bd,numz_bd,numv>::solve, this, nr, nz));
                }

            }

            calcThreads.join_all();

        }else{
            // solve serially

            for (unsigned int nr = 1; nr < numr-1; ++nr)
            {
                for (unsigned int nz = 0; nz < numz-1; ++nz)
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
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
void RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
    numv>::solve(
        unsigned int nr,
        unsigned int nz)
    {



        double r = nr*mDr+mMinr;
        double z = mZ[nz];

        // reset to zero

        mSolutionN[nr][nz] = 0;
        mSolutionVr[nr][nz] = 0;
        mSolutionVz[nr][nz] = 0;

        double I1 = 0;
        double I2 = 0;
        bd_intersect intersect;

        double int_fact;
        double N_factor = 2.0/pow(PI_CONST, 3.0/2.0);
        double V_factor = 2.0*mVth/pow(PI_CONST, 3.0/2.0);
        double theta, phi;

        for(unsigned int phi_i = 0; phi_i < total_numv; ++phi_i)
        {
            if (phi_i == 0)
            {
                int_fact = TWOPI_CONST*(1-sin(mDphi)/mDphi);

                integrate(r, z, 0, 0, I1, I2, intersect);

                mSolutionN[nr][nz] += int_fact*N_factor*I1;
                mSolutionVz[nr][nz] += int_fact*V_factor*I2;

            }else if (phi_i == total_numv-1){

                int_fact = TWOPI_CONST*(1-sin(mDphi)/mDphi);

                integrate(r, z, 0, PI_CONST, I1, I2, intersect);

                mSolutionN[nr][nz] += int_fact*N_factor*I1;
                mSolutionVz[nr][nz] += -int_fact*V_factor*I2;

            }else{


                phi = phi_i*mDphi;

                int_fact = 2*mDtheta*sin(phi)*(1-cos(mDphi))/mDphi;

                for (unsigned int theta_i = 0; theta_i < total_numv; ++theta_i)
                {
                    if (theta_i == 0)
                    {
                        integrate(r, z, 0, phi, I1, I2, intersect);

                        mSolutionN[nr][nz] += int_fact*N_factor*I1;
                        mSolutionVr[nr][nz] += sin(phi)*int_fact*V_factor*I2;
                        mSolutionVz[nr][nz] += cos(phi)*int_fact*V_factor*I2;

                    }else if (theta_i == total_numv-1){

                        integrate(r, z, PI_CONST, phi, I1, I2, intersect);

                        mSolutionN[nr][nz] += int_fact*N_factor*I1;
                        mSolutionVr[nr][nz] += -sin(phi)*int_fact*V_factor*I2;
                        mSolutionVz[nr][nz] += cos(phi)*int_fact*V_factor*I2;

                    }else{

                        theta = theta_i*mDtheta;

                        integrate(r, z, theta, phi, I1, I2, intersect);

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
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
void RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
    numv>::integrate(
        double r,
        double z,
        double theta,
        double phi,
        double& nM,
        double& nN,
        bd_intersect& intersect)
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

        double r_intersect = sqrt(r*r + L*L*sinphi*sinphi - 2*r*L*sinphi*costheta);
        double delta_r = r - r_intersect;
        double delta_r2 = delta_r*delta_r;
        double delta_z = -L*cosphi;
        double delta_z2 = delta_z*delta_z;
        double z_intersect = abs(z + delta_z); // is symmetric about +/- z, so just use positive z
        //int sign_z = z - L*cosphi > 0 ? 1 : -1;

        double frac1, frac2, fracA, fracB;
        unsigned int nri1, nri2, nzi1, nzi2;

        if (phi < phi1 || phi > phi2)
        {
            // floor or ceiling

            nri1 = floor((r_intersect - mMinr)/mDr_bd);
            nri2 = nri1+1;
            frac2 = (r_intersect - (mDr_bd*nri1 + mMinr))/mDr_bd;
            frac1 = 1-frac2;


            if (nri1 == 0)
            {
                //includes corner
                fracA = frac1*delta_z2/(delta_r2 + delta_z2);
                fracB = frac1*delta_r2/(delta_r2 + delta_z2);

                n0_wall = fracA*mSolN_BD_top[nri1] + fracB*mSolN_BD_inner[numz_bd-1] + frac2*mSolN_BD_top[nri2];

                intersect.bd_1 = 2;
                intersect.bd_2 = 2;
                intersect.bd_3 = 1;
                intersect.n_1 = nri1;
                intersect.n_2 = nri2;
                intersect.n_3 = numz_bd-1;
                intersect.fraction_1 = fracA;
                intersect.fraction_2 = frac2;
                intersect.fraction_3 = fracB;

            }else if (nri2 == numr_bd-1){

                //includes corner
                fracA = frac2*delta_z2/(delta_r2 + delta_z2);
                fracB = frac2*delta_r2/(delta_r2 + delta_z2);

                n0_wall = frac1*mSolN_BD_top[nri1] + fracA*mSolN_BD_top[nri2] + fracB*mSolN_BD_outer[numz_bd-1];

                intersect.bd_1 = 2;
                intersect.bd_2 = 2;
                intersect.bd_3 = 0;
                intersect.n_1 = nri1;
                intersect.n_2 = nri2;
                intersect.n_3 = numz_bd-1;
                intersect.fraction_1 = frac1;
                intersect.fraction_2 = fracA;
                intersect.fraction_3 = fracB;

            }else{
                // flat part

                n0_wall = frac1*mSolN_BD_top[nri1] + frac2*mSolN_BD_top[nri2];

                intersect.bd_1 = 2;
                intersect.bd_2 = 2;
                intersect.bd_3 = -1;
                intersect.n_1 = nri1;
                intersect.n_2 = nri2;
                intersect.fraction_1 = frac1;
                intersect.fraction_2 = frac2;

            }

        }else{
            // wall
            nzi1 = floor(z_intersect/mDz_bd);
            nzi2 = nzi1+1;
            frac2 = (z_intersect - mDz_bd*nzi1)/mDz_bd;
            frac1 = 1-frac2;

            if (theta <= theta_transition)
            {
                //inner wall


                if (nzi2 == numz_bd-1)
                {
                    fracA = frac2*delta_r2/(delta_r2 + delta_z2);
                    fracB = frac2*delta_z2/(delta_r2 + delta_z2);

                    // includes corner which is different

                    n0_wall = frac1*mSolN_BD_inner[nzi1] + fracA*mSolN_BD_inner[nzi2] + fracB*mSolN_BD_top[0];

                    intersect.bd_1 = 1;
                    intersect.bd_2 = 1;
                    intersect.bd_3 = 2;
                    intersect.n_1 = nzi1;
                    intersect.n_2 = nzi2;
                    intersect.n_3 = 0;
                    intersect.fraction_1 = frac1;
                    intersect.fraction_2 = fracA;
                    intersect.fraction_3 = fracB;

                }else{
                    // flat part
                    n0_wall = frac1*mSolN_BD_inner[nzi1] + frac2*mSolN_BD_inner[nzi2];

                    intersect.bd_1 = 1;
                    intersect.bd_2 = 1;
                    intersect.bd_3 = -1;
                    intersect.n_1 = nzi1;
                    intersect.n_2 = nzi2;
                    intersect.fraction_1 = frac1;
                    intersect.fraction_2 = frac2;
                }

            }else{

                // outer wall

                if (nzi2 == numz_bd-1)
                {
                    fracA = frac2*delta_r2/(delta_r2 + delta_z2);
                    fracB = frac2*delta_z2/(delta_r2 + delta_z2);

                    // includes corner which is different

                    n0_wall = frac1*mSolN_BD_outer[nzi1] + fracA*mSolN_BD_outer[nzi2] + fracB*mSolN_BD_top[numr_bd-1];

                    intersect.bd_1 = 0;
                    intersect.bd_2 = 0;
                    intersect.bd_3 = 2;
                    intersect.n_1 = nzi1;
                    intersect.n_2 = nzi2;
                    intersect.n_3 = numr_bd-1;
                    intersect.fraction_1 = frac1;
                    intersect.fraction_2 = fracA;
                    intersect.fraction_3 = fracB;

                }else{
                    // flat part
                    n0_wall = frac1*mSolN_BD_outer[nzi1] + frac2*mSolN_BD_outer[nzi2];

                    intersect.bd_1 = 0;
                    intersect.bd_2 = 0;
                    intersect.bd_3 = -1;
                    intersect.n_1 = nzi1;
                    intersect.n_2 = nzi2;
                    intersect.fraction_1 = frac1;
                    intersect.fraction_2 = frac2;
                }
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
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
unsigned int RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
    numv>::
    bd_index(
        unsigned int n,
        unsigned int bd)
    {
        if (bd == 0)
        {
            return n;
        }

        if (bd == 1)
        {
            return n + numz_bd;
        }

        if (bd == 2)
        {
            return n + 2*numz_bd;
        }

        throw -1;
    }




template<
    unsigned int numr,
    unsigned int numz,
    unsigned int numr_bd,
    unsigned int numz_bd,
    unsigned int numv>
double RadialSolver<
    numr,
    numz,
    numr_bd,
    numz_bd,
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
