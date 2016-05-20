#ifndef DYNARRAY_HPP_INCLUDED
#define DYNARRAY_HPP_INCLUDED


#include <algorithm>
#include <cmath>
#include <vector>

#include <string>
#include <iostream>
#include <fstream>


template<class T>
class dynarray
{
private:

    unsigned int mLength;

    T* mArrPtr;

    T mXmin;
    T mXmax;
    T mdx;

public:

    dynarray()
    : mLength(0), mArrPtr(NULL), mXmin(0), mXmax(0), mdx(0)
    {
    }

    /**
        constructor.
    */
    dynarray(unsigned int length, T xmin, T xmax)
        : mLength(length), mXmin(xmin), mXmax(xmax), mdx((xmax-xmin)/((T)(length-1)))
    {
        mArrPtr = new T[mLength];

        for(unsigned int i=0; i < mLength; ++i)
        {
            mArrPtr[i] = 0;
        }

    }

    /**
        copy constructor
    */
    dynarray(const dynarray<T>& arr2)
    {
        mLength = arr2.mLength;
        mXmin = arr2.mXmin;
        mXmax = arr2.mXmax;
        mdx = arr2.mdx;

        mArrPtr = new T[mLength];

        for(unsigned int i=0; i < mLength; ++i)
        {
            mArrPtr[i] = arr2.mArrPtr[i];
        }

    }

    ~dynarray()
    {
        delete [] mArrPtr;
    }

    void set_xminmax(T xmin, T xmax)
    {
        if (xmax <= xmin)
        {
            throw 1;
        }

        mXmin = xmin;
        mXmax = xmax;

        mdx = (mXmax-mXmin)/((T)(mLength-1));
    }

    T xmin () const
    {
        return mXmin;
    }

    T xmax() const
    {
        return mXmax;
    }

    T dx () const
    {
        return mdx;
    }

    T x(unsigned int i) const
    {
        return mXmin + ((T)i)*mdx;
    }

    void resize(unsigned int length)
    {
        T* tmpArrPtr = new T[length];


        for(unsigned int i=0; i < length; ++i)
        {
            tmpArrPtr[i] = 0;
        }

        if (mArrPtr != NULL)
        {

            unsigned int lmin = std::min(length, mLength);

            for(unsigned int i=0; i < lmin; ++i)
            {
                tmpArrPtr[i] = mArrPtr[i];
            }

            delete[] mArrPtr;
        }

        mLength = length;
        mArrPtr = tmpArrPtr;

    }

    unsigned int length() const
    {
        return mLength;
    }


    T& operator() (unsigned int i)
    {
        return mArrPtr[i];
    }

    T operator() (unsigned int i) const
    {
        return mArrPtr[i];
    }

    void derivative(dynarray<T>& out) const
    {
        if (out.length() != length())
        {
            out.resize(length());
        }

        if (mdx == 0)
        {
            throw 1;
        }


        out.mXmin = mXmin;
        out.mXmax = mXmax;
        out.mdx = mdx;

        // use 5 point method

        for (unsigned int i = 2; i < (mLength-2); ++i)
        {
            out.mArrPtr[i] = (-mArrPtr[i+2] + 8.0*mArrPtr[i+1] - 8.0*mArrPtr[i-1] + mArrPtr[i-2])/(12.0*mdx);
        }

        out.mArrPtr[0] = 0;
        out.mArrPtr[1] = (-mArrPtr[3] + 8.0*mArrPtr[2] - 8.0*mArrPtr[0] + mArrPtr[1])/(12.0*mdx);


        out.mArrPtr[mLength-2] = (3.0*mArrPtr[mLength-1] + 10.0*mArrPtr[mLength-2] -18.0*mArrPtr[mLength-3] + 6.0*mArrPtr[mLength-4] - mArrPtr[mLength-5])/(12.0*mdx);
        out.mArrPtr[mLength-1] = (25.0*mArrPtr[mLength-1] - 48.0*mArrPtr[mLength-2] + 36.0*mArrPtr[mLength-3] - 16.0*mArrPtr[mLength-4] + 3.0*mArrPtr[mLength-5])/(12.0*mdx);


    }

    // computes second derivative
    void derivative_2(dynarray<T>& out) const
    {
        if (out.mLength != mLength)
        {
            out.resize(mLength);
        }

        if (mdx == 0)
        {
            throw 1;
        }

        out.mXmin = mXmin;
        out.mXmax = mXmax;
        out.mdx = mdx;

        for(unsigned int i = 2; i < mLength-2; ++i)
        {
            out.mArrPtr[i] = (-mArrPtr[i+2] + 16.0*mArrPtr[i+1] - 30.0*mArrPtr[i] + 16.0*mArrPtr[i-1] - mArrPtr[i-2])/(12.0*mdx*mdx);
        }

        out.mArrPtr[1] = (-mArrPtr[3] + 16.0*mArrPtr[2] - 30.0*mArrPtr[1] + 16.0*mArrPtr[0] - mArrPtr[1])/(12.0*mdx*mdx);
        out.mArrPtr[0] = (-mArrPtr[2] + 16.0*mArrPtr[1] - 30.0*mArrPtr[0] + 16.0*mArrPtr[1] - mArrPtr[2])/(12.0*mdx*mdx);

        out.mArrPtr[mLength-2] = (11.0*mArrPtr[mLength-1] - 20.0*mArrPtr[mLength-2] + 6.0*mArrPtr[mLength-3] + 4.0*mArrPtr[mLength-4] - mArrPtr[mLength-5])/(12.0*mdx*mdx);
        out.mArrPtr[mLength-1] = (35.0*mArrPtr[mLength-1] - 104.0*mArrPtr[mLength-2] + 114.0*mArrPtr[mLength-3] - 56.0*mArrPtr[mLength-4] + 11.0*mArrPtr[mLength-5])/(12.0*mdx*mdx);
    }

    void integral(dynarray<T>& out) const
    {
        if (out.length() != mLength)
        {
            out.resize(mLength);
        }

        if (mdx == 0)
        {
            throw 1;
        }

        out.mXmin = mXmin;
        out.mXmax = mXmax;
        out.mdx = mdx;

        out.mArrPtr[0] = 0;

        for (unsigned int i = 1; i < mLength; i += 2)
        {
            // trap. rule for mid point
            out.mArrPtr[i] = out.mArrPtr[i-1] + mdx*0.5*(mArrPtr[i-1] + mArrPtr[i]);

            if (i+1 < mLength)
            {
                // simpsons rule for endpoint
                out.mArrPtr[i+1] = out.mArrPtr[i-1] + 2.0*mdx*(mArrPtr[i-1] + 4.0*mArrPtr[i] + mArrPtr[i+1])/6.0;
            }
        }
    }

    T integral_total() const
    {
        if (mdx == 0)
        {
            throw 1;
        }

        T integral = 0;

        for (unsigned int i = 1; i < mLength; i += 2)
        {

            if (i+1 < mLength)
            {
                // simpsons rule
                integral += 2.0*mdx*(mArrPtr[i-1] + 4.0*mArrPtr[i] + mArrPtr[i+1])/6.0;
            }
        }

        if (mLength%2 == 0)
        {
            // trap. rule for remaining point
            integral += mdx*0.5*(mArrPtr[mLength-2] + mArrPtr[mLength-1]);
        }

        return integral;

    }

    void smooth(dynarray<T> &out, unsigned int nk) const
    {
        if (out.length() != mLength)
        {
            out.resize(mLength);
        }

        out.mXmin = mXmin;
        out.mXmax = mXmax;
        out.mdx = mdx;

        // construct kernel

        std::vector<T> kernel(nk);
        T dx = 4.0/(nk-1);
        T x0 = -2.0;
        T norm = 0.0;

        for(unsigned int i=0; i < nk; ++i)
        {
            T x = x0 + i*dx;
            kernel[i] = exp(-x*x);
            norm += kernel[i];
        }

        for(unsigned int i=0; i < nk; ++i)
        {
            kernel[i] /= norm;
        }

        const unsigned int nk_m = (nk-1)/2;
        unsigned int ij = 0;

        for (unsigned int i = 0; i < mLength; ++i)
        {
            out.mArrPtr[i] = 0;

            for (unsigned int j = 0; j < nk; ++j)
            {
                ij = i-nk_m+j-1;

                if (ij < 0)
                {
                    ij = abs(ij)+1;
                }

                if (ij > mLength-1)
                {
                    ij = 2*(mLength-1)-ij;
                }

                out.mArrPtr[i] += mArrPtr[ij]*kernel[j];

            }

        }

    }

    T interp(T x) const
    {
        if (x > mXmax || x < mXmin)
        {
            throw -1;
        }

        unsigned int nx = floor((x - mXmin)/mdx);

        if (nx == mLength-1)
        {
            // for last point
            return mArrPtr[mLength-1];
        }

        //return mArrPtr[nx];
        //return (mArrPtr[nx+1] - mArrPtr[nx])*(x - (((T)nx)*mdx+mXmin))/mdx + mArrPtr[nx];

        if (nx == 0 || nx == mLength-2)
        {
            // linear interp for end segments
            return (mArrPtr[nx+1] - mArrPtr[nx])*(x - (((T)nx)*mdx+mXmin))/mdx + mArrPtr[nx];
        }else{

            // cubic interp for all middle segments

            T alpha = (x-(((T)nx)*mdx+mXmin))/mdx;


            T a3 = (3.0*(mArrPtr[nx] - mArrPtr[nx+1])+mArrPtr[nx+2]-mArrPtr[nx-1])/6.0;
            T a2 = (mArrPtr[nx+1]+mArrPtr[nx-1]-2.0*mArrPtr[nx])/2.0;
            T a1 = (6.0*mArrPtr[nx+1] - 3.0*mArrPtr[nx] - mArrPtr[nx+2] - 2.0*mArrPtr[nx-1])/6.0;

            return a3*alpha*alpha*alpha + a2*alpha*alpha + a1*alpha + mArrPtr[nx];


        }
    }

    // fills this array using another array using linear interpolation
    // assumes all values can be found in other array.
    void interp(dynarray<T>& in)
    {
        for (unsigned int i = 0; i < mLength; ++i)
        {
            mArrPtr[i] = in.interp(mXmin + mdx*((T)i) );
        }

    }

    void load_text_file(const std::string& filename)
    {
        std::ifstream inputfile(filename.c_str());


        if (inputfile.is_open())
        {
            unsigned int input_length = 0;

            inputfile >> input_length;

            resize(input_length);

            T input_min;
            T input_max;

            inputfile >> input_min;
            inputfile >> input_max;

            set_xminmax(input_min, input_max);

            for (unsigned int i=0; i < mLength; ++i)
            {
                inputfile >> mArrPtr[i];

            }

            inputfile.close();
        }


    }

};

#endif // ARRAY_HPP_INCLUDED
