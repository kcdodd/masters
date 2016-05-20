#ifndef DYNMATRIX_HPP_INCLUDED
#define DYNMATRIX_HPP_INCLUDED

#include <algorithm>

template<class T>
class dynmatrix
{
private:

    unsigned int mRows;
    unsigned int mCols;

    T** mPtr;

public:

    dynmatrix()
    : mRows(0), mCols(0),mPtr(NULL)
    {

    }

    dynmatrix(unsigned rows, unsigned int cols)
        : mRows(rows), mCols(cols)
    {
        mPtr = new T*[mRows];

        for(unsigned int r = 0; r < mRows; ++r)
        {
            mPtr[r] = new T[cols];

            for (unsigned int c = 0; c < mCols; ++c)
            {
                mPtr[r][c] = 0;
            }
        }
    }

    dynmatrix(const dynmatrix<T> &mat2)
    {
        mRows = mat2.mRows;
        mCols = mat2.mCols;

        mPtr = new T*[mRows];

        for(unsigned int r = 0; r < mRows; ++r)
        {
            mPtr[r] = new T[cols];

            for (unsigned int c = 0; c < mCols; ++c)
            {
                mPtr[r][c] = mat2.mPtr[r][c];
            }
        }
    }

    ~dynmatrix()
    {
        for(unsigned int r = 0; r < mRows; ++r)
        {
            delete [] mPtr[r];
        }

        delete [] mPtr;
    }

    void resize(unsigned int rows, unsigned int cols)
    {
        T** tmpPtr = new T*[rows];

        for(unsigned int r = 0; r < rows; ++r)
        {
            tmpPtr[r] = new T[cols];

            for (unsigned int c = 0; c < cols; ++c)
            {
                tmpPtr[r][c] = 0;
            }
        }

        if (mPtr != NULL)
        {

            unsigned int rmin = std::min(rows, mRows);
            unsigned int cmin = std::min(cols, mCols);


            for(unsigned int r = 0; r < rmin; ++r)
            {

                for (unsigned int c = 0; c < cmin; ++c)
                {
                    tmpPtr[r][c] = mPtr[r][c];
                }

                delete[] mPtr[r];
            }

            delete[] mPtr;
        }

        mRows = rows;
        mCols = cols;
        mPtr = tmpPtr;

    }

    unsigned int rows()
    {
        return mRows;
    }

    unsigned int cols()
    {
        return mCols;
    }

    T& operator() (unsigned int r, unsigned int c)
    {
        return mPtr[r][c];
    }



};
#endif // DYNMATRIX_HPP_INCLUDED
