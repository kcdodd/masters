/*

    Util
    Copyright (C) 2008 K. Carter Dodd (kcdodd@andromedaspace.com)
    http://www.andromedaspace.com

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

#ifndef __Util_TryCatchThrowError__
#define __Util_TryCatchThrowError__

#include <exception>
#include <string>
#include <sstream>

namespace Util {

    class TryCatchThrowError : public std::exception
    {
    private:
        std::string fullDescription;
    public:
        TryCatchThrowError(
            const char* func,
            const char* file) throw()
        {
            std::ostringstream os;

            os << "\tin function: " << func << ", in file: " << file;

            fullDescription = os.str();
        }

        virtual ~TryCatchThrowError() throw()
        {
        }

        virtual const char* what() const throw()
        {
            return fullDescription.c_str();
        }
    };
}
#endif
