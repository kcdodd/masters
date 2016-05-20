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

#ifndef __Util_GeneralError__
#define __Util_GeneralError__

#include <exception>
#include <string>
#include <sstream>

namespace Util
{

    class GeneralError : public std::exception
    {
    private:
        std::string fullDescription;
    public:
        GeneralError(
            const std::string& description,
            const char* file,
            long line) throw()
        {
            std::ostringstream os;

            os << description << ", in file: " << file << ", at line: " << line;

            fullDescription = os.str();
        }

        virtual ~GeneralError() throw()
        {
        }

        virtual const char* what() const throw()
        {
            return fullDescription.c_str();
        }
    };

}
#endif
