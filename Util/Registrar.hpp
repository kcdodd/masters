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

#ifndef __Util_Registrar__
#define __Util_Registrar__


namespace Util
{
    template<class T, class M>
    struct SetNumOnLoad {
        int classNum;
        SetNumOnLoad () {
            classNum = M::classTypeNum;
            M::classTypeNum++;
        }
    };


    template<class T, class M>
    class Registrar {
    public:

        static SetNumOnLoad<T, M> numRef;

        static int getNum () {
            return numRef.classNum;
        }
    };

    template<class T, class M>
    SetNumOnLoad<T, M> Registrar<T, M>::numRef;
}

#endif
