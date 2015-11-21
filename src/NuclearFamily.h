/* NuclearFamily.h - managing nuclear families

   Copyright (c) 2006 Frank Dudbridge
   MRC Biostatistics Unit
   Robinson Way
   Cambridge CB2 0SR, UK
   frank.dudbridge@mrc-bsu.cam.ac.uk

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
   USA.

*/

#include "stl.h"
#include <string>
#include "LinkageData.h"

#ifndef __NUCLEARFAMILY__

class NuclearFamily {
public:
    NuclearFamily();
    NuclearFamily(string, vector<Subject> &);
    ~NuclearFamily();
    string name;
    Subject father;
    Subject mother;
    vector<Subject> sibs;
    string status;
    void report();
    bool hasfather;
    bool hasmother;
    bool sibship;
protected:
};

class FamilyList: public vector<NuclearFamily> {
public:
    FamilyList();
    ~FamilyList();
    void extract(PedigreeSet &);
};

#endif
#define __NUCLEARFAMILY__
