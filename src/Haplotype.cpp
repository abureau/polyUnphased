/* Haplotype.cpp - haplotype class

   Copyright (c) 2006 Frank Dudbridge
   MRC Biostatistics Unit
   Robinson Way
   Cambridge CB2 2SR, UK
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

#include <iostream>
#include <strstream>
#include <cmath>
#include "Haplotype.h"

// convert an integer into a string
// gets round buggy things with ostrstreams
const string tostring(int i, bool ACGT) {
    string s;
    if (ACGT) {
        if (i == 1) {
            s = "A";
        }
        if (i == 2) {
            s = "C";
        }
        if (i == 3) {
            s = "G";
        }
        if (i == 4) {
            s = "T";
        }
    } else
        while (i) {
            char c = i % 10 + '0';
            s = c + s;
            i /= 10;
        }
    return(s);
}

// formatted haplotype data
const string Haplotype::str(int ncondition, bool genotype, bool ACGT) {
    string s = "";
    int ix = 0;
    for (int i = 0; i < ncondition; i++) {
        s += tostring((*this)[ix++], ACGT);
        s += "/" + tostring((*this)[ix++], ACGT);
        s += "-";
    }
    while (ix < size()) {
        s += tostring((*this)[ix++], ACGT);
        if (genotype) {
            s += "/" + tostring((*this)[ix++], ACGT);
        }
        if (ix < size()) {
            s += "-";
        }
    }
    return(s);
}

const string Haplotype::str(bool ACGT) {
    string s = "";
    int ix = 0;
    for (int i = 0; i < size(); i++) {
        s += tostring((*this)[ix++], ACGT);
        if (i < size() - 1) {
            s += "-";
        }
    }
    return(s);
}

const char *Haplotype::c_str(bool ACGT) {
    string s = this->str(ACGT);
    return(s.c_str());
}

// output a haplotype

// ostream & operator << (ostream & ostr, Haplotype &haplotype) {
//   ostr << haplotype.str();
//   return(ostr);
// }

// default constructor
Haplotype::Haplotype(): vector<int> () {}

// constructor with size {
Haplotype::Haplotype(int n): vector<int> (n, 0) {}

// destructor
Haplotype::~Haplotype() {}

// copy a vector into a Haplotype
void Haplotype::operator = (vector<int> &y) {
    resize(y.size(), 0);
    assign(y.begin(), y.end());
}

// concatenate two Haplotypes
Haplotype Haplotype::operator + (const Haplotype &haplotype) {
    Haplotype result = *this;
    for (int i = 0; i < haplotype.size(); i++) {
        result.push_back(haplotype[i]);
    }
    return(result);
}

// check whether haplotype has any zero elements
bool Haplotype::nozero() {
    for (int i = 0; i < size(); i++)
        if ((*this)[i] == 0) {
            return(false);
        }
    return(true);
}

// return the index of the first element of array matching this haplotype
int Haplotype::index(vector<Haplotype> &array) {
    return(find(array.begin(), array.end(), (*this)) - array.begin());
}

// enter this haplotype into a vector of haplotypes, if not already there
void Haplotype::enter(vector<Haplotype> &array) {
    if (find(array.begin(), array.end(), (*this)) == array.end()) {
        array.push_back(*this);
    }
}

// make a genotype from this haplotype and the input haplotype
Haplotype Haplotype::formgeno(Haplotype &hap, int ncondition, bool genotype) {
    if (hap.size() < size()) {
        hap.resize(size(), 0);
    }
    Haplotype geno(genotype ? (size() * 2) : (size() + ncondition));
    int ix = 0;
    for (int i = 0; i < ncondition; i++) {
        geno[ix++] = min((*this)[i], hap[i]);
        geno[ix++] = max((*this)[i], hap[i]);
    }
    for (int i = ncondition; i < size(); i++) {
        if (genotype) {
            geno[ix++] = min((*this)[i], hap[i]);
            geno[ix++] = max((*this)[i], hap[i]);
        } else {
            geno[ix++] = (*this)[i];
        }
    }
    return(geno);
}

