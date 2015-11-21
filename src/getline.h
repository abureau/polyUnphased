/* getline.h - read lines of text

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
#include <string>

using namespace std;

#ifndef __GETLINE__

// check whether string has any non-whitespace characters
static bool whitespace(string &s) {
    for (int i = 0; i < s.length(); i++)
        if (s[i] != ' ' && s[i] != '\t' && s[i] != '\r' && s[i] != '\f' && s[i] != '\n') {
            return(false);
        }
    return(true);
}

// read a nonempty line from an istream
static string getline(ifstream &istr) {
    const int buflen = 256;
    char buffer[buflen];
    string s = "";
    while (whitespace(s) && !istr.eof()) {
        istr.getline(buffer, buflen);
        s = buffer;
        while (istr.fail() && !istr.eof()) {
            istr.clear();
            istr.getline(buffer, buflen);
            s += buffer;
        }
    }
    return(s);
}

#endif
#define __GETLINE__
