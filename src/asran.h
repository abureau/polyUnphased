/* asran.h - Wichmall & Hill random number generator with additions

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

#ifndef __ASRAN__

class Random {
public:
    Random();
    Random(int, int, int);
    void init(int, int, int);
    double rand();
    int ranallele(double *);
    double normdev();
    void bivnormdev(double *, double);
protected:
    int
    ix,
    iy,
    iz;
};

#endif
#define __ASRAN__
