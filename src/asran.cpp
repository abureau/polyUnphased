/* asran.cpp - Wichmall & Hill random number generator with additions

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

#include <cmath>
#include "asran.h"

using namespace std;

// default starting seed
Random::Random() {
    ix = 1;
    iy = 1;
    iz = 1;
}

// specify starting seed
Random::Random(int i, int j, int k) {
    ix = i;
    iy = j;
    iz = k;
}

// specify starting seed
void Random::init(int i, int j, int k) {
    ix = i;
    iy = j;
    iz = k;
}

// random number in [0,1]
double Random::rand() {
    ix = (171 * ix) % 30269;
    iy = (172 * iy) % 30307;
    iz = (170 * iz) % 30323;
    double r  = (double)ix / 30269.0 + (double)iy / 30307.0 + (double)iz / 30323.0;
    return (r - (int)r);
}

// random "allele" with frequencies in f
int Random::ranallele(double *f) {
    double x = rand();
    double y = 0;
    int i = 0;
    for (; y < x; i++) {
        y += f[i];
    }
    return(i - 1);
}

// random deviate from N(0,1)
double Random::normdev() {
    double y = rand();
    double
    t1 = sqrt(log(1.0 / (y * y))),
    t2 = sqrt(log(1.0 / ((1 - y) * (1 - y)))),
    c0 = 2.515517,
    c1 = 0.802853,
    c2 = 0.010328,
    d1 = 1.432788,
    d2 = 0.189269,
    d3 = 0.001308;
    double
    X1 = t1 - (c0 + c1 * t1 + c2 * t1 * t1) / (1 + d1 * t1 + d2 * t1 * t1 + d3 * t1 * t1 * t1),
    X2 = t2 - (c0 + c1 * t2 + c2 * t2 * t2) / (1 + d1 * t2 + d2 * t2 * t2 + d3 * t2 * t2 * t2);
    double x = y < 0.5 ? X1 : -X2;
    return(x);
}

// Bivariate Gaussian deviate with correlation r
void Random::bivnormdev(double *x, double r) {
    x[0] = normdev();
    x[1] = x[0] * r + sqrt(1 - r * r) * normdev();
}
