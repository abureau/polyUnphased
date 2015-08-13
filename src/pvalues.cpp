/* pvalues.cpp - asymptotic p-values using Numerical Recipes

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
#include <cmath>

using namespace std;

static double betacf(double a, double b, double x) {
    const int MAXIT = 100;
    const double EPS = 3.0e-7;
    const double FPMIN = 1.0e-30;
    int m, m2;
    double aa, c, d, del, h, qab, qam, qap;

    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;
    c = 1.0;
    d = 1.0 - qab * x / qap;
    if (fabs(d) < FPMIN) {
        d = FPMIN;
    }
    d = 1.0 / d;
    h = d;
    for (m = 1; m <= MAXIT; m++) {
        m2 = 2 * m;
        aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) {
            c = FPMIN;
        }
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) {
            c = FPMIN;
        }
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del - 1.0) < EPS) {
            break;
        }
    }
    return(h);
}

double gammln(double xx) {
    double x, y, tmp, ser;
    static double cof[6] = {76.18009172947146, -86.50532032941677,
                            24.01409824083091, -1.231739572450155,
                            0.1208650973866179e-2, -0.5395239384953e-5
                           };
    int j;
    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++) {
        ser += cof[j] / ++y;
    }
    return -tmp + log(2.5066282746310005 * ser / x);
}

static double betai(double a, double b, double x) {
    double bt;
    if (x == 0.0 || x == 1.0) {
        bt = 0.0;
    } else {
        bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) + b * log(1.0 - x));
    }
    if (x < (a + 1.0) / (a + b + 2.0)) {
        return(bt * betacf(a, b, x) / a);
    } else {
        return(1.0 - bt * betacf(b, a, 1.0 - x) / b);
    }
}

double postudt(double t, double v) {
    double x;
    x = betai(v / 2.0, 0.5, v / (v + t * t));
    return(x > 0 ? x : 1e-45);
}

static void gcf(double &gammcf, double a, double x, double &gln) {
    const int ITMAX = 1000;
    const double EPS = 3.0e-7;
    const double FPMIN = 1.0e-30;
    int i;
    double an, b, c, d, del, h;

    gln = gammln(a);
    b = x + 1.0 - a;
    c = 1.0 / FPMIN;
    d = 1.0 / b;
    h = d;
    for (i = 1; i <= ITMAX; i++) {
        an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (fabs(d) < FPMIN) {
            d = FPMIN;
        }
        c = b + an / c;
        if (fabs(c) < FPMIN) {
            c = FPMIN;
        }
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del - 1.0) < EPS) {
            break;
        }
    }
    gammcf = exp(-x + a * log(x) - (gln)) * h;
}

static void gser(double &gamser, double a, double x, double &gln) {
    const int ITMAX = 1000;
    const double EPS = 3.0e-7;
    const double FPMIN = 1.0e-30;
    int n;
    double sum, del, ap;

    gln = gammln(a);
    ap = a;
    del = sum = 1.0 / a;
    for (n = 1; n <= ITMAX; n++) {
        ++ap;
        del *= x / ap;
        sum += del;
        if (fabs(del) < fabs(sum)*EPS) {
            gamser = sum * exp(-x + a * log(x) - (gln));
            return;
        }
    }
    return;
}

double pochisq(double chisq, double df) {
    double gamser, gammcf, gln;
    double x;

    if (chisq == 0) {
        return(1.0);
    }
    chisq /= 2.0;
    df /= 2.0;
    if (chisq < (df + 1.0)) {
        gser(gamser, df, chisq, gln);
        x = 1.0 - gamser;
    } else {
        gcf(gammcf, df, chisq, gln);
        x = gammcf;
    }
    return(x < 1.0 ? x : 1.0);
}

static double gammq(double a, double x) {
    double gamser, gammcf, gln;

    if (x < 0.0 || a <= 0.0) {
        cout << "Invalid arguments in routine gammq";
    }
    if (x < (a + 1.0)) {
        gser(gamser, a, x, gln);
        return(1.0 - gamser);
    } else {
        gcf(gammcf, a, x, gln);
        return(gammcf);
    }
}

void cntab1(int **nn, int ni, int nj, double &chisq, double &df, double &prob,
            double &cramrv, double &ccc) {
    const double TINY = 1.0e-30;
    int nnj, nni, j, i, minij;
    double sum = 0.0, expctd, sumi[100], sumj[100], temp;


    nni = ni;
    nnj = nj;
    for (i = 0; i < ni; i++) {
        sumi[i] = 0.0;
        for (j = 0; j < nj; j++) {
            sumi[i] += nn[i][j];
            sum += nn[i][j];
        }
        if (sumi[i] == 0.0) {
            --nni;
        }
    }
    for (j = 0; j < nj; j++) {
        sumj[j] = 0.0;
        for (i = 0; i < ni; i++) {
            sumj[j] += nn[i][j];
        }
        if (sumj[j] == 0.0) {
            --nnj;
        }
    }
    df = nni * nnj - nni - nnj + 1;
    chisq = 0.0;

    if ((df) < 0.5) {
        prob = 1.0;
        cramrv = ccc = 0.0;
    } else {
        for (i = 0; i < ni; i++) {
            for (j = 0; j < nj; j++) {
                expctd = sumj[j] * sumi[i] / sum;
                temp = nn[i][j] - expctd;
                chisq += temp * temp / (expctd + TINY);
            }
        }
        prob = gammq(0.5 * df, 0.5 * chisq);
        minij = nni < nnj ? nni - 1 : nnj - 1;
        cramrv = sqrt(chisq / (sum * minij));
        ccc = sqrt(chisq / (chisq + sum));
    }
}
