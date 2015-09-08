/* UnphasedAnalysisLikelihood.cpp - Association analysis in trios and unrelateds

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

// system includes
#include <iostream>
#include <cmath>

// local includes
#include "stl.h"
#include "UnphasedOptions.h"
#include "UnphasedAnalysis.h"
#include "asran.h"

static Random ran;

//  sizeOfGradient

// total number of parameters in the model
int UnphasedAnalysis::sizeOfGradient(UnphasedOptions &options) {

    int size = 0;
    int nhap = genoCode.size();
    int betasize;
    if (typeOfPhenotype == "polytomous") betasize = nhap*(K-1);
    else betasize = nhap;

    // beta
    else size += betasize;

    // freq
    size += nhap;

    // betaparent
    if (haveBetaParent(options)) {
        size += betasize + haveBetaParent0(options);
    }

    // alpha
    if (haveAlpha(options)) {
        size += nhap;
        if (haveAlpha0(options)) {
            size++;    // intercept
        }
    }

    // covariates
    for (int j = 0; j < betaCovariate.size(); j++)
        for (int k = 0; k < betaCovariate[j].size(); k++) {
            size += nhap; // betaCovariate
            if (!confounder[j] && haveFamilies && !options.hhrr) {
                size += nhap;    // betaparentCovariate
            }
            if (typeOfPhenotype == "quant") {
                size++; // intercept
                if (haveFamilies && !options.hhrr) {
                    size++;    // betaparent intercept
                }
            }
        }

    return(size);

}



//  freeParameters

// number of free variables in the likelihood
int UnphasedAnalysis::freeParameters(UnphasedOptions &options) {

    // groups of betas
    int ndim = 0;
    for (int i = 0; i < group.size(); i++) {
        int thisndim = 0;
        for (int j = 0; j < genoCode.size(); j++) if (!zero[j]) {
                thisndim = max(thisndim, group[i][j]);
            }
        if (typeOfPhenotype == "polytomous") thisndim *= K-1;
        ndim += thisndim;
    }

    // nuisance parameters
    int refIndex = reference.index(genoCode);
    int nhap = genoCode.size();
    // freq
    for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
            ndim++;
        }
    // betaparent
    if (haveBetaParent(options)) {
        for (int i = 0; i < nhap; i++)
            if (!zero[i] && i != refIndex) {
                if (typeOfPhenotype == "polytomous") ndim += K-1;
                else ndim++;
            }
        if (haveBetaParent0(options)) {
            ndim++;    // intercept
        }
    }
    // alpha
    if (haveAlpha(options)) {
        for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                ndim++;
            }
        if (haveAlpha0(options)) {
            ndim++;    // intercept
        }
    }
    // covariates
    for (int j = 0; j < betaCovariate.size(); j++)
        for (int k = 0; k < betaCovariate[j].size(); k++) {
            for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                    ndim++; // betaCovariate
                    if (haveFamilies && !confounder[j] && !options.hhrr) {
                        ndim++;    // betaparentCovariate
                    }
                }
            if (typeOfPhenotype == "quant") {
                ndim++; // intercept
                if (haveFamilies && !options.hhrr) {
                    if (typeOfPhenotype == "polytomous") ndim += K-1;
                    else ndim++;    // betaparent intercept
                }
            }
        }

    return(ndim);

}



//  gradientFreeParameters

// form the gradient of the free parameters
void UnphasedAnalysis::gradientFreeParameters(vector<double> &gradient,
        vector<double> &g,
        UnphasedOptions &options) {

    int nhap = genoCode.size();
    int refIndex = reference.index(genoCode);
    for (int i = 0; i < g.size(); i++) {
        g[i] = 0;
    }

    // beta
    vector<int> ndim(group.size(), 0);
    for (int i = 1; i < group.size(); i++) {
        ndim[i] = 0;
        for (int j = 0; j < nhap; j++) if (!zero[j]) {
                ndim[i] = max(ndim[i], group[i-1][j]);
            }
        ndim[i] += ndim[i-1];
    }
    for (int i = 0; i < group.size(); i++) {
        for (int j = 0; j < nhap; j++) if (!zero[j] && !(i && rare[j])) {
                if (options.model != "commonmain" ||
                        i == 0 && options.condition.size() > 0) {
                    if (group[i][j] < group[i][refIndex]) {
                        g[group[i][j] + ndim[i]] += gradient[j];
                        if (typeOfPhenotype == "polytomous")
                        	{
                        	for (int k = 1; k < K-1; k++)
                        		g[group[i][j] + ndim[i] + k*nhap] += gradient[j + k*nhap];
                        	}
                    }
                    if (group[i][j] > group[i][refIndex]) {
                        g[group[i][j] + ndim[i] - 1] += gradient[j];
                        if (typeOfPhenotype == "polytomous")
                        	{
                        	for (int k = 1; k < K-1; k++)
                        		g[group[i][j] + ndim[i] -1 + k*nhap] += gradient[j + k*nhap];
                        	}
                    }
                } else {
                    if (group[i][j] != group[i][refIndex]) {
                        g[ndim[options.condition.size()>0]] += gradient[j];
                        if (typeOfPhenotype == "polytomous")
                        	{
                        	for (int k = 1; k < K-1; k++)
                        		g[ndim[options.condition.size()>0] + k*nhap] += gradient[j + k*nhap];
                        	}
                    }
                }
            }
    }

    int ix = g.size() - 1;
    int iy;
    if (typeOfPhenotype == "polytomous") iy = nhap*(K-1);
    else iy = nhap;
    // freq
    for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
            g[ix--] = gradient[i+iy];
        }
    iy += nhap;
    // betaparent
    if (haveBetaParent(options)) {
        for (int i = 0; i < nhap; i++)
            if (!zero[i] && i != refIndex) {
                        if (typeOfPhenotype == "polytomous") {
                        	for (int k = 0; k < K-1; k++)
                        		g[ix - (K-2-k)*nhap] += gradient[j + k*nhap];
                        	}
                else g[ix] = gradient[i+iy];
                ix--;
            }
	    if (typeOfPhenotype == "polytomous") {
	    	iy += nhap*(K-1);
	    	ix -= nhap*(K-2);
	    	}
    	else iy += nhap;
        if (haveBetaParent0(options)) {
            g[ix--] = gradient[iy++];
        }
    }
    // alpha
    if (haveAlpha(options)) {
        for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                g[ix--] = gradient[i+iy];
            }
        iy += nhap;
        if (haveAlpha0(options)) {
            g[ix--] = gradient[iy++];
        }
    }
    // covariates
    for (int j = 0; j < betaCovariate.size(); j++)
        for (int k = 0; k < betaCovariate[j].size(); k++) {
            for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                    g[ix--] = gradient[i+iy];
                }
            iy += nhap;
            if (haveFamilies && !confounder[j] && !options.hhrr) {
                for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                        g[ix--] = gradient[i+iy];
                    }
                iy += nhap;
            }
            if (typeOfPhenotype == "quant") {
                g[ix--] = gradient[iy++];
                if (haveFamilies && !options.hhrr) {
                    g[ix--] = gradient[iy++];
                }
            }
        }
}



//  NelderMead

typedef vector<double> dvector;

/* Nelder-Mead maximisation of log-likelihood
   Hacked from code in Numerical Recipes in C (Press et al.) */

inline void get_psum(vector<dvector> &p, dvector &psum) {
    for (int j = 0; j < p.size() - 1; j++) {
        double sum = 0;
        for (int i = 0; i < p.size(); i++) {
            sum += p[i][j];
        }
        psum[j] = sum;
    }
}

double UnphasedAnalysis::amotry(vector<dvector> &p, dvector &y,
                                dvector &psum,
                                int ihi, double fac, UnphasedOptions &options,
                                int &neval) {
    int ndim = p.size() - 1;
    dvector ptry(ndim, 0);
    dvector g(ndim, 0);
    double fac1 = (1.0 - fac) / ndim;
    double fac2 = fac1 - fac;
    for (int j = 0; j < ndim; j++) {
        ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;
    }
    double ytry = evaluate(ptry, g, options, neval);
    if (ytry < y[ihi]) {
        y[ihi] = ytry;
        for (int j = 0; j < ndim; j++) {
            psum[j] += ptry[j] - p[ihi][j];
            p[ihi][j] = ptry[j];
        }
    }
    return ytry;
}

double UnphasedAnalysis::amoeba(UnphasedOptions &options, bool randomise, int &neval) {

    int ndim = freeParameters(options);
    dvector g(ndim, 0);
    vector<dvector> p;
    p.resize(ndim + 1);
    dvector y(ndim + 1, 0);
    for (int i = 0; i < ndim + 1; i++) {
        p[i].resize(ndim, 0);
        if (randomise) for (int j = 0; j < ndim; j++) {
                p[i][j] = ran.rand();
            }
        if (i) {
            p[i][i-1] += 1.0;
        }
        y[i] = evaluate(p[i], g, options, neval);
    }

    const int NMAX = 50000;
    vector<double> psum(ndim, 0);
    int nfunk = 0;
    int mpts = ndim + 1;
    get_psum(p, psum);
    while (true) {
        int ilo = 0;
        int ihi = ndim ? (y[0] > y[1] ? 0 : 1) : 0;
        int inhi = ndim ? (y[0] > y[1] ? 1 : 0) : 0;
        for (int i = 0; i < mpts; i++) {
            if (y[i] <= y[ilo]) {
                ilo = i;
            }
            if (y[i] > y[ihi]) {
                inhi = ihi;
                ihi = i;
            } else if (y[i] > y[inhi] && i != ihi) {
                inhi = i;
            }
        }

        double rtol = 2.0 * fabs(y[ihi] - y[ilo]) / (fabs(y[ihi]) + fabs(y[ilo]));
        if (fabs(y[ihi]) + fabs(y[ilo]) == 0) {
            rtol = 0;
        }
        if (rtol < options.epsilon) {
            swap(y[0], y[ilo]);
            for (int i = 0; i < ndim; i++) {
                swap(p[0][i], p[ilo][i]);
            }
            break;
        }
        if (nfunk >= NMAX) {
            *outStream << "WARNING: did not maximise likelihood after " << NMAX << " evaluations" << endl;
            return(-y[0]);
        }
        nfunk += 2;
        double ytry = amotry(p, y, psum, ihi, -1.0, options, neval);
        if (ytry <= y[ilo]) {
            ytry = amotry(p, y, psum, ihi, 2.0, options, neval);
        } else if (ytry >= y[inhi]) {
            double ysave = y[ihi];
            ytry = amotry(p, y, psum, ihi, 0.5, options, neval);
            if (ytry >= ysave) {
                for (int i = 0; i < mpts; i++) {
                    if (i != ilo) {
                        for (int j = 0; j < ndim; j++) {
                            p[i][j] = psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
                        }
                        y[i] = evaluate(psum, g, options, neval);
                    }
                }
                nfunk += ndim;
                get_psum(p, psum);
            }
        } else {
            nfunk--;
        }
    }
    return(-y[0]);
}



//  DPF minimisation

/* Routines for Davidon-Powell-Fletcher maximisation of log-likelihood
   Based on code in Numerical Recipes in C (Press et al.) */

//  lnsrch

double UnphasedAnalysis::lnsrch(vector<double> &xold, double fold,
                                vector<double> &g,
                                vector<double> &p,
                                vector<double> &x,
                                double stpmax, int &check,
                                UnphasedOptions &options, int &neval) {

    const double ALF = 1.0e-4;
    const double TOLX = 1.0e-7;

    check = 0;
    double sum = 0;
    int n = p.size();
    for (int i = 0; i < n; i++) {
        sum += p[i] * p[i];
    }
    sum = sqrt(sum);
    if (sum > stpmax)
        for (int i = 0; i < n; i++) {
            p[i] *= stpmax / sum;
        }
    double slope = 0;
    for (int i = 0; i < n; i++) {
        slope += g[i] * p[i];
    }

    if (slope >= 0.0) {
        //      *outStream << "Roundoff problem in lnsearch 1 " << slope << endl;
        //      exit(-1);
        check = -1;
        //      cout << "roundoff 1" << " " << slope << endl;
        return(fold);
    }

    double test = 0.0;
    for (int i = 0; i < n; i++) {
        double temp = fabs(p[i]) / max(fabs(xold[i]), 1.0);
        if (temp > test) {
            test = temp;
        }
    }
    double alamin = TOLX / test;
    double alam = 1.0;

    // my modification to save time, only it doesn't always work
    //  for (int its=0;its<100;its++) {

    for (;;) { // original code
        double alam2;
        double tmplam;
        double fold2;
        double f, f2;
        for (int i = 0; i < n; i++) {
            x[i] = xold[i] + alam * p[i];
        }
        f = evaluate(x, g, options, neval);

        if (alam < alamin) {
            for (int i = 0; i < n; i++) {
                x[i] = xold[i];
            }
            check = 1;
            return(f);
        } else if (f <= fold + ALF * alam * slope) {
            return(f);
        } else {
            if (alam == 1.0) {
                tmplam = -slope / (2.0 * (f - fold - slope));
            } else {
                double rhs1 = f - fold - alam * slope;
                double rhs2 = f2 - fold2 - alam2 * slope;
                double a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
                double b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
                if (a == 0.0) {
                    tmplam = -slope / (2.0 * b);
                } else {
                    double disc = b * b - 3.0 * a * slope;
                    if (disc < 0.0) {
                        //        *outStream << "Roundoff problem in lnsrch 2" << disc<<endl;
                        //       exit(-1);
                        check = -1;
                        //       cout << "roundoff 2" << endl;
                        return(fold);
                    }
                    if (disc < 0.0) {
                        tmplam = 0.5 * alam;
                    } else if (b <= 0.0) {
                        tmplam = (-b + sqrt(disc)) / (3.0 + a);
                    } else {
                        tmplam = -slope / (b + sqrt(disc));
                    }
                }
                if (tmplam > 0.5 * alam) {
                    tmplam = 0.5 * alam;
                }
            }
        }
        alam2 = alam;
        f2 = f;
        fold2 = fold;
        alam = max(tmplam, 0.1 * alam);
    }
}



double UnphasedAnalysis::mydfpmin(UnphasedOptions &options, bool randomise, int &neval) {
    const int ITMAX = 1000;
    const double EPS = 3.0e-8;
    const double TOLX = 4 * EPS;
    const double STPMX = 100.0; // was 100 - FD 16/3/05

    int ndim = freeParameters(options);

    vector<double> p(ndim, 0);
    if (randomise) for (int i = 0; i < ndim; i++) {
            p[i] = ran.rand();
        }
    vector<double> dg(ndim, 0);
    vector<double> g(ndim, 0);
    vector<double> hdg(ndim, 0);
    vector<vector<double> > hessin;
    vector<double> pnew(ndim, 0);
    vector<double> xi(ndim, 0);
    double fp = evaluate(p, g, options, neval);
    double sum = 0;
    double sumg = 0;
    hessin.resize(ndim);
    for (int i = 0; i < ndim; i++) {
        hessin[i].resize(ndim);
        for (int j = 0; j < ndim; j++) {
            hessin[i][j] = 0.0;
        }
        hessin[i][i] = 1.0;
        xi[i] = -g[i];
        sum += p[i] * p[i];
        sumg += g[i] * g[i];
    }
    if (sumg == 0) {
        return(fp);
    }

    double fret = 0; // return value
    double stpmax = STPMX * max(sqrt(sum), (double)ndim);

    for (int its = 0; its < ITMAX; its++) {
        dg = g;
        int check = 0;
        fret = lnsrch(p, fp, g, xi, pnew, stpmax, check, options, neval);

        // use NelderMead if there was a roundoff problem
        if (check == -1) {
            //neval=-neval;
            return(-fret);
        }

        if (fabs(fp - fret) < options.epsilon) {
            goto returnPoint;
        }
        fp = fret;
        for (int i = 0; i < ndim; i++) {
            xi[i] = pnew[i] - p[i];
        }
        p = pnew;
        double test = 0.0;
        for (int i = 0; i < ndim; i++) {
            double temp = fabs(xi[i]) / max(fabs(p[i]), 1.0);
            if (temp > test) {
                test = temp;
            }
        }
        if (test < TOLX) {
            goto returnPoint;
        }

        evaluate(p, g, options, neval);
        test = 0.0;
        double den = max(fret, 1.0);
        for (int i = 0; i < ndim; i++) {
            double temp = fabs(g[i]) * max(fabs(p[i]), 1.0) / den;
            if (temp > test) {
                test = temp;
            }
        }
        if (test < options.epsilon) {
            goto returnPoint;
        }
        for (int i = 0; i < ndim; i++) {
            dg[i] = g[i] - dg[i];
        }
        for (int i = 0; i < ndim; i++) {
            hdg[i] = 0.0;
            for (int j = 0; j < ndim; j++) {
                hdg[i] += hessin[i][j] * dg[j];
            }
        }
        double fac = 0, fad = 0, fae = 0, sumdg = 0, sumxi = 0.0;
        for (int i = 0; i < ndim; i++) {
            fac += dg[i] * xi[i];
            fae += dg[i] * hdg[i];
            sumdg += dg[i] * dg[i];
            sumxi += xi[i] * xi[i];
        }
        if (fac > sqrt(EPS * sumdg * sumxi)) {
            fac = 1.0 / fac;
            fad = 1.0 / fae;
            for (int i = 0; i < ndim; i++) {
                dg[i] = fac * xi[i] - fad * hdg[i];
            }
            for (int i = 0; i < ndim; i++) {
                for (int j = i; j < ndim; j++) {
                    hessin[i][j] += fac * xi[i] * xi[j]
                                    - fad * hdg[i] * hdg[j] + fae * dg[i] * dg[j];
                    hessin[j][i] = hessin[i][j];
                }
            }
        }
        for (int i = 0; i < ndim; i++) {
            xi[i] = 0.0;
            for (int j = 0; j < ndim; j++) {
                xi[i] -= hessin[i][j] * g[j];
            }
        }
    }

    // failed to converge in these iterations, return anyway
    //neval=-neval;
    return(-fret);

returnPoint:
    if (0) { // output numerical gradients for debugging
        for (int i = 0; i < g.size(); i++) {
            cout << g[i] << " ";
        }
        cout << endl;
        double ss = fret;
        for (int i = 0; i < p.size(); i++) {
            p[i] += 0.000001;
            double tt = evaluate(p, g, options, neval);
            cout << "numerical gradient " << (tt - ss) / 0.000001 << endl;
            p[i] -= 0.000001;
        }
    }
    return(-fret);
}



//  evaluate

/* evaluate the log-likelihood for the free parameters in y
   returns score vector in gradient
   Instead of using summary counts, we now loop through all the subjects
   so we can include individual covariate measurements
*/

double UnphasedAnalysis::evaluate(vector<double> &y, vector<double> &g,
                                  UnphasedOptions &options, int &neval) {
    int nhap = genoCode.size();
    int refIndex = reference.index(genoCode);

    vector<int> ndim(group.size(), 0);
    for (int i = 1; i < group.size(); i++) {
        ndim[i] = 0;
        for (int j = 0; j < nhap; j++) if (!zero[j]) {
                ndim[i] = max(ndim[i], group[i-1][j]);
            }
        ndim[i] += ndim[i-1];
    }

    //    for(int i=0;i<ndim.size();i++) cout << ndim[i] << endl;
    // set up the betas
    frequency = 0;
    betaparent = 0;
    betaparent0 = 0;
    alpha = 0;
    alpha0 = 0;
    beta = 0;
    for (int i = 0; i < betaCovariate.size(); i++)
        for (int j = 0; j < betaCovariate[i].size(); j++) {
            betaCovariate[i][j] = 0;
            betaCovariate0[i][j] = 0;
            betaparentCovariate[i][j] = 0;
            betaparentCovariate0[i][j] = 0;
        }

    for (int i = 0; i < group.size(); i++) {
        for (int j = 0; j < nhap; j++) if (!zero[j] && !(i && rare[j])) {
                if (options.model != "commonmain" ||
                        i == 0 && options.condition.size() > 0) {
                    if (group[i][j] < group[i][refIndex]) {
                        beta[j] += y[group[i][j] + ndim[i]];
                    }
                    if (group[i][j] > group[i][refIndex]) {
                        beta[j] += y[group[i][j] + ndim[i] - 1];
                    }
                } else {
                    if (group[i][j] != group[i][refIndex]) {
                        beta[j] += y[ndim[options.condition.size()>0]];
                    }
                }
            }
    }

    // set up the nuisance parameters
    int ix = y.size() - 1;
    // freq
    for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
            frequency[i] = y[ix--];
        }
    // betaparent
    if (haveBetaParent(options)) {
        for (int i = 0; i < nhap; i++)
            if (!zero[i] && i != refIndex) {
                betaparent[i] = y[ix--];
            }
        if (haveBetaParent0(options)) {
            betaparent0 = y[ix--];
        }
    }
    // alpha
    if (haveAlpha(options)) {
        for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                alpha[i] = y[ix--];
            }
        if (haveAlpha0(options)) {
            alpha0 = y[ix--];
        } else {
            alpha0 = 1;
        }
    }
    // covariates
    for (int j = 0; j < betaCovariate.size(); j++)
        for (int k = 0; k < betaCovariate[j].size(); k++) {
            for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                    betaCovariate[j][k][i] = y[ix--];
                }
            if (haveFamilies && !confounder[j] && !options.hhrr)
                for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                        betaparentCovariate[j][k][i] = y[ix--];
                    }
            if (typeOfPhenotype == "quant") {
                betaCovariate0[j][k] = y[ix--];
                if (haveFamilies && !options.hhrr) {
                    betaparentCovariate0[j][k] = y[ix--];
                }
            }
        }

    //  haplotype-based haplotype relative risk (unmatched case/control)
    if (options.hhrr) {
        betaparent = beta;
        betaparentCovariate = betaCovariate;
        betaparent0 = alpha0;
        betaparentCovariate0 = betaCovariate0;
    }

    // working score vector
    vector<double> gradient(sizeOfGradient(options), 0);

    // calculate log-likelihood and score
    double llhd = 0;
    score(options, llhd, gradient);
    if (options.follow) {
        cout << llhd << endl;
    }
    if (options.llhd) {
        for (int i = 0; i < gradient.size(); i++) {
            cout << gradient[i] << " ";
        }
        cout << endl;
    }
    neval++;
    gradientFreeParameters(gradient, g, options);

    // compute numerical gradients for debugging purposes
    if (options.checkgradient) {
        cout << endl;
        double ss = llhd;
        cout << "likelihood " << ss << endl;
        ix = 0;
        vector<double> temp = gradient;
        cout << "beta" << endl;
        for (int i = 0; i < nhap; i++) {
            cout << gradient[ix++] << " " << flush;
        }
        cout << endl;
        for (int i = 0; i < beta.size(); i++) {
            beta[i] += options.epsilon;
            if (options.hhrr) {
                betaparent = beta;
            }
            score(options, llhd, temp);
            cout << (ss - llhd) / options.epsilon << " " << flush;
            beta[i] -= options.epsilon;
            if (options.hhrr) {
                betaparent = beta;
            }
        }
        cout << endl;
        cout << "frequency" << endl;
        for (int i = 0; i < nhap; i++) {
            cout << gradient[ix++] << " " << flush;
        }
        cout << endl;
        for (int i = 0; i < frequency.size(); i++) {
            frequency[i] += options.epsilon;
            score(options, llhd, temp);
            cout << (ss - llhd) / options.epsilon << " " << flush;
            frequency[i] -= options.epsilon;
        }
        cout << endl;
        if (haveBetaParent(options)) {
            cout << "betaparent" << endl;
            for (int i = 0; i < nhap; i++) {
                cout << gradient[ix++] << " " << flush;
            }
            cout << endl;
            for (int i = 0; i < betaparent.size(); i++) {
                betaparent[i] += options.epsilon;
                score(options, llhd, temp);
                cout << (ss - llhd) / options.epsilon << " " << flush;
                betaparent[i] -= options.epsilon;
            }
            cout << endl;
            if (haveBetaParent0(options)) {
                cout << "betaparent0" << endl;
                cout << gradient[ix++] << " " << flush;
                cout << endl;
                betaparent0 += options.epsilon;
                score(options, llhd, temp);
                cout << (ss - llhd) / options.epsilon << " " << flush;
                betaparent0 -= options.epsilon;
                cout << endl;
            }
        }
        if (haveAlpha(options)) {
            cout << "alpha" << endl;
            for (int i = 0; i < nhap; i++) {
                cout << gradient[ix++] << " " << flush;
            }
            cout << endl;
            for (int i = 0; i < alpha.size(); i++) {
                alpha[i] += options.epsilon;
                score(options, llhd, temp);
                cout << (ss - llhd) / options.epsilon << " " << flush;
                alpha[i] -= options.epsilon;
            }
            cout << endl;
            if (haveAlpha0(options)) {
                cout << "alpha0" << endl;
                cout << gradient[ix++] << " " << flush;
                cout << endl;
                alpha0 += options.epsilon;
                if (options.hhrr) {
                    betaparent0 = alpha0;
                }
                score(options, llhd, temp);
                cout << (ss - llhd) / options.epsilon << " " << flush;
                alpha0 -= options.epsilon;
                if (options.hhrr) {
                    betaparent0 = alpha0;
                }
                cout << endl;
            }
        }
        for (int cov = 0; cov < betaCovariate.size(); cov++)
            for (int level = 0; level < betaCovariate[cov].size(); level++) {
                cout << "betaCovariate " << cov << " level " << level << endl;
                for (int i = 0; i < nhap; i++) {
                    cout << gradient[ix++] << " " << flush;
                }
                cout << endl;
                for (int j = 0; j < betaCovariate[cov][level].size(); j++) {
                    betaCovariate[cov][level][j] += options.epsilon;
                    if (!confounder[cov] && options.hhrr) {
                        betaparentCovariate = betaCovariate;
                    }
                    score(options, llhd, temp);
                    cout << (ss - llhd) / options.epsilon << " " << flush;
                    betaCovariate[cov][level][j] -= options.epsilon;
                    if (!confounder[cov] && options.hhrr) {
                        betaparentCovariate = betaCovariate;
                    }
                }
                cout << endl;
                if (haveFamilies && !confounder[cov] && !options.hhrr) {
                    cout << "betaparentCovariate " << cov << " level " << level << endl;
                    for (int i = 0; i < nhap; i++) {
                        cout << gradient[ix++] << " " << flush;
                    }
                    cout << endl;
                    for (int j = 0; j < betaparentCovariate[cov][level].size(); j++) {
                        betaparentCovariate[cov][level][j] += options.epsilon;
                        score(options, llhd, temp);
                        cout << (ss - llhd) / options.epsilon << " " << flush;
                        betaparentCovariate[cov][level][j] -= options.epsilon;
                    }
                    cout << endl;
                }
                if (typeOfPhenotype == "quant") {
                    cout << "betaCovariate0 " << cov << " level " << level << endl;
                    cout << gradient[ix++] << " " << flush;
                    cout << endl;
                    betaCovariate0[cov][level] += options.epsilon;
                    if (options.hhrr) {
                        betaparentCovariate0 = betaCovariate0;
                    }
                    score(options, llhd, temp);
                    cout << (ss - llhd) / options.epsilon << " " << flush;
                    cout << endl;
                    betaCovariate0[cov][level] -= options.epsilon;
                    if (options.hhrr) {
                        betaparentCovariate0 = betaCovariate0;
                    }
                    if (haveFamilies && !options.hhrr) {
                        cout << "betaparentCovariate0 " << cov << " level " << level << endl;
                        cout << gradient[ix++] << " " << flush;
                        cout << endl;
                        betaparentCovariate0[cov][level] += options.epsilon;
                        score(options, llhd, temp);
                        cout << (ss - llhd) / options.epsilon << " " << flush;
                        betaparentCovariate0[cov][level] -= options.epsilon;
                        cout << endl;
                    }
                }
            }
        cout << endl;
        exit(0);
    }

    return(-llhd);
}



//  getloglikelihood

double UnphasedAnalysis::getloglikelihood(const string &title,
        UnphasedOptions &options,
        const string &which,
        bool null) {

    // scan the data and build log-likelihood
    double loglikelihood;
    double maxllhd;
    int size = genoCode.size();
    valarray<double> bestfreq(size), bestbetaparent(size),
             bestalpha(size), bestbeta(size);
    valarray<double> bestFamilyCount[2], bestUnrelatedCount[2];
    bestFamilyCount[0].resize(size);
    bestFamilyCount[1].resize(size);
    bestUnrelatedCount[0].resize(size);
    bestUnrelatedCount[1].resize(size);
    vector<vector<valarray<double> > >
    bestbetaCovariate, bestbetaparentCovariate;
    bestbetaCovariate.resize(betaCovariate.size());
    bestbetaparentCovariate.resize(betaCovariate.size());
    for (int i = 0; i < betaCovariate.size(); i++) {
        bestbetaCovariate[i].resize(betaCovariate[i].size());
        bestbetaparentCovariate[i].resize(betaCovariate[i].size());
        for (int j = 0; j < betaCovariate[i].size(); j++) {
            bestbetaCovariate[i][j].resize(size);
            bestbetaparentCovariate[i][j].resize(size);
        }
    }

    if (which == "asymptotic") {
        *outStream << title << " likelihood..." << flush;
    }

    // allow multiple restarts
    int neval = 0; // number of likelihood evaluations
    for (int restart = 0; restart <= options.restarts; restart++) {

        loglikelihood = 0;
        bool randomise = restart || options.checkgradient;

        // solve estimating equations
        if (options.neldermead) {
            amoeba(options, randomise, neval);
        } else {
            mydfpmin(options, randomise, neval);
            if (neval < 0) {
                if (which == "asymptotic") {
                    *outStream << "(DFP algorithm failed, using Nelder-Mead, try -epsilon 1e-6)..." << flush;
                }
                neval = -neval;
                amoeba(options, restart, neval);
            }
        }

        // calculate likelihood for the converged solution
        score(options, loglikelihood, null);

        // save MLEs
        if (restart == 0 || loglikelihood > maxllhd) {
            bestfreq = frequency;
            bestbetaparent = betaparent;
            bestFamilyCount[0] = familyCount[0];
            bestFamilyCount[1] = familyCount[1];
            bestUnrelatedCount[0] = unrelatedCount[0];
            bestUnrelatedCount[1] = unrelatedCount[1];
            bestalpha = alpha;
            bestbeta = beta;
            for (int i = 0; i < betaCovariate.size(); i++)
                for (int j = 0; j < betaCovariate[i].size(); j++) {
                    bestbetaCovariate[i][j] = betaCovariate[i][j];
                    bestbetaparentCovariate[i][j] = betaparentCovariate[i][j];
                }
            maxllhd = loglikelihood;
        }
    }

    // restore MLEs
    loglikelihood = maxllhd;
    frequency = bestfreq;
    betaparent = bestbetaparent;
    familyCount[0] = bestFamilyCount[0];
    familyCount[1] = bestFamilyCount[1];
    unrelatedCount[0] = bestUnrelatedCount[0];
    unrelatedCount[1] = bestUnrelatedCount[1];
    alpha = bestalpha;
    beta = bestbeta;
    for (int i = 0; i < betaCovariate.size(); i++)
        for (int j = 0; j < betaCovariate[i].size(); j++) {
            betaCovariate[i][j] = bestbetaCovariate[i][j];
            betaparentCovariate[i][j] = bestbetaparentCovariate[i][j];
        }
    if (which == "asymptotic" || options.outputPermutation) {
        *outStream << "done (" << neval << " evaluations)" << endl;
    }

    return(loglikelihood);

}



