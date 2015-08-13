/* UnphasedScoreUnrelated.cpp - Association analysis in unrelateds

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

#include <iostream>
#include <cmath>
#include <string>
#include "stl.h"
#include "UnphasedAnalysis.h"
#include "UnphasedOptions.h"
#include "pvalues.h"
#include "asran.h"
#include "matinv.h"

using namespace std;

inline bool seq(int x, int y) {
    return(x != 0 && y != 0 && x == y);
}
inline bool weq(int x, int y) {
    return(x == 0 || y == 0 || x == y);
}

template<class Type>
void debug(Type t) {
    for (int i = 0; i < t.size(); i++) {
        cout << t[i] << " ";
    }
    cout << endl;
}

//  scoreUnrelated

void UnphasedAnalysis::scoreUnrelated(Subject &subject, int nsubject,
                                      double &subjectTrait,
                                      vector<double> &subjectCovariate,
                                      UnphasedOptions &options, double &prob,
                                      vector<double> &gradient,
                                      bool *haveGlobal, double *globalDenom,
                                      vector<double> *globalGrad) {

    int nhap = genoCode.size(); // number of haplotypes
    // return if no haplotypes have positive frequency
    bool somehap = false;
    for (int i = 0; i < nhap; i++) if (!zero[i]) {
            somehap = true;
        }
    if (!somehap) {
        return;
    }

    int nconfounder = options.confounder.size(); // number of confounders
    int nmodifier = options.modifier.size(); // number of modifiers

    realisationProb[nsubject].clear(); // probs of different realisations
    bool haveFamilies = familylist.size() > 0;

    bool normal = options.normal && typeOfPhenotype == "quant";

    // situations where we can use the fast algorithm
    bool fast = !normal && !options.slow &&
                !options.condgenotype && !options.genotype;

    // situations where only homozygotes will be counted
    bool onlyHomoz = options.certain && options.condition.size() + options.window + options.tag.size() > 1 ||
                     options.chrX && subject.sex == MALE ||
                     options.chrY;

    // situations where we would only count one chromosome
    bool oneChr = onlyHomoz || options.genotype;

    // denominator first
    vector<double> denom;
    double denomSum;
    bool useSavedDenom[2];

    // if a control, use saved denominator
    useSavedDenom[0] = (!options.chrX && subjectTrait == 0 && !subjectCovariate.size());
    // if a case, use saved denominator
    useSavedDenom[1] = (!options.chrX && subjectTrait == 1 && !subjectCovariate.size());
    for (int i = 0; i < 2; i++)
        if (haveGlobal[i] && useSavedDenom[i]) {
            denomSum = globalDenom[i];
        }

    // otherwise calculate it from scratch
    if (!(haveGlobal[0] && useSavedDenom[0]) &&
            !(haveGlobal[1] && useSavedDenom[1])) {


        if (fast) {
            for (int i = 0; i < nhap; i++) if (!zero[i]) {
                    double freq = 0, linear = 0, alphaterm = 0;
                    freq = frequency[i];
                    linear = beta[i];
                    alphaterm = alpha0;
                    for (int cov = 0; cov < subjectCovariate.size(); cov++) {
                        int level = 0;
                        double weight = 0;
                        if (factor[cov] && subjectCovariate[cov] != 0) {
                            level = (int)subjectCovariate[cov] - 1;
                            weight = 1.0;
                        } else {
                            weight = subjectCovariate[cov];
                        }
                        double x = weight * betaCovariate[cov][level][i];
                        if (confounder[cov]) {
                            freq += x;
                        } else {
                            linear += x;
                        }
                        alphaterm += betaCovariate0[cov][level];
                    }
                    denom.push_back(freq + (subjectTrait * linear - alphaterm*(linear + alpha[i])) /
                                    options.variance);
                }
        } else {
            for (int i = 0; i < haploCode.size(); i++) if (!zero[i])
                    for (int j = 0; j < haploCode.size(); j++) if (!zero[j] && (!onlyHomoz || j == i)) {
                            int ncondition = options.condgenotype * options.condition.size();
                            Haplotype hap1 = haploCode[i].formgeno(haploCode[j], ncondition, options.genotype);
                            Haplotype hap2 = haploCode[j].formgeno(haploCode[i], ncondition, options.genotype);
                            int which1 = hap1.index(genoCode);
                            int which2 = hap2.index(genoCode);
                            if (which1 < nhap && !zero[which1] &&
                                    which2 < nhap && !zero[which2]) {
                                double freq = 0, linear = 0, alphaterm = 0;
                                freq += frequency[which1];
                                linear += beta[which1];
                                alphaterm += alpha0;
                                if (!oneChr) {
                                    freq += frequency[which2];
                                    linear += beta[which2];
                                }

                                for (int cov = 0; cov < subjectCovariate.size(); cov++) {
                                    int level = 0;
                                    double weight = 0;
                                    if (factor[cov] && subjectCovariate[cov] != 0) {
                                        level = (int)subjectCovariate[cov] - 1;
                                        weight = 1.0;
                                    } else {
                                        weight = subjectCovariate[cov];
                                    }
                                    double x1 = weight * betaCovariate[cov][level][which1];
                                    double x2 = weight * betaCovariate[cov][level][which2];
                                    if (confounder[cov]) {
                                        freq += x1;
                                    } else {
                                        linear += x1;
                                    }
                                    alphaterm += betaCovariate0[cov][level];
                                    if (!oneChr) {
                                        if (confounder[cov]) {
                                            freq += x2;
                                        } else {
                                            linear += x2;
                                        }
                                    }
                                }

                                double beta1 = subjectTrait * linear - alphaterm * (linear + alpha[which1] + (!oneChr) * alpha[which2]);
                                if (normal) {
                                    beta1 -= linear * linear / 2.0;
                                }

                                denom.push_back(freq + beta1 / options.variance);
                            }
                        }
        }

        if (denom.size() > 0) {
            denomSum = 0;
            double maxbeta = *(max_element(denom.begin(), denom.end()));
            for (int i = 0; i < denom.size(); i++) {
                denomSum += exp(denom[i] - maxbeta);
            }
            // denomSum must be at least 1
            denomSum = log(denomSum) + maxbeta;
            if (!oneChr && fast) {
                denomSum *= 2;
            }
        }
    }

    // now the numerator
    vector<double> numer;
    for (int i = 0; i < consistentHaps[nsubject].size(); i += 2) {
        int h1 = consistentHaps[nsubject][i];
        int h2 = consistentHaps[nsubject][i+1];
        if (!zero[h1] && (oneChr || !zero[h2])) {
            double freq = 0, linear = 0, alphaterm = 0;
            freq += frequency[h1];
            linear += beta[h1];
            alphaterm += alpha0;
            if (!oneChr) {
                freq += frequency[h2];
                linear += beta[h2];
            }
            for (int cov = 0; cov < subjectCovariate.size(); cov++) {
                int level = 0;
                double weight = 0;
                if (factor[cov] && subjectCovariate[cov] != 0) {
                    level = (int)subjectCovariate[cov] - 1;
                    weight = 1.0;
                } else {
                    weight = subjectCovariate[cov];
                }
                double x = weight * betaCovariate[cov][level][h1];
                if (confounder[cov]) {
                    freq += x;
                } else {
                    linear += x;
                }
                alphaterm += betaCovariate0[cov][level];
                if (!oneChr) {
                    x = weight * betaCovariate[cov][level][h2];
                    if (confounder[cov]) {
                        freq += x;
                    } else {
                        linear += x;
                    }
                }
            }

            double beta1 = subjectTrait * linear - alphaterm * (linear + alpha[h1] + (!oneChr) * alpha[h2]);
            if (normal) {
                beta1 -= linear * linear / 2.0;
            }

            numer.push_back(freq + beta1 / options.variance - denomSum);
        }
    }

    // form the numerator (which is already divided by denominator :)
    double numerSum = 0;
    if (numer.size()) {
        double maxnumer = *(max_element(numer.begin(), numer.end()));
        for (int i = 0; i < numer.size(); i++) {
            if (options.llhd) {
                cout << "numer " << numer[i] << " denom " << denomSum << endl;
            }
            numerSum += exp(numer[i] - maxnumer);
        }

        // the log-likelihood contribution
        numerSum = log(numerSum) + maxnumer;
        prob = numerSum;
    }
    if (!numer.size()) {
        return;
    }

    // gradients
    vector<double> betaGradient(nhap, 0);
    vector<double> freqGradient(nhap, 0);
    vector<double> alphaGradient(nhap, 0);
    double alpha0Gradient = 0;
    vector<vector<vector<double> > > betaCovariateGradient(betaCovariate.size());
    vector<vector<double> > betaCovariate0Gradient(betaCovariate.size());
    for (int i = 0; i < betaCovariate.size(); i++) {
        betaCovariateGradient[i].resize(betaCovariate[i].size());
        betaCovariate0Gradient[i].resize(betaCovariate0[i].size(), 0);
        for (int j = 0; j < betaCovariate[i].size(); j++) {
            betaCovariateGradient[i][j].resize(nhap, 0);
        }
    }

    // contributions to gradient from denominator
    for (int i = 0; i < 2; i++)
        if (haveGlobal[i] && useSavedDenom[i]) {
            gradient = globalGrad[i];
        }

    if (!(haveGlobal[0] && useSavedDenom[0]) &&
            !(haveGlobal[1] && useSavedDenom[1])) {

        int ix = 0;
        if (fast) {
            for (int i = 0; i < nhap; i++) if (!zero[i]) {
                    double x = exp(denom[ix++] - denomSum / (2 - oneChr)) * (2 - oneChr);
                    double linear = 0, alphaterm = 0;
                    freqGradient[i] += x;
                    betaGradient[i] += subjectTrait * x / options.variance;
                    linear += beta[i];
                    alphaterm = alpha0;
                    for (int cov = 0; cov < subjectCovariate.size(); cov++) {
                        int level = 0;
                        double weight = 0;
                        if (factor[cov] && subjectCovariate[cov] != 0) {
                            level = (int)subjectCovariate[cov] - 1;
                            weight = 1.0;
                        } else {
                            weight = subjectCovariate[cov];
                        }
                        if (confounder[cov]) {
                            betaCovariateGradient[cov][level][i] += weight * x;
                        } else {
                            betaCovariateGradient[cov][level][i] += subjectTrait * weight * x / options.variance;
                            linear += weight * betaCovariate[cov][level][i];
                        }
                        alphaterm += betaCovariate0[cov][level];
                    }
                    betaGradient[i] -= alphaterm * x / options.variance;
                    alphaGradient[i] -= alphaterm * x / options.variance;
                    alpha0Gradient -= (linear + alpha[i]) * x / options.variance;
                    for (int cov = 0; cov < subjectCovariate.size(); cov++) {
                        int level = 0;
                        double weight = 0;
                        if (factor[cov] && subjectCovariate[cov] != 0) {
                            level = (int)subjectCovariate[cov] - 1;
                            weight = 1.0;
                        } else {
                            weight = subjectCovariate[cov];
                        }
                        if (!confounder[cov]) {
                            betaCovariateGradient[cov][level][i] -= alphaterm * weight * x / options.variance;
                        }
                        betaCovariate0Gradient[cov][level] -= (linear + alpha[i]) * x / options.variance;
                    }
                }
        } else {
            for (int i = 0; i < haploCode.size(); i++) if (!zero[i])
                    for (int j = 0; j < haploCode.size(); j++) if (!zero[j] && (!onlyHomoz || j == i)) {
                            int ncondition = options.condgenotype * options.condition.size();
                            Haplotype hap1 = haploCode[i].formgeno(haploCode[j], ncondition, options.genotype);
                            Haplotype hap2 = haploCode[j].formgeno(haploCode[i], ncondition, options.genotype);
                            int which1 = hap1.index(genoCode);
                            int which2 = hap2.index(genoCode);
                            if (which1 < nhap && !zero[which1] &&
                                    which2 < nhap && !zero[which2]) {
                                double x = exp(denom[ix++] - denomSum);
                                double linear = 0, alphaterm = 0;
                                freqGradient[which1] += x;
                                betaGradient[which1] += subjectTrait * x / options.variance;
                                linear += beta[which1];
                                alphaterm += alpha0;
                                if (!oneChr) {
                                    freqGradient[which2] += x;
                                    betaGradient[which2] += subjectTrait * x / options.variance;
                                    linear += beta[which2];
                                }
                                for (int cov = 0; cov < subjectCovariate.size(); cov++) {
                                    int level = 0;
                                    double weight = 0;
                                    if (factor[cov] && subjectCovariate[cov] != 0) {
                                        level = (int)subjectCovariate[cov] - 1;
                                        weight = 1.0;
                                    } else {
                                        weight = subjectCovariate[cov];
                                    }
                                    if (confounder[cov]) {
                                        betaCovariateGradient[cov][level][which1] += weight * x;
                                    } else {
                                        betaCovariateGradient[cov][level][which1] += subjectTrait * weight * x / options.variance;
                                        linear += weight * betaCovariate[cov][level][which1];
                                    }
                                    alphaterm += betaCovariate0[cov][level];
                                    if (!oneChr) {
                                        if (confounder[cov]) {
                                            betaCovariateGradient[cov][level][which2] += weight * x;
                                        } else {
                                            betaCovariateGradient[cov][level][which2] += subjectTrait * weight * x / options.variance;
                                            linear += weight * betaCovariate[cov][level][which2];
                                        }
                                    }
                                }
                                betaGradient[which1] -= alphaterm * x / options.variance;
                                alphaGradient[which1] -= alphaterm * x / options.variance;
                                alpha0Gradient -= (linear + alpha[which1]) * x / options.variance;
                                if (!oneChr) {
                                    betaGradient[which2] -= alphaterm * x / options.variance;
                                    alphaGradient[which2] -= alphaterm * x / options.variance;
                                    alpha0Gradient -= alpha[which2] * x / options.variance;
                                }
                                for (int cov = 0; cov < subjectCovariate.size(); cov++) {
                                    int level = 0;
                                    double weight = 0;
                                    if (factor[cov] && subjectCovariate[cov] != 0) {
                                        level = (int)subjectCovariate[cov] - 1;
                                        weight = 1.0;
                                    } else {
                                        weight = subjectCovariate[cov];
                                    }
                                    if (!confounder[cov]) {
                                        betaCovariateGradient[cov][level][which1] -= alphaterm * weight * x / options.variance;
                                        if (!oneChr) {
                                            betaCovariateGradient[cov][level][which2] -= alphaterm * weight * x / options.variance;
                                        }
                                    }
                                    betaCovariate0Gradient[cov][level] -= (linear + alpha[which1]) * x / options.variance;
                                    if (!oneChr) {
                                        betaCovariate0Gradient[cov][level] -= alpha[which2] * x / options.variance;
                                    }
                                }
                                if (normal) {
                                    betaGradient[which1] -= linear * x / options.variance;
                                    if (!oneChr) {
                                        betaGradient[which2] -= linear * x / options.variance;
                                    }
                                    for (int cov = 0; cov < subjectCovariate.size(); cov++) {
                                        int level = 0;
                                        double weight = 0;
                                        if (factor[cov] && subjectCovariate[cov] != 0) {
                                            level = (int)subjectCovariate[cov] - 1;
                                            weight = 1.0;
                                        } else {
                                            weight = subjectCovariate[cov];
                                        }
                                        betaCovariateGradient[cov][level][which1] -= linear * weight * x / options.variance;
                                        if (!oneChr) {
                                            betaCovariateGradient[cov][level][which2] -= linear * weight * x / options.variance;
                                        }
                                    }
                                }
                            }
                        }
        }
    }

    // save global denominator
    for (int i = 0; i < 2; i++)
        if (!haveGlobal[i] && useSavedDenom[i]) {
            haveGlobal[i] = true;
            globalDenom[i] = denomSum;
            globalGrad[i].resize(gradient.size(), 0);
            int ix = 0;
            for (int j = 0; j < nhap; j++) {
                globalGrad[i][ix++] = betaGradient[j];
            }
            for (int j = 0; j < nhap; j++) {
                globalGrad[i][ix++] = freqGradient[j];
            }
            if (haveBetaParent(options)) {
                ix += nhap + haveBetaParent0(options);
            }
            if (haveAlpha(options)) {
                for (int j = 0; j < nhap; j++) {
                    globalGrad[i][ix++] = alphaGradient[j];
                }
                if (haveAlpha0(options)) {
                    globalGrad[i][ix++] = alpha0Gradient;
                }
            }
            for (int j = 0; j < betaCovariate.size(); j++)
                for (int k = 0; k < betaCovariate[j].size(); k++) {
                    for (int l = 0; l < nhap; l++) {
                        globalGrad[i][ix++] = betaCovariateGradient[j][k][l];
                    }
                    if (!confounder[i] && haveFamilies && !options.hhrr) {
                        ix += nhap;
                    }
                    if (typeOfPhenotype == "quant") {
                        globalGrad[i][ix++] = betaCovariate0Gradient[j][k];
                        if (haveFamilies && !options.hhrr) {
                            ix++;
                        }
                    }
                }
        }
    // contributions to gradient from numerator
    int ix = 0;
    for (int i = 0; i < consistentHaps[nsubject].size(); i += 2) {
        int h1 = consistentHaps[nsubject][i];
        int h2 = consistentHaps[nsubject][i+1];
        if (h1 >= 0 && !zero[h1] && (oneChr || h2 >= 0 && !zero[h2])) {
            double x = 0;
            double linear = 0, alphaterm = 0;
            x = exp(numer[ix++] - numerSum);
            realisationProb[nsubject].push_back(h1);
            realisationProb[nsubject].push_back(h2);
            realisationProb[nsubject].push_back(x);
            betaGradient[h1] -= subjectTrait * x / options.variance;
            freqGradient[h1] -= x;
            linear += beta[h1];
            alphaterm += alpha0;
            for (int cov = 0; cov < subjectCovariate.size(); cov++) {
                int level = 0;
                double weight = 0;
                if (factor[cov] && subjectCovariate[cov] != 0) {
                    level = (int)subjectCovariate[cov] - 1;
                    weight = 1.0;
                } else {
                    weight = subjectCovariate[cov];
                }
                if (confounder[cov]) {
                    betaCovariateGradient[cov][level][h1] -= weight * x;
                } else {
                    betaCovariateGradient[cov][level][h1] -= subjectTrait * weight * x / options.variance;
                    linear += weight * betaCovariate[cov][level][h1];
                }
                alphaterm += betaCovariate0[cov][level];
            }
            unrelatedCount[typeOfPhenotype!="quant" && subjectTrait][h1] += x;
            if (!oneChr) {
                betaGradient[h2] -= subjectTrait * x / options.variance;
                freqGradient[h2] -= x;
                linear += beta[h2];
                for (int cov = 0; cov < subjectCovariate.size(); cov++) {
                    int level = 0;
                    double weight = 0;
                    if (factor[cov] && subjectCovariate[cov] != 0) {
                        level = (int)subjectCovariate[cov] - 1;
                        weight = 1.0;
                    } else {
                        weight = subjectCovariate[cov];
                    }
                    if (confounder[cov]) {
                        betaCovariateGradient[cov][level][h2] -= weight * x;
                    } else {
                        betaCovariateGradient[cov][level][h2] -= subjectTrait * weight * x / options.variance;
                        linear += weight * betaCovariate[cov][level][h2];
                    }
                }
                unrelatedCount[typeOfPhenotype!="quant" && subjectTrait][h2] += x;
            }
            betaGradient[h1] += alphaterm * x / options.variance;
            alphaGradient[h1] += alphaterm * x / options.variance;
            alpha0Gradient += (linear + alpha[h1]) * x / options.variance;
            if (!oneChr) {
                betaGradient[h2] += alphaterm * x / options.variance;
                alphaGradient[h2] += alphaterm * x / options.variance;
                alpha0Gradient += alpha[h2] * x / options.variance;
            }
            for (int cov = 0; cov < subjectCovariate.size(); cov++) {
                int level = 0;
                double weight = 0;
                if (factor[cov] && subjectCovariate[cov] != 0) {
                    level = (int)subjectCovariate[cov] - 1;
                    weight = 1.0;
                } else {
                    weight = subjectCovariate[cov];
                }
                if (!confounder[cov]) {
                    betaCovariateGradient[cov][level][h1] += alphaterm * weight * x / options.variance;
                    if (!oneChr) {
                        betaCovariateGradient[cov][level][h2] += alphaterm * weight * x / options.variance;
                    }
                }
                betaCovariate0Gradient[cov][level] += (linear + alpha[h1]) * x / options.variance;
                if (!oneChr) {
                    betaCovariate0Gradient[cov][level] += alpha[h2] * x / options.variance;
                }
            }
            if (normal) {
                betaGradient[h1] += linear * x / options.variance;
                if (!oneChr) {
                    betaGradient[h2] += linear * x / options.variance;
                }
                for (int cov = 0; cov < subjectCovariate.size(); cov++) {
                    int level = 0;
                    double weight = 0;
                    if (factor[cov] && subjectCovariate[cov] != 0) {
                        level = (int)subjectCovariate[cov] - 1;
                        weight = 1.0;
                    } else {
                        weight = subjectCovariate[cov];
                    }
                    betaCovariateGradient[cov][level][h1] += linear * weight * x / options.variance;
                    if (!oneChr) {
                        betaCovariateGradient[cov][level][h2] += linear * weight * x / options.variance;
                    }
                }
            }

        } // if ((h1<0 || !zero[h1]) && (h2<0 || !zero[h2]))
    } // for(int i=0;i<consistentHaps[nsubject].size();i+=2)

    ix = 0;
    for (int i = 0; i < nhap; i++) {
        gradient[ix++] += betaGradient[i];
    }
    for (int i = 0; i < nhap; i++) {
        gradient[ix++] += freqGradient[i];
    }
    if (haveBetaParent(options)) {
        ix += nhap + haveBetaParent0(options);
    }
    if (haveAlpha(options)) {
        for (int i = 0; i < nhap; i++) {
            gradient[ix++] += alphaGradient[i];
        }
        if (haveAlpha0(options)) {
            gradient[ix++] += alpha0Gradient;
        }
    }
    for (int i = 0; i < betaCovariate.size(); i++)
        for (int j = 0; j < betaCovariate[i].size(); j++) {
            for (int k = 0; k < nhap; k++) {
                gradient[ix++] += betaCovariateGradient[i][j][k];
            }
            if (!confounder[i] && haveFamilies && !options.hhrr) {
                ix += nhap;
            }
            if (typeOfPhenotype == "quant") {
                gradient[ix++] += betaCovariate0Gradient[i][j];
                if (haveFamilies && !options.hhrr) {
                    ix++;
                }
            }
        }
}


