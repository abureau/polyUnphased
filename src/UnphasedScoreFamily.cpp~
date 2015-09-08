/* UnphasedScoreLinkage.cpp - Association analysis in nuclear families

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

//  scoreFamily

void UnphasedAnalysis::scoreFamily(NuclearFamily &family, int nfamily,
                                   vector<double> &sibTrait,
                                   vector<int> &sibSex,
                                   vector<vector<double> > &sibCovariate,
                                   UnphasedOptions &options, double &prob,
                                   vector<double> &gradient,
                                   bool *haveGlobal, double *globalDenom,
                                   vector<double> *globalGrad) {

    // count the number of siblings used in this family
    int nsib = sibTrait.size();

    // no sibs to score in this family
    if (nsib == 0) {
        return;
    }

    // no consistent solutions for all sibs in this family
    // could happen when there are recombinant haplotypes
    if (consistentHaps[nfamily][consistentHaps[nfamily].size() - 1] == -1) {
        return;
    }

    // using normal likelihood model
    bool normal = options.normal && typeOfPhenotype == "quant";

    // use two virtual controls instead of four
    bool twoVirtualControls = homix[nfamily] || options.chrX || options.onefbc;

    // covariance and variance terms from inverted var/cov matrix
    double invCov = options.covariance / ((nsib - 1) * options.covariance * options.covariance - options.variance * (options.variance + (nsib - 2) * options.covariance));
    double invVar = (1 - (nsib - 1) * options.covariance * invCov) / options.variance;
    // no covariance if no linkage
    if (options.nolinkage) {
        invCov = 0;
        invVar = 1;
    }

    // sum of intercepts
    double sumIntercept = 0;
    sumIntercept = alpha0;
    for (int cov = 0; cov < betaCovariate.size(); cov++)
        for (int level = 0; level < betaCovariate[cov].size(); level++) {
            sumIntercept += betaCovariate0[cov][level];
        }
    double parentSumIntercept = 0;
    parentSumIntercept = betaparent0;
    for (int cov = 0; cov < betaCovariate.size(); cov++)
        for (int level = 0; level < betaCovariate[cov].size(); level++) {
            parentSumIntercept += betaparentCovariate0[cov][level];
        }

    int nphase = 1 << (nsib - 1); // 2^(nsib-1) = number of phases in sibship
    int nhap = genoCode.size(); // number of haplotypes
    // return if no haplotypes have positive frequency
    int nhapz = 0;
    for (int i = 0; i < nhap; i++) if (!zero[i]) {
            nhapz++;
        }
    if (!nhapz) {
        return;
    }

    int betasize;
    int nconfounder = options.confounder.size(); // number of confounders
    int nmodifier = options.modifier.size(); // number of modifiers

    realisationProb[nfamily].clear(); // probs of different realisations

    // list of whether configurations have a positive probability
    deque<bool> positiveProbQueue;

    // these arrays hold linear models for single chromosomes
    // to be used in rapid calculation of the likelihood denominator
    // store one for each phase
    vector<vector<double> > dft, dfnt, dmt, dmnt;
    dft.resize(nphase);
    dfnt.resize(nphase);
    dmt.resize(nphase);
    dmnt.resize(nphase);
    for (int i = 0; i < nphase; i++) {
        dft[i].resize(nhap, 0);
        dfnt[i].resize(nhap, 0);
        dmt[i].resize(nhap, 0);
        dmnt[i].resize(nhap, 0);
    }

    // sums of linear models for single chromosomes
    vector<double> sumft(nphase, 0);
    vector<double> sumfnt(nphase, 0);
    vector<double> summt(nphase, 0);
    vector<double> summnt(nphase, 0);

    // denominator first
    double denomSum = 0;
    vector<double> denom;

    bool useSavedDenom[2];
    // Dans ce qui suit, on semble supposer que le trait est dichotomique sans l'avoir expressément vérifié
    // Cela suppose qu'un trait quantitatif ne prend jamais exactement la valeur 1.
    // if a single affected, use saved denominator
    useSavedDenom[0] = (nsib == 1 && ((typeOfPhenotype == "binary" && sibTrait[0] == 1) || (typeOfPhenotype == "polytomous" && sibTrait[0] < K))  && !sibCovariate[0].size() && !options.chrX);
    // if an affected sib pair, use saved denominator
    useSavedDenom[1] = (nsib == 2 && ((typeOfPhenotype == "binary" && sibTrait[0] == 1 && sibTrait[1] == 1) || (typeOfPhenotype == "polytomous" && sibTrait[0] < K && sibTrait[0] == sibTrait[1])) && !sibCovariate[0].size() && !options.chrX);
    for (int i = 0; i < 2; i++)
        if (haveGlobal[i] && useSavedDenom[i]) {
            denomSum = globalDenom[i];
        }

    // otherwise calculate it from scratch
    if (!(haveGlobal[0] && useSavedDenom[0]) &&
            !(haveGlobal[1] && useSavedDenom[1])) {

        //  fast algorithm

        if (!normal && !options.slow &&
                !options.condgenotype && !options.genotype) {

            bool virgin = true; // the first time
            vector<double> maxbeta0(nphase, 0);
            vector<double> maxbeta1(nphase, 0);
            vector<double> maxbeta2(nphase, 0);
            vector<double> maxbeta3(nphase, 0);
            vector<vector<double> > maxbetaq(nphase);
            for (int i = 0; i < nphase; i++) {
                maxbetaq[i].resize(nphase, 0);
            }

            for (int i = 0; i < nhap; i++) if (!zero[i]) {
                    for (int phase = 0; phase < nphase; phase++)
                        for (int parsex = MALE; parsex <= FEMALE; parsex++) {
			    double linear = 0;
			    double sibTerm = 0;
                            double freq = 0;
                            double beta1 = 0;
                            double beta2 = 0;
                            for (int sib = 0; sib < nsib; sib++)
                                if (!(options.chrX && parsex == MALE && sibSex[sib] == MALE)) {
                                    freq = frequency[i] / nsib * (family.sibship ? 2 : 1);
                                    if (typeOfPhenotype == "polytomous") {
                                    	if (sibTrait[sib] < K) linear = betaparent[i + sibTrait[sib]*nhap];
                                    	}
                                    else linear = betaparent[i];
                                    for (int cov = 0; cov < covariateName.size(); cov++) {
                                        int level = 0;
                                        double weight = 0;
                                        if (factor[cov] && sibCovariate[sib][cov] != 0) {
                                            level = (int)sibCovariate[sib][cov] - 1;
                                            weight = 1.0;
                                        } else {
                                            weight = sibCovariate[sib][cov];
                                            if (!family.sibship && covariateName[cov] == "parsex") {
                                                weight = (baseline[cov] != parsex);
                                            }
                                        }
                                        if (confounder[cov]) {
                                            freq += weight * betaCovariate[cov][level][i] / nsib * (family.sibship ? 2 : 1);
                                        } else {
                                            linear += weight * betaparentCovariate[cov][level][i];
                                        }
                                    }
                                    if (typeOfPhenotype == "polytomous") sibTerm = linear;
                                    else sibTerm = invVar * (sibTrait[sib] - parentSumIntercept) * linear;
                                    if (options.hhrr) {
                                        sibTerm -= invVar * parentSumIntercept * alpha[i];
                                    }
                                    if (invCov != 0)
                                        for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                sibTerm += invCov * (sibTrait[sib2] - parentSumIntercept) * linear;
                                                if (options.hhrr) {
                                                    sibTerm -= invCov * parentSumIntercept * alpha[i];
                                                }
                                            }
                                    if (!options.chrX || parsex == FEMALE) {
                                        if ((phase & 1 << sib) == 0) {
                                            beta1 += freq + sibTerm;
                                            if (!family.sibship) {
                                                beta2 += freq;
                                            }
                                        } else {
                                            beta2 += freq + sibTerm;
                                            if (!family.sibship) {
                                                beta1 += freq;
                                            }
                                        }
                                    } else {
                                        beta1 += freq + sibTerm;
                                        beta2 += freq + sibTerm;
                                    }
                                } // for sib
                            if (parsex == MALE) {
                                double x = beta1;
                                if (virgin || x > maxbeta0[phase]) {
                                    maxbeta0[phase] = x;
                                }
                                dft[phase][i] = x;
                                x = beta2;
                                if (virgin || x > maxbeta1[phase]) {
                                    maxbeta1[phase] = x;
                                }
                                dfnt[phase][i] = x;
                            } else {
                                double x = beta1;
                                if (virgin || x > maxbeta2[phase]) {
                                    maxbeta2[phase] = x;
                                }
                                dmt[phase][i] = x;
                                x = beta2;
                                if (virgin || x > maxbeta3[phase]) {
                                    maxbeta3[phase] = x;
                                }
                                dmnt[phase][i] = x;
                            }
                        }
                    virgin = false;
                } // for i

            // set up sumt, sumnt
            sumft.resize(nphase, 0);
            sumfnt.resize(nphase, 0);
            summt.resize(nphase, 0);
            summnt.resize(nphase, 0);
            for (int i = 0; i < nhap; i++) if (!zero[i]) {
                    for (int phase = 0; phase < nphase; phase++) {
                        sumft[phase] += exp(dft[phase][i] - maxbeta0[phase]);
                        sumfnt[phase] += exp(dfnt[phase][i] - maxbeta1[phase]);
                        summt[phase] += exp(dmt[phase][i] - maxbeta2[phase]);
                        summnt[phase] += exp(dmnt[phase][i] - maxbeta3[phase]);
                    }
                }
            for (int i = 0; i < nphase; i++) {
                sumft[i] = log(sumft[i]) + maxbeta0[i];
                sumfnt[i] = log(sumfnt[i]) + maxbeta1[i];
                summt[i] = log(summt[i]) + maxbeta2[i];
                summnt[i] = log(summnt[i]) + maxbeta3[i];
            }

            // form the denominator sums
            // sum the denomSums over all paternal and maternal phases
            for (int i = 0; i < nphase; i++)
                for (int j = 0; j < nphase; j++) {
                    denom.push_back(sumft[i] + (!options.chrX)*sumfnt[i] + summt[j] + summnt[j]);
                }

            denomSum = 0;
            if (denom.size() > 0) {
                double maxDenom = *(max_element(denom.begin(), denom.end()));
                for (int i = 0; i < denom.size(); i++) {
                    denomSum += exp(denom[i] - maxDenom);
                }
                denomSum = log(denomSum) + maxDenom;
                if (options.llhd) {
                    cout << "denomSum " << denomSum << endl;
                }
            }

        } // fast algorithm



        //  slow algorithm, used for genotype tests

        if (normal || options.slow ||
                options.condgenotype || options.genotype) {
            denom.clear();
            for (int mt = 0; mt < haploCode.size(); mt++)
                for (int mnt = 0; mnt < haploCode.size(); mnt++)
                    for (int ft = 0; ft < haploCode.size(); ft++)
                        for (int fnt = 0; fnt < haploCode.size(); fnt++) if (!options.chrX || ft == fnt) {

                                for (int i = 0; i < nphase; i++)
                                    for (int j = 0; j < nphase; j++) {
                                        double freq = 0; // frequency
                                        bool positiveProb = true; // true if all haplotypes have +ve frequency
                                        vector<double> linear(nsib, 0); // linear models in sibs
                                        vector<double> alphaterm(nsib, 0); // alpha terms in sibs

                                        for (int sib = 0; sib < nsib; sib++) {
                                            bool fphase = (i & 1 << sib) == 0;
                                            bool mphase = (j & 1 << sib) == 0;
                                            int whichft = ft;
                                            int whichfnt = fnt;
                                            int whichmt = mt;
                                            int whichmnt = mnt;
                                            if (!fphase) {
                                                whichft = fnt;
                                                whichfnt = ft;
                                            }
                                            if (!mphase) {
                                                whichmt = mnt;
                                                whichmnt = mt;
                                            }

                                            bool MchrX = options.chrX && sibSex[sib] == MALE;

                                            // get the genotype code for this child
                                            if (options.condgenotype || options.genotype) {
                                                int ncondition = options.condgenotype * options.condition.size();
                                                Haplotype fthap = haploCode[whichft];
                                                Haplotype fnthap = haploCode[whichfnt];
                                                Haplotype mthap = haploCode[whichmt];
                                                Haplotype mnthap = haploCode[whichmnt];
                                                Haplotype ftgeno = fthap.formgeno(mthap, ncondition, options.genotype);
                                                Haplotype fntgeno = fnthap.formgeno(mnthap, ncondition, options.genotype);
                                                Haplotype mtgeno = mthap.formgeno(fthap, ncondition, options.genotype);
                                                Haplotype mntgeno = mnthap.formgeno(fnthap, ncondition, options.genotype);
                                                if (MchrX) {
                                                    ftgeno = mthap.formgeno(mthap, ncondition, options.genotype);
                                                    fntgeno = mnthap.formgeno(mnthap, ncondition, options.genotype);
                                                    mtgeno = mthap.formgeno(mthap, ncondition, options.genotype);
                                                    mntgeno = mnthap.formgeno(mnthap, ncondition, options.genotype);
                                                }
                                                whichft = ftgeno.index(genoCode);
                                                whichfnt = fntgeno.index(genoCode);
                                                whichmt = mtgeno.index(genoCode);
                                                whichmnt = mntgeno.index(genoCode);
                                            }

                                            if ((options.genotype ||
                                                    whichft < nhap && !zero[whichft] &&
                                                    whichfnt < nhap && !zero[whichfnt]) &&
                                                    whichmt < nhap && !zero[whichmt] &&
                                                    whichmnt < nhap && !zero[whichmnt]) {} // this is ok
                                            else {
                                                positiveProb = false;
                                            }

                                            if (positiveProb) {
                                                if (!options.genotype) {
                                                    if (!MchrX) {
                                                        freq += frequency[whichft] / nsib * (family.sibship ? 2 : 1);
                                                    }
                                                    if (!options.chrX && !family.sibship) {
                                                        freq += frequency[whichfnt] / nsib * (family.sibship ? 2 : 1);
                                                    }
                                                }
                                                freq += (frequency[whichmt] + (!family.sibship) * frequency[whichmnt]) / nsib * (family.sibship ? 2 : 1);
                                                if (!options.genotype && !MchrX) {
                                                    linear[sib] += betaparent[whichft];
                                                    if (options.hhrr) {
                                                        alphaterm[sib] += alpha[whichft];
                                                    }
                                                }
                                                linear[sib] += betaparent[whichmt];
                                                if (options.hhrr) {
                                                    alphaterm[sib] += alpha[whichmt];
                                                }
                                                for (int cov = 0; cov < covariateName.size(); cov++) {
                                                    double fweight = 0, mweight = 0;
                                                    int level = 0;
                                                    if (factor[cov] && sibCovariate[sib][cov] != 0) {
                                                        level = (int)sibCovariate[sib][cov] - 1;
                                                        fweight = 1.0;
                                                        mweight = 1.0;
                                                    } else {
                                                        fweight = sibCovariate[sib][cov];
                                                        mweight = sibCovariate[sib][cov];
                                                        if (!family.sibship && covariateName[cov] == "parsex") {
                                                            fweight = (baseline[cov] != MALE);
                                                            mweight = (baseline[cov] != FEMALE);
                                                        }
                                                    }
                                                    if (confounder[cov]) {
                                                        if (!options.genotype) {
                                                            if (!MchrX) {
                                                                freq += fweight * betaCovariate[cov][level][whichft] / nsib * (family.sibship ? 2 : 1);
                                                            }
                                                            if (!options.chrX && !family.sibship) {
                                                                freq += fweight * betaCovariate[cov][level][whichfnt] / nsib * (family.sibship ? 2 : 1);
                                                            }
                                                        }
                                                        freq += mweight * betaCovariate[cov][level][whichmt] / nsib;
                                                        if (!family.sibship) {
                                                            freq += mweight * betaCovariate[cov][level][whichmnt] / nsib * (family.sibship ? 2 : 1);
                                                        }
                                                    } else {
                                                        if (!options.genotype && !MchrX) {
                                                            linear[sib] += fweight * betaparentCovariate[cov][level][whichft];
                                                        }
                                                        linear[sib] += mweight * betaparentCovariate[cov][level][whichmt];
                                                    }
                                                }
                                            }
                                        } // for sib

                                        // total contribution from all sibs
                                        double beta1 = 0;
                                        for (int sib = 0; sib < nsib; sib++) {
                                            beta1 += invVar * (sibTrait[sib] - parentSumIntercept) * linear[sib];
                                            if (options.hhrr) {
                                                beta1 -= invVar * parentSumIntercept * alphaterm[sib];
                                            }
                                            if (normal) {
                                                beta1 -= invVar * 0.5 * linear[sib] * linear[sib];
                                            }
                                            if (invCov != 0) {
                                                for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                        beta1 += invCov * (sibTrait[sib2] - parentSumIntercept) * linear[sib];
                                                        if (options.hhrr) {
                                                            beta1 -= invCov * parentSumIntercept * alphaterm[sib];
                                                        }
                                                        if (normal) {
                                                            beta1 -= invCov * 0.5 * linear[sib] * linear[sib2];
                                                        }
                                                    } // for sib2
                                            } // if invCov!=0
                                        } // for sib

                                        if (positiveProb) {
                                            denom.push_back(freq + beta1);
                                        }
                                        positiveProbQueue.push_back(positiveProb);
                                    } // for phase


                            } // for ft, fnt, mt, mnt

            denomSum = 0;
            if (denom.size() > 0) {
                double maxbeta1 = *(max_element(denom.begin(), denom.end()));
                for (int j = 0; j < denom.size(); j++) {
                    denomSum += exp(denom[j] - maxbeta1);
                }
                denomSum = log(denomSum) + maxbeta1;
                if (options.llhd) {
                    cout << "denomSum " << denomSum << endl;
                }
            }

        } // options.slow



    } // !useSavedDenom

    // now the numerator

    //  numerator

    vector<double> numer;
    vector<double> sibNumer;

    // if allowing for linkage, loop through all combinations in all sibs
    for (int posn = nsib + 1; posn < consistentHaps[nfamily].size(); posn += 3 + nsib * 8) {

        // check for zero frequency haplotypes
        int sumZero = 0;
        for (int sib = 0; sib < nsib; sib++) {
            for (int i = 0; i < 8; i++) {
                int hap = consistentHaps[nfamily][posn+3+sib*8+i];
                if (!(options.genotype && i < 4)) {
                    sumZero += zero[hap];
                }
            }
        }

        if (!sumZero) {
            // virtual controls
            vector<vector<double> > virtualc(nsib);
            vector<vector<double> > parentvirtualc(nsib);
            for (int i = 0; i < nsib; i++) {
                virtualc[i].resize(4, 0);
                parentvirtualc[i].resize(4, 0);
            }

            // components of the linear models
            vector<vector<double> > freq(nsib);
            vector<vector<double> > linear(nsib);
            vector<vector<double> > alphaterm(nsib);
            vector<vector<double> > parentlinear(nsib);
            vector<vector<double> > parentalpha(nsib);
            for (int i = 0; i < nsib; i++) {
                freq[i].resize(4, 0);
                linear[i].resize(4, 0);
                alphaterm[i].resize(4, 0);
                parentlinear[i].resize(4, 0);
                parentalpha[i].resize(4, 0);
            }

            // note whether both haplotypes are known in the parents
            // relevant for the sibship test
            int knowFather = consistentHaps[nfamily][posn+1];
            int knowMother = consistentHaps[nfamily][posn+2];

            for (int sib = 0; sib < nsib; sib++) {
                int theHaps[2][2][2];
                int whichHaps[2][2][2];
                int ix = posn + 3 + sib * 8;
                for (int i = 0; i < 2; i++) // father, mother
                    for (int j = 0; j < 2; j++)
                        for (int k = 0; k < 2; k++) {
                            theHaps[i][j][k] = consistentHaps[nfamily][ix++];
                        }
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++)
                        for (int k = 0; k < 2; k++)
                            if (permute[nfamily]) {
                                whichHaps[i][j][k] = theHaps[i][j][k];
                            } else {
                                whichHaps[i][j][k] = theHaps[i][!j][!k];
                            }

                bool MchrX = options.chrX && sibSex[sib] == MALE;

                // frequency
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++) {
                        if (!options.genotype) {
                            if (!options.chrX) {
                                freq[sib][i*2+j] += frequency[whichHaps[0][i][j]] + (!family.sibship) * frequency[whichHaps[0][!i][!j]];
                            } else if (!MchrX) {
                                freq[sib][i*2+j] += frequency[whichHaps[0][0][j]];
                            }
                        }
                        freq[sib][i*2+j] += frequency[whichHaps[1][i][j]] + (!family.sibship) * frequency[whichHaps[1][!i][!j]];
                    }

                // haplotype effects
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++) {
                        if (typeOfPhenotype == "polytomous") {
                        	if (sibTrait[sib] < K) {
                        	    if (!options.genotype && !MchrX) { // On compte l'haplotype paternel
	                        		linear[sib][i*2+j] += beta[whichHaps[0][i][j] + sibTrait[sib]*nhap];
									parentlinear[sib][i*2+j] += betaparent[whichHaps[0][i][j] + sibTrait[sib]*nhap];
									}
								// ensuite on compte l'haplotype maternel
		                        linear[sib][i*2+j] += beta[whichHaps[1][i][j] + sibTrait[sib]*nhap];
		                        parentlinear[sib][i*2+j] += betaparent[whichHaps[1][i][j] + sibTrait[sib]*nhap];
		                        }
							}
                        else { // Phenotype is not polytomous
                        	if (!options.genotype && !MchrX) {                        
	                            linear[sib][i*2+j] += beta[whichHaps[0][i][j]];
	                            alphaterm[sib][i*2+j] += alpha[whichHaps[0][i][j]];
								parentlinear[sib][i*2+j] += betaparent[whichHaps[0][i][j]];
	                        }
	                        linear[sib][i*2+j] += beta[whichHaps[1][i][j]];
	                        alphaterm[sib][i*2+j] += alpha[whichHaps[1][i][j]];
	                        parentlinear[sib][i*2+j] += betaparent[whichHaps[1][i][j]];
	                        }
                    }

                // covariate effects
                for (int cov = 0; cov < covariateName.size(); cov++) {
                    int level = 0;
                    vector<double> weight(2, 0);
                    if (factor[cov] && sibCovariate[sib][cov] != 0) {
                        level = (int)sibCovariate[sib][cov] - 1;
                        weight[0] = 1.0;
                        weight[1] = 1.0;
                    } else {
                        weight[0] = sibCovariate[sib][cov];
                        weight[1] = sibCovariate[sib][cov];
                        if (!family.sibship && covariateName[cov] == "parsex") {
                            weight[0] = (baseline[cov] != MALE);
                            weight[1] = (baseline[cov] != FEMALE);
                        }
                    }
                    double covTerm[2][2][2];
                    double parentCovTerm[2][2][2];
                    for (int i = 0; i < 2; i++)
                        for (int j = 0; j < 2; j++)
                            for (int k = 0; k < 2; k++) {
                                covTerm[i][j][k] = weight[i] * betaCovariate[cov][level][whichHaps[i][j][k]];
                                parentCovTerm[i][j][k] = weight[i] * betaparentCovariate[cov][level][whichHaps[i][j][k]];
                            }

                    if (confounder[cov]) {
                        for (int i = 0; i < 2; i++)
                            for (int j = 0; j < 2; j++) {
                                if (!options.genotype) {
                                    if (!options.chrX) {
                                        freq[sib][i*2+j] += covTerm[0][i][j] + (!family.sibship) * covTerm[0][!i][!j];
                                    } else if (!MchrX) {
                                        freq[sib][i*2+j] += covTerm[0][0][j];
                                    }
                                }
                                freq[sib][i*2+j] += covTerm[1][i][j] + (!family.sibship) * covTerm[1][!i][!j];
                            }
                    } else {
                        for (int i = 0; i < 2; i++)
                            for (int j = 0; j < 2; j++) {
                                if (!options.genotype && !MchrX) {
                                    linear[sib][i*2+j] += covTerm[0][i][j];
                                    parentlinear[sib][i*2+j] += parentCovTerm[0][i][j];
                                }
                                linear[sib][i*2+j] += covTerm[1][i][j];
                                parentlinear[sib][i*2+j] += parentCovTerm[1][i][j];
                            }
                    }
                } // for covariates

            } // for sib

            // the virtual controls
            for (int sib = 0; sib < nsib; sib++)
                for (int v = 0; v < 4; v++) {
                	if (typeOfPhenotype == "polytomous") {
                		virtualc[sib][v] += linear[sib][v];
                		parentvirtualc[sib][v] += parentlinear[sib][v];
                		}
                    else {
                    	virtualc[sib][v] += invVar * (sibTrait[sib] - sumIntercept) * linear[sib][v];
                    	virtualc[sib][v] -= invVar * sumIntercept * alphaterm[sib][v];
                    	parentvirtualc[sib][v] += invVar * (sibTrait[sib] - parentSumIntercept) * parentlinear[sib][v];
                    }
                    if (normal) {
                        virtualc[sib][v] -= invVar * 0.5 * linear[sib][v] * linear[sib][v];
                        parentvirtualc[sib][v] -= invVar * 0.5 * parentlinear[sib][v] * parentlinear[sib][v];
                    }
                    if (invCov != 0)
                        for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                virtualc[sib][v] += invCov * (sibTrait[sib2] - sumIntercept) * linear[sib][v];
                                virtualc[sib][v] -= invCov * sumIntercept * alphaterm[sib][v];
                                parentvirtualc[sib][v] += invCov * (sibTrait[sib2] - parentSumIntercept) * parentlinear[sib][v];
                                if (normal) {
                                    virtualc[sib][v] -= invCov * 0.5 * linear[sib][v] * linear[sib2][v];
                                    parentvirtualc[sib][v] -= invCov * 0.5 * parentlinear[sib][v] * parentlinear[sib2][v];
                                }
                            } // for sib2
                    // cancel out virtual controls in HHRR analysis
                    if (options.hhrr) {
                        parentvirtualc[sib][v] = virtualc[sib][v];
                    }
                } // for v, sib

            // number of realisations consistent with missing haplotypes in parents
            double knowParents = 1 + !knowFather * (nhapz - 1 + (1 << nsib) / 2 - 1);
            knowParents *= 1 + !knowMother * (nhapz - 1 + (1 << nsib) / 2 - 1);
            knowParents = log(knowParents);

            // number of equivalent phases from homozygous parents
            double nhomphase = log((double)consistentHaps[nfamily][posn]);

            // total contribution for the family
            double beta0 = 0;

            // if no linkage, do sibs independently
            if (options.nolinkage) {
                // loglikelihood omitting one sib
                vector<double> betaDropOne(nsib, 0);
                // contributions of each sib and virtual control
                vector<vector<double> > sibContrib(nsib);
                for (int i = 0; i < nsib; i++) {
                    sibContrib[i].resize(4, 0);
                }

                for (int sib = 0; sib < nsib; sib++) {
                    vector<double> beta1(4, 0);
                    for (int i = 0; i < 4; i++) {
                        beta1[i] = freq[sib][i] / nsib * (family.sibship ? 2 : 1) + virtualc[sib][0];
                    }
                    double maxbeta = max(virtualc[sib][0], virtualc[sib][3]);
                    if (!twoVirtualControls) {
                        maxbeta = max(maxbeta, virtualc[sib][1]);
                        maxbeta = max(maxbeta, virtualc[sib][2]);
                    }
                    for (int i = 0; i < 4; i += (twoVirtualControls ? 3 : 1)) {
                        double x = exp(virtualc[sib][0] - maxbeta) + exp(virtualc[sib][3] - maxbeta);
                        if (!twoVirtualControls) {
                            x += exp(virtualc[sib][1] - maxbeta) + exp(virtualc[sib][2] - maxbeta);
                        }
                        beta1[i] -= log(x) + maxbeta;
                        sibContrib[sib][i] = beta1[i] + parentvirtualc[sib][i];
                    }
                    maxbeta = max(sibContrib[sib][0], sibContrib[sib][3]);
                    if (!twoVirtualControls) {
                        maxbeta = max(maxbeta, sibContrib[sib][1]);
                        maxbeta = max(maxbeta, sibContrib[sib][2]);
                    }
                    double x = exp(sibContrib[sib][0] - maxbeta) + exp(sibContrib[sib][3] - maxbeta);
                    if (!twoVirtualControls) {
                        x += exp(sibContrib[sib][1] - maxbeta) + exp(sibContrib[sib][2] - maxbeta);
                    }
                    beta0 += log(x) + maxbeta;
                    for (int sib2 = 0; sib2 < nsib; sib2++)
                        if (sib2 != sib) {
                            betaDropOne[sib2] += log(x) + maxbeta;
                        }

                } // for sib

                for (int sib = 0; sib < nsib; sib++)
                    for (int v = 0; v < 4; v += (twoVirtualControls ? 3 : 1)) {
                        sibNumer.push_back(sibContrib[sib][v] + betaDropOne[sib] + knowParents + nhomphase - denomSum);
                    }

            } // if nolinkage

            else { // linkage, do sibs together
                vector<double> sumFreq(4, 0);
                vector<double> sumVirtualc(4, 0);
                vector<double> parentSumVirtualc(4, 0);
                for (int sib = 0; sib < nsib; sib++)
                    for (int v = 0; v < 4; v++) {
                        sumFreq[v] += freq[sib][v];
                        sumVirtualc[v] += virtualc[sib][v];
                        parentSumVirtualc[v] += parentvirtualc[sib][v];
                    }
                vector<double> beta1(4, 0);
                vector<double> contrib(4, 0);
                for (int i = 0; i < 4; i++) {
                    beta1[i] = sumFreq[i] / nsib * (family.sibship ? 2 : 1) + sumVirtualc[0];
                }
                double maxbeta = max(sumVirtualc[0], sumVirtualc[3]);
                if (!twoVirtualControls) {
                    maxbeta = max(maxbeta, sumVirtualc[1]);
                    maxbeta = max(maxbeta, sumVirtualc[2]);
                }
                double x = exp(sumVirtualc[0] - maxbeta) + exp(sumVirtualc[3] - maxbeta);
                if (!twoVirtualControls) {
                    x += exp(sumVirtualc[1] - maxbeta) + exp(sumVirtualc[2] - maxbeta);
                }
                for (int i = 0; i < 4; i += (twoVirtualControls ? 3 : 1)) {
                    beta1[i] -= log(x) + maxbeta;
                    contrib[i] = beta1[i] + parentSumVirtualc[i];
                }
                maxbeta = max(contrib[0], contrib[3]);
                if (!twoVirtualControls) {
                    maxbeta = max(maxbeta, contrib[1]);
                    maxbeta = max(maxbeta, contrib[2]);
                }
                x = exp(contrib[0] - maxbeta) + exp(contrib[3] - maxbeta);
                if (!twoVirtualControls) {
                    x += exp(contrib[1] - maxbeta) + exp(contrib[2] - maxbeta);
                }
                beta0 = log(x) + maxbeta;
                for (int v = 0; v < 4; v += (twoVirtualControls ? 3 : 1)) {
                    sibNumer.push_back(contrib[v] + knowParents + nhomphase - denomSum);
                }
            }
            // total family likelihood contribution for this completion
            numer.push_back(beta0 + knowParents + nhomphase - denomSum);


        } // if positiveProb
    } // for posn



    // form the numerator
    double numerSum = 0;
    if (numer.size() > 0) {
        double maxnumer = *(max_element(numer.begin(), numer.end()));
        for (int i = 0; i < numer.size(); i++) {
            if (options.llhd) {
                cout << numer[i] << " ";
            }
            numerSum += exp(numer[i] - maxnumer);
        }
        if (options.llhd) {
            cout << endl;
        }

        // the log-likelihood contribution
        numerSum = log(numerSum) + maxnumer;
        prob += numerSum;
    }

    // if (!numer.size()) return;

    // gradients
    vector<double> freqGradient(nhap, 0);
    vector<double> betaGradient(nhap*(K-1), 0);
    vector<double> betaparentGradient(nhap*(K-1), 0);
    vector<double> betaparent0Gradient(K-1, 0);
    vector<double> alphaGradient(nhap, 0);
    double alpha0Gradient = 0;
    vector<vector<vector<double> > > betaCovariateGradient(betaCovariate.size());
    vector<vector<vector<double> > > betaparentCovariateGradient(betaparentCovariate.size());
    vector<vector<double> > betaCovariate0Gradient(betaCovariate0.size());
    vector<vector<double> > betaparentCovariate0Gradient(betaparentCovariate0.size());
    for (int i = 0; i < betaCovariate.size(); i++) {
        betaCovariateGradient[i].resize(betaCovariate[i].size());
        for (int j = 0; j < betaCovariate[i].size(); j++) {
            betaCovariateGradient[i][j].resize(nhap, 0);
        }
    }
    for (int i = 0; i < betaparentCovariate.size(); i++) {
        betaparentCovariateGradient[i].resize(betaparentCovariate[i].size());
        for (int j = 0; j < betaparentCovariate[i].size(); j++) {
            betaparentCovariateGradient[i][j].resize(nhap, 0);
        }
    }
    for (int i = 0; i < betaCovariate0.size(); i++) {
        betaCovariate0Gradient[i].resize(betaCovariate0[i].size(), 0);
    }
    for (int i = 0; i < betaparentCovariate0.size(); i++) {
        betaparentCovariate0Gradient[i].resize(betaparentCovariate0[i].size(), 0);
    }

    // contributions to gradient from denominator
    for (int i = 0; i < 2; i++)
        if (haveGlobal[i] && useSavedDenom[i]) {
            gradient = globalGrad[i];
        }

    if (!(haveGlobal[0] && useSavedDenom[0]) &&
            !(haveGlobal[1] && useSavedDenom[1])) {

        //  fast algorithm

        // fast algorithm
        if (!normal && !options.slow &&
                !options.condgenotype && !options.genotype) {

            for (int i = 0; i < nhap; i++) if (!zero[i]) {
                    // loop through parental phases
                    for (int fphase = 0; fphase < nphase; fphase++)
                        for (int mphase = 0; mphase < nphase; mphase++) {
                            double xf = sumft[fphase] + (!options.chrX) * sumfnt[fphase] - denomSum;
                            double xm = summt[mphase] + summnt[mphase] - denomSum;
                            for (int sib = 0; sib < nsib; sib++) {
                                bool MchrX = options.chrX && sibSex[sib] == MALE;
                                double derivTermft = 0, derivTermfnt = 0, derivTermmt = 0, derivTermmnt = 0;
                                if ((fphase & 1 << sib) == 0) {
                                    derivTermft = exp(dft[fphase][i] + (!options.chrX) * sumfnt[fphase] + xm);
                                    derivTermfnt = exp(dfnt[fphase][i] + sumft[fphase] + xm);
                                } else {
                                    derivTermft = exp(dfnt[fphase][i] + (!options.chrX) * sumft[fphase] + xm);
                                    derivTermfnt = exp(dft[fphase][i] + sumfnt[fphase] + xm);
                                }
                                if ((mphase & 1 << sib) == 0) {
                                    derivTermmt = exp(dmt[mphase][i] + summnt[mphase] + xf);
                                    derivTermmnt = exp(dmnt[mphase][i] + summt[mphase] + xf);
                                } else {
                                    derivTermmt = exp(dmnt[mphase][i] + summt[mphase] + xf);
                                    derivTermmnt = exp(dmt[mphase][i] + summnt[mphase] + xf);
                                }

                                // paternal transmissions
                                if (!MchrX) {
                                    freqGradient[i] += derivTermft / nsib * (family.sibship ? 2 : 1);
                                    if (typeOfPhenotype == "polytomous") {
                                    	if (sibTrait[sib] < K) betaparentGradient[i + sibTrait[sib]*nhap] += derivTermft;
                                    	}
                                    else betaparentGradient[i] += invVar * (sibTrait[sib] - parentSumIntercept) * derivTermft;
                                    if (options.hhrr) {
                                        alphaGradient[i] -= invVar * parentSumIntercept * derivTermft;
                                    }
                                    if (invCov != 0)
                                        for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                betaparentGradient[i] += invCov * (sibTrait[sib2] - parentSumIntercept) * derivTermft;
                                                if (options.hhrr) {
                                                    alphaGradient[i] -= invCov * parentSumIntercept * derivTermft;
                                                }

                                            }
                                    if (typeOfPhenotype == "polytomous") {
                                    	if (sibTrait[sib] < K) betaparent0Gradient[sibTrait[sib]] -= betaparent[i + sibTrait[sib]*nhap] * derivTermft;
                                    	}
                                    else betaparent0Gradient[0] -= (invVar + (nsib - 1) * invCov) * betaparent[i] * derivTermft;
                                    if (options.hhrr) {
                                        betaparent0Gradient[0] -= (invVar + (nsib - 1) * invCov) * alpha[i] * derivTermft;
                                    }
                                    for (int cov = 0; cov < betaCovariate.size(); cov++)
                                        for (int level = 0; level < betaCovariate[cov].size(); level++) {
                                            betaparentCovariate0Gradient[cov][level] -= (invVar + (nsib - 1) * invCov) * betaparent[i] * derivTermft;
                                            if (options.hhrr) {
                                                betaparentCovariate0Gradient[cov][level] -= (invVar + (nsib - 1) * invCov) * alpha[i] * derivTermft;
                                            }
                                        }
                                }
                                if (!options.chrX && !family.sibship) {
                                    freqGradient[i] += derivTermfnt / nsib * (family.sibship ? 2 : 1);
                                }

                                // maternal transmissions
                                freqGradient[i] += derivTermmt / nsib * (family.sibship ? 2 : 1);
                                if (!family.sibship) {
                                    freqGradient[i] += derivTermmnt / nsib * (family.sibship ? 2 : 1);
                                }
                                if (typeOfPhenotype == "polytomous") {
                                    if (sibTrait[sib] < K) betaparentGradient[i + sibTrait[sib]*nhap] += derivTermmt;
                                    }
                                else betaparentGradient[i] += invVar * (sibTrait[sib] - parentSumIntercept) * derivTermmt;
                                if (options.hhrr) {
                                    alphaGradient[i] -= invVar * parentSumIntercept * derivTermmt;
                                }
                                if (invCov != 0)
                                    for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                            betaparentGradient[i] += invCov * (sibTrait[sib2] - parentSumIntercept) * derivTermmt;
                                            if (options.hhrr) {
                                                alphaGradient[i] -= invCov * parentSumIntercept * derivTermmt;
                                            }
                                        }
                                if (typeOfPhenotype == "polytomous") {
                                    if (sibTrait[sib] < K) betaparent0Gradient[sibTrait[sib]] -= betaparent[i + sibTrait[sib]*nhap] * derivTermmt;
                                    }
                                else betaparent0Gradient[0] -= (invVar + (nsib - 1) * invCov) * betaparent[i] * derivTermmt;
                                if (options.hhrr) {
                                    betaparent0Gradient[0] -= (invVar + (nsib - 1) * invCov) * alpha[i] * derivTermmt;
                                }
                                for (int cov = 0; cov < betaCovariate.size(); cov++)
                                    for (int level = 0; level < betaCovariate[cov].size(); level++) {
                                        betaparentCovariate0Gradient[cov][level] -= (invVar + (nsib - 1) * invCov) * betaparent[i] * derivTermmt;
                                        if (options.hhrr) {
                                            betaparentCovariate0Gradient[cov][level] -= (invVar + (nsib - 1) * invCov) * alpha[i] * derivTermmt;
                                        }
                                    }

                                // covariates
                                for (int cov = 0; cov < covariateName.size(); cov++) {
                                    int level = 0;
                                    double fweight = 0, mweight = 0;
                                    if (factor[cov] && sibCovariate[sib][cov] != 0) {
                                        level = (int)sibCovariate[sib][cov] - 1;
                                        fweight = 1.0;
                                        mweight = 1.0;
                                    } else {
                                        fweight = sibCovariate[sib][cov];
                                        mweight = sibCovariate[sib][cov];
                                        if (covariateName[cov] == "parsex") {
                                            fweight = (baseline[cov] != MALE);
                                            mweight = (baseline[cov] != FEMALE);
                                        }
                                    }
                                    if (!MchrX) {
                                        if (confounder[cov]) {
                                            betaCovariateGradient[cov][level][i] += fweight * derivTermft / nsib * (family.sibship ? 2 : 1);
                                            if (!options.chrX && !family.sibship) {
                                                betaCovariateGradient[cov][level][i] += fweight * derivTermfnt / nsib * (family.sibship ? 2 : 1);
                                            }
                                        } else {
                                            betaparentCovariateGradient[cov][level][i] += invVar * (sibTrait[sib] - parentSumIntercept) * fweight * derivTermft;
                                            if (invCov != 0)
                                                for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                        betaparentCovariateGradient[cov][level][i] += invCov * (sibTrait[sib2] - parentSumIntercept) * fweight * derivTermft;
                                                    }
                                            betaparent0Gradient[0] -= (invVar + (nsib - 1) * invCov) * betaparentCovariate[cov][level][i] * fweight * derivTermft;
                                            for (int cov2 = 0; cov2 < betaCovariate.size(); cov2++)
                                                for (int level2 = 0; level2 < betaCovariate[cov2].size(); level2++) {
                                                    betaparentCovariate0Gradient[cov2][level2] -= (invVar + (nsib - 1) * invCov) * betaparentCovariate[cov][level][i] * fweight * derivTermft;
                                                }
                                        }
                                    }
                                    if (confounder[cov]) {
                                        betaCovariateGradient[cov][level][i] += mweight * derivTermmt / nsib * (family.sibship ? 2 : 1);
                                        if (!family.sibship) {
                                            betaCovariateGradient[cov][level][i] += mweight * derivTermmnt / nsib * (family.sibship ? 2 : 1);
                                        }
                                    } else {
                                        betaparentCovariateGradient[cov][level][i] += invVar * (sibTrait[sib] - parentSumIntercept) * mweight * derivTermmt;
                                        if (invCov != 0)
                                            for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                    betaparentCovariateGradient[cov][level][i] += invCov * (sibTrait[sib2] - parentSumIntercept) * mweight * derivTermmt;
                                                }
                                        betaparent0Gradient[0] -= (invVar + (nsib - 1) * invCov) * betaparentCovariate[cov][level][i] * mweight * derivTermmt;
                                        for (int cov2 = 0; cov2 < betaCovariate.size(); cov2++)
                                            for (int level2 = 0; level2 < betaCovariate[cov2].size(); level2++) {
                                                betaparentCovariate0Gradient[cov2][level2] -= (invVar + (nsib - 1) * invCov) * betaparentCovariate[cov][level][i] * mweight * derivTermmt;
                                            }
                                    }

                                } // for cov
                            } // for sib
                        } // for mphase
                } // for i

        } // if (!options.slow)



        //  slow version

        if (normal || options.slow ||
                options.condgenotype || options.genotype)  {
            // slow version for development
            int ix = 0;
            for (int mt = 0; mt < haploCode.size(); mt++)
                for (int mnt = 0; mnt < haploCode.size(); mnt++)
                    for (int ft = 0; ft < haploCode.size(); ft++)
                        for (int fnt = 0; fnt < haploCode.size(); fnt++) if (!options.chrX || ft == fnt) {

                                for (int i = 0; i < nphase; i++)
                                    for (int j = 0; j < nphase; j++) {
                                        if (positiveProbQueue[0]) {
                                            double x = exp(denom[ix++] - denomSum);

                                            vector<double> linear(nsib, 0);

                                            for (int sib = 0; sib < nsib; sib++) {
                                                bool fphase = (i & 1 << sib) == 0;
                                                bool mphase = (j & 1 << sib) == 0;

                                                int whichft = ft;
                                                int whichfnt = fnt;
                                                int whichmt = mt;
                                                int whichmnt = mnt;
                                                if (!fphase) {
                                                    whichft = fnt;
                                                    whichfnt = ft;
                                                }
                                                if (!mphase) {
                                                    whichmt = mnt;
                                                    whichmnt = mt;
                                                }

                                                bool MchrX = options.chrX && sibSex[sib] == MALE;

                                                // get the genotype code for this child
                                                if (options.condgenotype || options.genotype) {
                                                    Haplotype fthap = haploCode[whichft];
                                                    Haplotype fnthap = haploCode[whichfnt];
                                                    Haplotype mthap = haploCode[whichmt];
                                                    Haplotype mnthap = haploCode[whichmnt];
                                                    int ncondition = options.condgenotype * options.condition.size();
                                                    Haplotype ftgeno = fthap.formgeno(mthap, ncondition, options.genotype);
                                                    Haplotype fntgeno = fnthap.formgeno(mnthap, ncondition, options.genotype);
                                                    Haplotype mtgeno = mthap.formgeno(fthap, ncondition, options.genotype);
                                                    Haplotype mntgeno = mnthap.formgeno(fnthap, ncondition, options.genotype);
                                                    if (options.chrX && sibSex[sib] == MALE) {
                                                        ftgeno = mthap.formgeno(mthap, ncondition, options.genotype);
                                                        fntgeno = mnthap.formgeno(mnthap, ncondition, options.genotype);
                                                        mtgeno = mthap.formgeno(mthap, ncondition, options.genotype);
                                                        mntgeno = mnthap.formgeno(mnthap, ncondition, options.genotype);
                                                    }
                                                    whichft = ftgeno.index(genoCode);
                                                    whichfnt = fntgeno.index(genoCode);
                                                    whichmt = mtgeno.index(genoCode);
                                                    whichmnt = mntgeno.index(genoCode);
                                                }

                                                if (!options.genotype) {
                                                    if (!options.chrX) {
                                                        freqGradient[whichft] += x / nsib * (family.sibship ? 2 : 1);
                                                    }
                                                    if (!MchrX && !family.sibship) {
                                                        freqGradient[whichfnt] += x / nsib * (family.sibship ? 2 : 1);
                                                    }
                                                }
                                                freqGradient[whichmt] += x / nsib * (family.sibship ? 2 : 1);
                                                if (!family.sibship) {
                                                    freqGradient[whichmnt] += x / nsib * (family.sibship ? 2 : 1);
                                                }

                                                // paternal transmission
                                                if (!options.genotype && !MchrX) {
                                                    betaparentGradient[whichft] += invVar * (sibTrait[sib] - parentSumIntercept) * x;
                                                    if (options.hhrr) {
                                                        alphaGradient[whichft] -= invVar * parentSumIntercept * x;
                                                    }
                                                    if (invCov != 0)
                                                        for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                                betaparentGradient[whichft] += invCov * (sibTrait[sib2] - parentSumIntercept) * x;
                                                                if (options.hhrr) {
                                                                    alphaGradient[whichft] -= invCov * parentSumIntercept * x;
                                                                }
                                                            }
                                                    linear[sib] += betaparent[whichft];
                                                    betaparent0Gradient[0] -= (invVar + (nsib - 1) * invCov) * betaparent[whichft] * x;
                                                    if (options.hhrr) {
                                                        betaparent0Gradient[0] -= (invVar + (nsib - 1) * invCov) * alpha[whichft] * x;
                                                    }
                                                    for (int cov = 0; cov < betaCovariate.size(); cov++)
                                                        for (int level = 0; level < betaCovariate[cov].size(); level++) {
                                                            betaparentCovariate0Gradient[cov][level] -= (invVar + (nsib - 1) * invCov) * betaparent[whichft] * x;
                                                            if (options.hhrr) {
                                                                betaparentCovariate0Gradient[cov][level] -= (invVar + (nsib - 1) * invCov) * alpha[whichft] * x;
                                                            }
                                                        }
                                                }
                                                // maternal transmission
                                                betaparentGradient[whichmt] += invVar * (sibTrait[sib] - parentSumIntercept) * x;
                                                if (options.hhrr) {
                                                    alphaGradient[whichmt] -= invVar * parentSumIntercept * x;
                                                }
                                                if (invCov != 0)
                                                    for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                            betaparentGradient[whichmt] += invCov * (sibTrait[sib2] - parentSumIntercept) * x;
                                                            if (options.hhrr) {
                                                                alphaGradient[whichmt] -= invCov * parentSumIntercept * x;
                                                            }
                                                        }
                                                linear[sib] += betaparent[whichmt];
                                                betaparent0Gradient[0] -= (invVar + (nsib - 1) * invCov) * betaparent[whichmt] * x;
                                                if (options.hhrr) {
                                                    betaparent0Gradient[0] -= (invVar + (nsib - 1) * invCov) * alpha[whichmt] * x;
                                                }
                                                for (int cov = 0; cov < betaCovariate.size(); cov++)
                                                    for (int level = 0; level < betaCovariate[cov].size(); level++) {
                                                        betaparentCovariate0Gradient[cov][level] -= (invVar + (nsib - 1) * invCov) * betaparent[whichmt] * x;
                                                        if (options.hhrr) {
                                                            betaparentCovariate0Gradient[cov][level] -= (invVar + (nsib - 1) * invCov) * alpha[whichmt] * x;
                                                        }
                                                    }

                                                // covariates
                                                for (int cov = 0; cov < covariateName.size(); cov++) {
                                                    double fweight = 0, mweight = 0;
                                                    int level = 0;
                                                    if (factor[cov] && sibCovariate[sib][cov] != 0) {
                                                        level = (int)sibCovariate[sib][cov] - 1;
                                                        fweight = 1.0;
                                                        mweight = 1.0;
                                                    } else {
                                                        fweight = sibCovariate[sib][cov];
                                                        mweight = sibCovariate[sib][cov];
                                                        if (covariateName[cov] == "parsex") {
                                                            fweight = (baseline[cov] != MALE);
                                                            mweight = (baseline[cov] != FEMALE);
                                                        }
                                                    }

                                                    if (confounder[cov]) {
                                                        if (!options.genotype) {
                                                            if (!options.chrX) {
                                                                betaCovariateGradient[cov][level][whichft] += fweight * x / nsib * (family.sibship ? 2 : 1);
                                                            }
                                                            if (!MchrX && !family.sibship) {
                                                                betaCovariateGradient[cov][level][whichfnt] += fweight * x / nsib * (family.sibship ? 2 : 1);
                                                            }
                                                        }
                                                        betaCovariateGradient[cov][level][whichmt] += mweight * x / nsib * (family.sibship ? 2 : 1);
                                                        if (!family.sibship) {
                                                            betaCovariateGradient[cov][level][whichmnt] += mweight * x / nsib * (family.sibship ? 2 : 1);
                                                        }
                                                    } else {
                                                        if (!options.genotype && !MchrX) {
                                                            betaparentCovariateGradient[cov][level][whichft] += invVar * (sibTrait[sib] - parentSumIntercept) * fweight * x;
                                                            if (invCov != 0)
                                                                for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                                        betaparentCovariateGradient[cov][level][whichft] += invCov * (sibTrait[sib2] - parentSumIntercept) * fweight * x;
                                                                    }
                                                            linear[sib] += fweight * betaparentCovariate[cov][level][whichft];
                                                            betaparent0Gradient[0] -= (invVar + (nsib - 1) * invCov) * betaparentCovariate[cov][level][whichft] * fweight * x;
                                                            for (int cov2 = 0; cov2 < betaCovariate.size(); cov2++)
                                                                for (int level2 = 0; level2 < betaCovariate[cov2].size(); level2++) {
                                                                    betaparentCovariate0Gradient[cov2][level2] -= (invVar + (nsib - 1) * invCov) * betaparentCovariate[cov][level][whichft] * fweight * x;
                                                                }
                                                        }
                                                        betaparentCovariateGradient[cov][level][whichmt] += invVar * (sibTrait[sib] - parentSumIntercept) * mweight * x;
                                                        if (invCov != 0)
                                                            for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                                    betaparentCovariateGradient[cov][level][whichmt] += invCov * (sibTrait[sib2] - parentSumIntercept) * mweight * x;
                                                                }
                                                        linear[sib] += mweight * betaparentCovariate[cov][level][whichmt];
                                                        betaparent0Gradient[0] -= (invVar + (nsib - 1) * invCov) * betaparentCovariate[cov][level][whichmt] * mweight * x;
                                                        for (int cov2 = 0; cov2 < betaCovariate.size(); cov2++)
                                                            for (int level2 = 0; level2 < betaCovariate[cov2].size(); level2++) {
                                                                betaparentCovariate0Gradient[cov2][level2] -= (invVar + (nsib - 1) * invCov) * betaparentCovariate[cov][level][whichmt] * mweight * x;
                                                            }
                                                    }
                                                } // for cov
                                            } // for sib

                                            // additions for normal model
                                            if (normal) {
                                                for (int sib = 0; sib < nsib; sib++) {
                                                    bool fphase = (i & 1 << sib) == 0;
                                                    bool mphase = (j & 1 << sib) == 0;

                                                    int whichft = ft;
                                                    int whichfnt = fnt;
                                                    int whichmt = mt;
                                                    int whichmnt = mnt;
                                                    if (!fphase) {
                                                        whichft = fnt;
                                                        whichfnt = ft;
                                                    }
                                                    if (!mphase) {
                                                        whichmt = mnt;
                                                        whichmnt = mt;
                                                    }

                                                    bool MchrX = options.chrX && sibSex[sib] == MALE;

                                                    // get the genotype code for this child
                                                    if (options.condgenotype || options.genotype) {
                                                        Haplotype fthap = haploCode[whichft];
                                                        Haplotype fnthap = haploCode[whichfnt];
                                                        Haplotype mthap = haploCode[whichmt];
                                                        Haplotype mnthap = haploCode[whichmnt];
                                                        int ncondition = options.condgenotype * options.condition.size();
                                                        Haplotype ftgeno = fthap.formgeno(mthap, ncondition, options.genotype);
                                                        Haplotype fntgeno = fnthap.formgeno(mnthap, ncondition, options.genotype);
                                                        Haplotype mtgeno = mthap.formgeno(fthap, ncondition, options.genotype);
                                                        Haplotype mntgeno = mnthap.formgeno(fnthap, ncondition, options.genotype);
                                                        if (options.chrX && sibSex[sib] == MALE) {
                                                            ftgeno = mthap.formgeno(mthap, ncondition, options.genotype);
                                                            fntgeno = mnthap.formgeno(mnthap, ncondition, options.genotype);
                                                            mtgeno = mthap.formgeno(mthap, ncondition, options.genotype);
                                                            mntgeno = mnthap.formgeno(mnthap, ncondition, options.genotype);
                                                        }
                                                        whichft = ftgeno.index(genoCode);
                                                        whichfnt = fntgeno.index(genoCode);
                                                        whichmt = mtgeno.index(genoCode);
                                                        whichmnt = mntgeno.index(genoCode);
                                                    }

                                                    if (!options.genotype && !MchrX) {
                                                        betaparentGradient[whichft] -= invVar * linear[sib] * x;
                                                    }
                                                    betaparentGradient[whichmt] -= invVar * linear[sib] * x;
                                                    if (invCov != 0)
                                                        for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                                if (!options.genotype && !MchrX) {
                                                                    betaparentGradient[whichft] -= invCov * linear[sib2] * x;
                                                                }
                                                                betaparentGradient[whichmt] -= invCov * linear[sib2] * x;
                                                            }
                                                    for (int cov = 0; cov < covariateName.size(); cov++) {
                                                        double fweight = 0, mweight = 0;
                                                        int level = 0;
                                                        if (factor[cov] && sibCovariate[sib][cov] != 0) {
                                                            level = (int)sibCovariate[sib][cov] - 1;
                                                            fweight = 1.0;
                                                            mweight = 1.0;
                                                        } else {
                                                            fweight = sibCovariate[sib][cov];
                                                            mweight = sibCovariate[sib][cov];
                                                            if (covariateName[cov] == "parsex") {
                                                                fweight = (baseline[cov] != MALE);
                                                                mweight = (baseline[cov] != FEMALE);
                                                            }
                                                        }
                                                        if (!confounder[cov]) {
                                                            if (!options.genotype && !MchrX) {
                                                                betaparentCovariateGradient[cov][level][whichft] -= invVar * linear[sib] * fweight * x;
                                                            }
                                                            betaparentCovariateGradient[cov][level][whichmt] -= invVar * linear[sib] * mweight * x;
                                                            if (invCov != 0)
                                                                for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                                        if (!options.genotype && !MchrX) {
                                                                            betaparentCovariateGradient[cov][level][whichft] -= invCov * linear[sib2] * fweight * x;
                                                                        }
                                                                        betaparentCovariateGradient[cov][level][whichmt] -= invCov * linear[sib2] * mweight * x;
                                                                    }
                                                        }
                                                    } // for cov
                                                } // for sib
                                            } // if normal

                                        } // if positiveProb
                                        positiveProbQueue.pop_front();
                                    } // for j

                            } // for fnt
        }



    } // if !useSavedDenom

    if (options.hhrr) {
        for (int i = 0; i < nhap; i++) {
            betaGradient[i] += betaparentGradient[i];
            for (int j = 0; j < betaCovariate.size(); j++)
                for (int k = 0; k < betaCovariate[j].size(); k++) {
                    betaCovariateGradient[j][k][i] += betaparentCovariateGradient[j][k][i];
                }
        }
        alpha0Gradient += betaparent0Gradient[0];
        for (int j = 0; j < betaCovariate.size(); j++)
            for (int k = 0; k < betaCovariate[j].size(); k++) {
                betaCovariate0Gradient[j][k] += betaparentCovariate0Gradient[j][k];
            }
    }

    // save global denominator
    // i = 0 pour un seul enfant atteint, i = 1 pour une paire d'enfants atteints
    for (int i = 0; i < 2; i++)
        if (!haveGlobal[i] && useSavedDenom[i]) {
            haveGlobal[i] = true;
            globalDenom[i] = denomSum;
            globalGrad[i].resize(gradient.size(), 0);
            int ix = 0;
            if (typeOfPhenotype == "polytomous") betasize = nhap*(K-1);
            else betasize = nhap;
            for (int j = 0; j < betasize; j++) {
                globalGrad[i][ix++] += betaGradient[j];
            }
            for (int j = 0; j < nhap; j++) {
                globalGrad[i][ix++] += freqGradient[j];
            }
            if (haveBetaParent(options)) {
                for (int j = 0; j < betasize; j++) {
                    globalGrad[i][ix++] += betaparentGradient[j];
                }
                if (haveBetaParent0(options)) {
                    globalGrad[i][ix++] += betaparent0Gradient[0];
                }
            }
            if (haveAlpha(options)) {
                for (int j = 0; j < nhap; j++) {
                    globalGrad[i][ix++] += alphaGradient[j];
                }
                if (haveAlpha0(options)) {
                    globalGrad[i][ix++] += alpha0Gradient;
                }
            }
            for (int j = 0; j < betaCovariate.size(); j++)
                for (int k = 0; k < betaCovariate[j].size(); k++) {
                    for (int l = 0; l < nhap; l++) {
                        globalGrad[i][ix++] += betaCovariateGradient[j][k][l];
                    }
                    if (!confounder[j] && haveFamilies && !options.hhrr)
                        for (int l = 0; l < nhap; l++) {
                            globalGrad[i][ix++] += betaparentCovariateGradient[j][k][l];
                        }
                    if (typeOfPhenotype == "quant") {
                        globalGrad[i][ix++] += betaCovariate0Gradient[j][k];
                        if (haveFamilies && !options.hhrr) {
                            globalGrad[i][ix++] += betaparentCovariate0Gradient[j][k];
                        }
                    }
                }
        }

    // contributions to gradient from numerator

    //  numerator

    // if allowing for linkage, loop through all combinations in all sibs
    int numerix = 0;
    int sibNumerix = 0;

    for (int posn = nsib + 1; posn < consistentHaps[nfamily].size(); posn += 3 + 8 * nsib) {

        int sumZero = 0;
        for (int sib = 0; sib < nsib; sib++)
            for (int i = 0; i < 8; i++) {
                int hap = consistentHaps[nfamily][posn+3+sib*8+i];
                if (!(options.genotype && i < 4)) {
                    sumZero += zero[hap];
                }
            }

        if (!sumZero) {
            vector<double> x(4, 0);
            double sumX = exp(numer[numerix++] - numerSum);

            vector<vector<double> > linear(nsib);
            vector<vector<double> > parentlinear(nsib);
            vector<vector<double> > alphaterm(nsib);
            for (int i = 0; i < nsib; i++) {
                linear[i].resize(4, 0);
                parentlinear[i].resize(4, 0);
                alphaterm[i].resize(4, 0);
            }

            int saveNumerix = sibNumerix;
            for (int sib = 0; sib < nsib; sib++) {

                if (sib == 0 || options.nolinkage)
                    for (int i = 0; i < 4; i += (twoVirtualControls ? 3 : 1)) {
                        x[i] = exp(sibNumer[sibNumerix++] - numerSum);
                    }

                int theHaps[2][2][2];
                int whichHaps[2][2][2];
                int ix = posn + 3 + sib * 8;
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++)
                        for (int k = 0; k < 2; k++) {
                            theHaps[i][j][k] = consistentHaps[nfamily][ix++];
                        }
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++)
                        for (int k = 0; k < 2; k++) {
                            if (permute[nfamily]) {
                                whichHaps[i][j][k] = theHaps[i][j][k];
                            } else {
                                whichHaps[i][j][k] = theHaps[i][!j][!k];
                            }
                        }

                realisationProb[nfamily].push_back(theHaps[0][0][0]);
                realisationProb[nfamily].push_back(theHaps[1][0][0]);
                realisationProb[nfamily].push_back(theHaps[0][1][1]);
                realisationProb[nfamily].push_back(theHaps[1][1][1]);

                bool MchrX = options.chrX && sibSex[sib] == MALE;

                // frequencies
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++) {
                        if (!options.genotype) {
                            if (!options.chrX) {
                                freqGradient[whichHaps[0][i][j]] -= x[i*2+j] / nsib * (family.sibship ? 2 : 1);
                            }
                            if (!MchrX && !family.sibship) {
                                freqGradient[whichHaps[0][!i][!j]] -= x[i*2+j] / nsib * (family.sibship ? 2 : 1);
                            }
                        }
                        freqGradient[whichHaps[1][i][j]] -= x[i*2+j] / nsib * (family.sibship ? 2 : 1);
                        if (!family.sibship) {
                            freqGradient[whichHaps[1][!i][!j]] -= x[i*2+j] / nsib * (family.sibship ? 2 : 1);
                        }
                    }

                // counts
                if (!options.genotype) {
                    if (!MchrX) {
                        familyCount[typeOfPhenotype!="quant" && sibTrait[sib]][whichHaps[0][0][0]] += sumX;
                    }
                    if (!options.chrX)
                        if (!family.sibship) {
                            familyCount[0][whichHaps[0][1][1]] += sumX;
                        }
                }
                familyCount[typeOfPhenotype!="quant" && sibTrait[sib]][whichHaps[1][0][0]] += sumX;
                if (!family.sibship) {
                    familyCount[0][whichHaps[1][1][1]] += sumX;
                }

                // haplotype effects
                for (int i = 0; i < 2; i++)
                    if (i == 1 || !options.genotype && !MchrX) {
                    	if (typeOfPhenotype == "polytomous") {
                    		if (sibTrait[sib] < K) betaGradient[whichHaps[i][0][0] + sibTrait[sib]*nhap] -= sumX;
                    		}
                        else betaGradient[whichHaps[i][0][0]] -= invVar * (sibTrait[sib] - sumIntercept) * sumX;
                        alphaGradient[whichHaps[i][0][0]] += invVar * sumIntercept * sumX;
                        if (invCov != 0)
                            for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                    betaGradient[whichHaps[i][0][0]] -= invCov * (sibTrait[sib2] - sumIntercept) * sumX;
                                    alphaGradient[whichHaps[i][0][0]] += invCov * sumIntercept * sumX;
                                }
                        double tmp = (invVar + (nsib - 1) * invCov) * (beta[whichHaps[i][0][0]] + alpha[whichHaps[i][0][0]]) * sumX;
                        alpha0Gradient += tmp;
                        for (int cov = 0; cov < betaCovariate.size(); cov++)
                            for (int level = 0; level < betaCovariate[cov].size(); level++) {
                                betaCovariate0Gradient[cov][level] += tmp;
                            }

                        // update linear models
                        if (typeOfPhenotype == "polytomous") {
                    		if (sibTrait[sib] < K) {
	                        for (int j = 0; j < 2; j++)
	                            for (int k = 0; k < 2; k++) {
	                                linear[sib][j*2+k] += beta[whichHaps[i][j][k] + sibTrait[sib]*nhap];
	                                parentlinear[sib][j*2+k] += betaparent[whichHaps[i][j][k] + sibTrait[sib]*nhap];
								}
							}
						}
						else {                    		
                        	for (int j = 0; j < 2; j++)
                            	for (int k = 0; k < 2; k++) {
                                	linear[sib][j*2+k] += beta[whichHaps[i][j][k]];
                                	alphaterm[sib][j*2+k] += alpha[whichHaps[i][j][k]];
                                	parentlinear[sib][j*2+k] += betaparent[whichHaps[i][j][k]];
                            	}
                        }
                    }

                // covariates
                for (int cov = 0; cov < covariateName.size(); cov++) {
                    vector<double> weight(2, 0);
                    int level = 0;
                    if (factor[cov] && sibCovariate[sib][cov] != 0) {
                        level = (int)sibCovariate[sib][cov] - 1;
                        weight[0] = 1.0;
                        weight[1] = 1.0;
                    } else {
                        weight[0] = sibCovariate[sib][cov];
                        weight[1] = sibCovariate[sib][cov];
                        if (covariateName[cov] == "parsex") {
                            weight[0] = (baseline[cov] != MALE);
                            weight[1] = (baseline[cov] != FEMALE);
                        }
                    }

                    if (confounder[cov]) {
                        for (int i = 0; i < 2; i++)
                            for (int j = 0; j < 2; j++) {
                                if (!options.genotype) {
                                    if (!options.chrX) {
                                        betaCovariateGradient[cov][level][whichHaps[0][i][j]] -= weight[0] * x[i*2+j] / nsib * (family.sibship ? 2 : 1);
                                        if (!family.sibship) {
                                            betaCovariateGradient[cov][level][whichHaps[0][!i][!j]] -= weight[0] * x[i*2+j] / nsib * (family.sibship ? 2 : 1);
                                        }
                                    } else {
                                        if (!MchrX && !family.sibship) {
                                            betaCovariateGradient[cov][level][whichHaps[0][0][j]] -= weight[0] * x[i*2+j] / nsib * (family.sibship ? 2 : 1);
                                        }
                                    }
                                }
                                betaCovariateGradient[cov][level][whichHaps[1][i][j]] -= weight[1] * x[i*2+j] / nsib * (family.sibship ? 2 : 1);
                                if (!family.sibship) {
                                    betaCovariateGradient[cov][level][whichHaps[1][!i][!j]] -= weight[1] * x[i*2+j] / nsib * (family.sibship ? 2 : 1);
                                }
                            }
                    } else {
                        for (int i = 0; i < 2; i++)
                            if (i == 1 || !options.genotype && !MchrX) {
                                betaCovariateGradient[cov][level][whichHaps[i][0][0]] -= invVar * (sibTrait[sib] - sumIntercept) * weight[i] * sumX;
                                if (invCov != 0)
                                    for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                            betaCovariateGradient[cov][level][whichHaps[i][0][0]] -= invCov * (sibTrait[sib2] - sumIntercept) * weight[i] * sumX;
                                        }

                                double tmp = (invVar + (nsib - 1) * invCov) * betaCovariate[cov][level][whichHaps[i][0][0]] * weight[i] * sumX;
                                alpha0Gradient += tmp;
                                for (int cov2 = 0; cov2 < betaCovariate.size(); cov2++)
                                    for (int level2 = 0; level2 < betaCovariate[cov2].size(); level2++) {
                                        betaCovariate0Gradient[cov2][level2] += tmp;
                                    }

                                for (int j = 0; j < 2; j++)
                                    for (int k = 0; k < 2; k++) {
                                        linear[sib][j*2+k] += weight[i] * betaCovariate[cov][level][whichHaps[i][j][k]];
                                        parentlinear[sib][j*2+k] += weight[i] * betaparentCovariate[cov][level][whichHaps[i][j][k]];
                                    }
                            }
                    }

                } // for cov

            } // for sib

            if (normal) {
                for (int sib = 0; sib < nsib; sib++) {
                    int theHaps[2][2][2];
                    int whichHaps[2][2][2];
                    int ix = posn + 3 + sib * 8;
                    for (int i = 0; i < 2; i++)
                        for (int j = 0; j < 2; j++)
                            for (int k = 0; k < 2; k++) {
                                theHaps[i][j][k] = consistentHaps[nfamily][ix++];
                            }
                    for (int i = 0; i < 2; i++)
                        for (int j = 0; j < 2; j++)
                            for (int k = 0; k < 2; k++) {
                                if (permute[nfamily]) {
                                    whichHaps[i][j][k] = theHaps[i][j][k];
                                } else {
                                    whichHaps[i][j][k] = theHaps[i][!j][!k];
                                }
                            }

                    bool MchrX = options.chrX && sibSex[sib] == MALE;
                    for (int i = 0; i < 2; i++)
                        if (i == 1 || !options.genotype && !MchrX) {
                            betaGradient[whichHaps[i][0][0]] += invVar * linear[sib][0] * sumX;
                            if (invCov != 0)
                                for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                        betaGradient[whichHaps[i][0][0]] += invCov * linear[sib2][0] * sumX;
                                    }
                        }
                    for (int cov = 0; cov < covariateName.size(); cov++) {
                        vector<double> weight(2, 0);
                        weight[0] = 0;
                        weight[1] = 0;
                        int level = 0;
                        if (factor[cov] && sibCovariate[sib][cov] != 0) {
                            level = (int)sibCovariate[sib][cov] - 1;
                            weight[0] = 1.0;
                            weight[1] = 1.0;
                        } else {
                            weight[0] = sibCovariate[sib][cov];
                            weight[1] = sibCovariate[sib][cov];
                            if (covariateName[cov] == "parsex") {
                                weight[0] = (baseline[cov] != MALE);
                                weight[1] = (baseline[cov] != FEMALE);
                            }
                        }
                        if (!confounder[cov]) {
                            for (int i = 0; i < 2; i++)
                                if (i == 1 || !options.genotype && !MchrX) {
                                    betaCovariateGradient[cov][level][whichHaps[i][0][0]] += invVar * linear[sib][0] * weight[i] * sumX;
                                    if (invCov != 0)
                                        for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                betaCovariateGradient[cov][level][whichHaps[i][0][0]] += invCov * linear[sib2][0] * weight[i] * sumX;
                                            }
                                }
                        }
                    } // for cov
                } // for sib
            } // if normal

            // form the virtual controls
            vector<vector<double> > virtualc(nsib);
            for (int sib = 0; sib < nsib; sib++) {
                virtualc[sib].resize(4, 0);
            }
            for (int sib = 0; sib < nsib; sib++)
                for (int v = 0; v < 4; v++) {
                	if (typeOfPhenotype == "polytomous") virtualc[sib][v] += linear[sib][v];
                    else virtualc[sib][v] += invVar * (sibTrait[sib] - sumIntercept) * linear[sib][v];
                    virtualc[sib][v] -= invVar * sumIntercept * alphaterm[sib][v];
                    if (normal) {
                        virtualc[sib][v] -= invVar * 0.5 * linear[sib][v] * linear[sib][v];
                    }
                    if (invCov != 0)
                        for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                virtualc[sib][v] += invCov * (sibTrait[sib2] - sumIntercept) * linear[sib][v];
                                virtualc[sib][v] -= invCov * sumIntercept * alphaterm[sib][v];
                                if (normal) {
                                    virtualc[sib][v] -= invCov * 0.5 * linear[sib2][v] * linear[sib][v];
                                }
                            } // for sib2
                } // for v, sib

            // if linkage, make all virtual controls equal
            if (!options.nolinkage)
                for (int v = 0; v < 4; v++) {
                    double sumVirtualc = 0;
                    for (int sib = 0; sib < nsib; sib++) {
                        sumVirtualc += virtualc[sib][v];
                    }
                    for (int sib = 0; sib < nsib; sib++) {
                        virtualc[sib][v] = sumVirtualc;
                    }
                }

            realisationProb[nfamily].push_back(sumX);

            // gradients from denominator of conditional likelihood
            if (!options.hhrr) {
                sibNumerix = saveNumerix;
                for (int sib = 0; sib < nsib; sib++) {

                    if (sib == 0 || options.nolinkage)
                        for (int i = 0; i < 4; i += (twoVirtualControls ? 3 : 1)) {
                            x[i] = exp(sibNumer[sibNumerix++] - numerSum);
                        }

                    int theHaps[2][2][2];
                    int whichHaps[2][2][2];
                    int ix = posn + 3 + sib * 8;
                    for (int i = 0; i < 2; i++)
                        for (int j = 0; j < 2; j++)
                            for (int k = 0; k < 2; k++) {
                                theHaps[i][j][k] = consistentHaps[nfamily][ix++];
                            }
                    for (int i = 0; i < 2; i++)
                        for (int j = 0; j < 2; j++)
                            for (int k = 0; k < 2; k++) {
                                if (permute[nfamily]) {
                                    whichHaps[i][j][k] = theHaps[i][j][k];
                                } else {
                                    whichHaps[i][j][k] = theHaps[i][!j][!k];
                                }
                            }

                    bool MchrX = options.chrX && sibSex[sib] == MALE;

                    double maxbeta2 = max(virtualc[sib][0], virtualc[sib][3]);
                    if (!twoVirtualControls) {
                        maxbeta2 = max(maxbeta2, virtualc[sib][1]);
                        maxbeta2 = max(maxbeta2, virtualc[sib][2]);
                    }
                    double condLhd[4];
                    double condDenom = 0;
                    if (twoVirtualControls) {
                        condDenom = exp(virtualc[sib][0] - maxbeta2) + exp(virtualc[sib][3] - maxbeta2);
                    } else {
                        condDenom = exp(virtualc[sib][0] - maxbeta2) + exp(virtualc[sib][1] - maxbeta2) + exp(virtualc[sib][2] - maxbeta2) + exp(virtualc[sib][3] - maxbeta2);
                    }
                    for (int i = 0; i < 4; i++)
                        if (!twoVirtualControls || i == 0 || i == 3) {
                            condLhd[i] = exp(virtualc[sib][i] - maxbeta2) / condDenom;
                        } else {
                            condLhd[i] = 0;
                        }

                    for (int i = 0; i < 2; i++)
                        for (int j = 0; j < 2; j++)
                            for (int k = 0; k < 2; k++) if (!twoVirtualControls || j == k) {
                                    if (i == 1 || !options.genotype && !MchrX) {
                                        // haplotype effects
                                        if (typeOfPhenotype == "polytomous") {
                                        	if (sibTrait[sib] < K) betaGradient[whichHaps[i][j][k] + sibTrait[sib]*nhap] += sumX * condLhd[j*2+k];
                                        	}
                                        else betaGradient[whichHaps[i][j][k]] += invVar * sumX * (sibTrait[sib] - sumIntercept) * condLhd[j*2+k];
                                        alphaGradient[whichHaps[i][j][k]] -= invVar * sumX * sumIntercept * condLhd[j*2+k];

                                        if (normal) {
                                            betaGradient[whichHaps[i][j][k]] -= invVar * sumX * linear[sib][j*2+k] * condLhd[j*2+k];
                                        }
                                        if (invCov != 0)
                                            for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                    betaGradient[whichHaps[i][j][k]] += invCov * sumX * (sibTrait[sib2] - sumIntercept) * condLhd[j*2+k];
                                                    alphaGradient[whichHaps[i][j][k]] -= invCov * sumX * sumIntercept * condLhd[j*2+k];
                                                    if (normal) {
                                                        betaGradient[whichHaps[i][j][k]] -= invCov * sumX * linear[sib2][j*2+k] * condLhd[j*2+k];
                                                    }
                                                }

                                        alpha0Gradient -= (invVar + (nsib - 1) * invCov) * sumX * (beta[whichHaps[i][j][k]] + alpha[whichHaps[i][j][k]]) * condLhd[j*2+k];
                                        for (int cov = 0; cov < betaCovariate.size(); cov++)
                                            for (int level = 0; level < betaCovariate[cov].size(); level++) {
                                                betaCovariate0Gradient[cov][level] -= (invVar + (nsib - 1) * invCov) * sumX * (beta[whichHaps[i][j][k]] + alpha[whichHaps[i][j][k]]) * condLhd[j*2+k];
                                            }

                                        if (typeOfPhenotype == "polytomous") {
                                        	if (sibTrait[sib] < K) betaparentGradient[whichHaps[i][j][k] + sibTrait[sib]*nhap] -= x[j*2+k];
                                        	}
                                        else betaparentGradient[whichHaps[i][j][k]] -= invVar * x[j*2+k] * (sibTrait[sib] - parentSumIntercept);
                                        if (options.hhrr) {
                                            alphaGradient[whichHaps[i][j][k]] += invVar * x[j*2+k] * parentSumIntercept;
                                        }
                                        if (normal) {
                                            betaparentGradient[whichHaps[i][j][k]] += invVar * x[j*2+k] * parentlinear[sib][j*2+k];
                                        }
                                        if (invCov != 0)
                                            for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                    betaparentGradient[whichHaps[i][j][k]] -= invCov * x[j*2+k] * (sibTrait[sib2] - parentSumIntercept);
                                                    if (options.hhrr) {
                                                        alphaGradient[whichHaps[i][j][k]] += invCov * x[j*2+k] * parentSumIntercept;
                                                    }
                                                    if (normal) {
                                                        betaparentGradient[whichHaps[i][j][k]] += invCov * x[j*2+k] * parentlinear[sib2][j*2+k];
                                                    }
                                                }

                                        betaparent0Gradient[0] += (invVar + (nsib - 1) * invCov) * betaparent[whichHaps[i][j][k]] * x[j*2+k];
                                        if (options.hhrr) {
                                            betaparent0Gradient[0] += (invVar + (nsib - 1) * invCov) * alpha[whichHaps[i][j][k]] * x[j*2+k];
                                        }
                                        for (int cov = 0; cov < betaCovariate.size(); cov++)
                                            for (int level = 0; level < betaCovariate[cov].size(); level++) {
                                                betaparentCovariate0Gradient[cov][level] += (invVar + (nsib - 1) * invCov) * betaparent[whichHaps[i][j][k]] * x[j*2+k];
                                                if (options.hhrr) {
                                                    betaparentCovariate0Gradient[cov][level] += (invVar + (nsib - 1) * invCov) * alpha[whichHaps[i][j][k]] * x[j*2+k];
                                                }
                                            }

                                        // covariate effects
                                        for (int cov = 0; cov < covariateName.size(); cov++) {
                                            vector<double> weight(2, 0);
                                            weight[0] = 0;
                                            weight[1] = 0;
                                            int level = 0;
                                            if (factor[cov] && sibCovariate[sib][cov] != 0) {
                                                level = (int)sibCovariate[sib][cov] - 1;
                                                weight[0] = 1.0;
                                                weight[1] = 1.0;
                                            } else {
                                                weight[0] = sibCovariate[sib][cov];
                                                weight[1] = sibCovariate[sib][cov];
                                                if (covariateName[cov] == "parsex") {
                                                    weight[0] = (baseline[cov] != MALE);
                                                    weight[1] = (baseline[cov] != FEMALE);
                                                }
                                            }
                                            if (!confounder[cov]) {
                                                betaCovariateGradient[cov][level][whichHaps[i][j][k]] += invVar * sumX * (sibTrait[sib] - sumIntercept) * weight[i] * condLhd[j*2+k];
                                                betaparentCovariateGradient[cov][level][whichHaps[i][j][k]] -= invVar * x[j*2+k] * (sibTrait[sib] - parentSumIntercept) * weight[i];
                                                if (normal) {
                                                    betaCovariateGradient[cov][level][whichHaps[i][j][k]] -= invVar * sumX * linear[sib][j*2+k] * condLhd[j*2+k] * weight[i];
                                                    betaparentCovariateGradient[cov][level][whichHaps[i][j][k]] += invVar * x[j*2+k] * parentlinear[sib][j*2+k] * weight[i];
                                                }

                                                double tmp = (invVar + (nsib - 1) * invCov) * sumX * betaCovariate[cov][level][whichHaps[i][j][k]] * weight[i] * condLhd[j*2+k];
                                                alpha0Gradient -= tmp;
                                                for (int cov2 = 0; cov2 < betaCovariate.size(); cov2++)
                                                    for (int level2 = 0; level2 < betaCovariate[cov2].size(); level2++) {
                                                        betaCovariate0Gradient[cov2][level2] -= tmp;
                                                    }
                                                tmp = (invVar + (nsib - 1) * invCov) * betaparentCovariate[cov][level][whichHaps[i][j][k]] * weight[i] * x[j*2+k];
                                                betaparent0Gradient[0] += tmp;
                                                for (int cov2 = 0; cov2 < betaCovariate.size(); cov2++)
                                                    for (int level2 = 0; level2 < betaCovariate[cov2].size(); level2++) {
                                                        betaparentCovariate0Gradient[cov2][level2] += tmp;
                                                    }

                                                if (invCov != 0)
                                                    for (int sib2 = 0; sib2 < nsib; sib2++) if (sib2 != sib) {
                                                            betaCovariateGradient[cov][level][whichHaps[i][j][k]] += invCov * sumX * (sibTrait[sib2] - sumIntercept) * weight[i] * condLhd[j*2+k];
                                                            betaparentCovariateGradient[cov][level][whichHaps[i][j][k]] -= invCov * x[j*2+k] * (sibTrait[sib2] - parentSumIntercept) * weight[i];
                                                            if (normal) {
                                                                betaCovariateGradient[cov][level][whichHaps[i][j][k]] -= invCov * sumX * linear[sib2][j*2+k] * condLhd[j*2+k] * weight[i];
                                                                betaparentCovariateGradient[cov][level][whichHaps[i][j][k]] += invCov * x[j*2+k] * parentlinear[sib2][j*2+k] * weight[i];
                                                            }
                                                        }  // for sib2
                                            } // if modifier
                                        } // for cov
                                    }
                                } // for k

                } // for sib
            } // if !hhrr
        } // if !sumZero
    } // for posn



    int ix = 0;
    for (int i = 0; i < betasize; i++) {
        gradient[ix++] += betaGradient[i];
    }
    for (int i = 0; i < nhap; i++) {
        gradient[ix++] += freqGradient[i];
    }
    if (haveBetaParent(options)) {
        for (int i = 0; i < betasize; i++) {
            gradient[ix++] += betaparentGradient[i];
        }
        if (haveBetaParent0(options)) {
            gradient[ix++] += betaparent0Gradient[0];
        	if (typeOfPhenotype == "polytomous")
        		for (int i = 1; i < K-1; i++)
            		gradient[ix++] += betaparent0Gradient[i];
        }
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
            if (!confounder[i] && haveFamilies && !options.hhrr)
                for (int k = 0; k < nhap; k++) {
                    gradient[ix++] += betaparentCovariateGradient[i][j][k];
                }
            if (typeOfPhenotype == "quant") {
                gradient[ix++] += betaCovariate0Gradient[i][j];
                if (haveFamilies && !options.hhrr) {
                    gradient[ix++] += betaparentCovariate0Gradient[i][j];
                }
            }
        }

}


