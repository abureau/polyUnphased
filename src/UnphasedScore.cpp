/* UnphasedScore.cpp - Association analysis in trios and unrelated

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

//  score

// score the data without returning the gradient, but with standard errors
void UnphasedAnalysis::score(UnphasedOptions &options, double &loglikelihood,
                             bool null) {

    vector<double> gradient(sizeOfGradient(options), 0);
    double llhd = 0;

    score(options, loglikelihood, gradient);
    if (null) {
        return;
    }

    int nhap = genoCode.size();
    int refIndex = reference.index(genoCode);
    int betasize;
    if (typeOfPhenotype == "polytomous") betasize = nhap*(K-1);
    else betasize = nhap;

    // calculate standard errors
    // numerical variance matrix
    vector<vector<double> > variance(gradient.size());
    for (int i = 0; i < variance.size(); i++) {
        variance[i].resize(variance.size(), 0);
    }
    numericalHessian(options, variance, gradient, llhd);

    // reduce to the free parameters
    vector<vector<double> > v(freeParameters(options));
    for (int i = 0; i < v.size(); i++) {
        v[i].resize(v.size(), 0);
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
    for (int i = 0; i < group.size(); i++)
        for (int j = 0; j < nhap; j++) if (!zero[j] && !(i && rare[j])) {
                vector<double> vv(freeParameters(options), 0);
                gradientFreeParameters(variance[j], vv, options);
                if (options.model != "commonmain" ||
                        i == 0 && options.condition.size() > 0) {
                    if (group[i][j] < group[i][refIndex])
                        for (int k = 0; k < v.size(); k++) {
                            v[group[i][j] + ndim[i]][k] += vv[k];
                        }
                    if (group[i][j] > group[i][refIndex])
                        for (int k = 0; k < v.size(); k++) {
                            v[group[i][j] + ndim[i] - 1][k] += vv[k];
                        }
                } else {
                    if (group[i][j] != group[i][refIndex]) {
                        for (int k = 0; k < v.size(); k++) {
                            v[ndim[options.condition.size()>0]][k] += vv[k];
                        }
                    }
                }
            }
    int ix = v.size() - 1;
    int iy = betasize;
    // frequency
    for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
            gradientFreeParameters(variance[i+iy], v[ix--], options);
        }
    iy += nhap;
    // betaparent
    if (haveBetaParent(options)) {
        for (int i = 0; i < nhap; i++)
            if (!zero[i] && i != refIndex) {
                if (typeOfPhenotype == "polytomous") {
                    for (int k = 0; k < K-1; k++)
                        gradientFreeParameters(variance[i+iy + k*nhap], v[ix - (K-2-k)*nhap], options);
                    }
                else gradientFreeParameters(variance[i+iy], v[ix], options);
                ix--;
            }
	    if (typeOfPhenotype == "polytomous") {
	    	ix -= nhap*(K-2);
	    	}
    	iy += betasize;
        if (haveBetaParent0(options)) {
            gradientFreeParameters(variance[iy++], v[ix--], options);
        }
    }
    // alpha
    if (haveAlpha(options)) {
        for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                gradientFreeParameters(variance[i+iy], v[ix--], options);
            }
        iy += nhap;
        if (haveAlpha0(options)) {
            gradientFreeParameters(variance[iy++], v[ix--], options);
        }
    }
    // covariates
    for (int j = 0; j < betaCovariate.size(); j++)
        for (int k = 0; k < betaCovariate[j].size(); k++) {
            for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                    gradientFreeParameters(variance[i+iy], v[ix--], options);
                }
            iy += nhap;
            if (!confounder[j] && haveFamilies && !options.hhrr) {
                for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                        gradientFreeParameters(variance[i+iy], v[ix--], options);
                    }
                iy += nhap;
            }
            if (typeOfPhenotype == "quant") {
                gradientFreeParameters(variance[iy++], v[ix--], options);

                if (haveFamilies && !options.hhrr) {
                    gradientFreeParameters(variance[iy++], v[ix--], options);
                }
            }
        }

    // invert variance matrix
    svdinvert(v);

    // standard errors of beta parameters
    stderror = 0;
    for (int i = 0; i < group.size(); i++)
        for (int j = 0; j < nhap; j++) if (!zero[j] && !(i && rare[j])) {
                if (options.model != "commonmain" ||
                        i == 0 && options.condition.size() > 0) {
                    if (group[i][j] < group[i][refIndex]) {
                        stderror[j] += v[group[i][j] + ndim[i]][group[i][j] + ndim[i]];
			if (typeOfPhenotype == "polytomous") 
				for (int k = 1; k < K-1; k++)
					stderror[j + k*nhap] += v[group[i][j] + ndim[i] + k*thisndim][group[i][j] + ndim[i] + k*thisndim];
                    }
                    if (group[i][j] > group[i][refIndex]) {
                        stderror[j] += v[group[i][j] + ndim[i] - 1][group[i][j] + ndim[i] - 1];
			if (typeOfPhenotype == "polytomous") 
				for (int k = 1; k < K-1; k++)
					stderror[j + k*nhap] += v[group[i][j] -1 + ndim[i] + k*thisndim][group[i][j] + ndim[i] -1 + k*thisndim];
                    }
                } else if (group[i][j] != group[i][refIndex]) {
                    stderror[j] += v[ndim[options.condition.size()>0]][ndim[options.condition.size()>0]];
			if (typeOfPhenotype == "polytomous") 
				for (int k = 1; k < K-1; k++)
					stderror[j + k*nhap] += v[ndim[options.condition.size()>0] + k*thisndim][ndim[options.condition.size()>0] + k*thisndim];
                }
            }
    stderror = sqrt(abs(stderror));

    // standard errors of covariate effects
    // used for testing whether an effect is zero
    for (int i = 0; i < betaCovariate.size(); i++)
        for (int j = 0; j < betaCovariate[i].size(); j++) {
            stderrorCovariate[i][j] = 0;
        }

    ix = v.size() - 1;
    // frequency
    for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
            ix--;
        }
    // betaparent
    if (haveBetaParent(options)) {
        for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                ix-=(K-1);
            }
        if (haveBetaParent0(options)) {
            ix--;
        }
    }
    // alpha
    if (haveAlpha(options)) {
        for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                ix--;
            }
        if (haveAlpha0(options)) {
            ix--;
        }
    }
    // covariates
    for (int j = 0; j < betaCovariate.size(); j++)
        for (int k = 0; k < betaCovariate[j].size(); k++) {
            for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                    stderrorCovariate[j][k][i] = sqrt(abs(v[ix][ix]));
                    ix--;
                }
            // betaparentCovariate
            if (haveFamilies && !confounder[j] && !options.hhrr)
                for (int i = 0; i < nhap; i++) if (!zero[i] && i != refIndex) {
                        ix--;
                    }
            // baseline
            if (typeOfPhenotype == "quant") {
                ix--;
                if (haveFamilies && !options.hhrr) {
                    ix--;
                }
            }
        }
}


// score the data and return loglikelihood and gradient
void UnphasedAnalysis::score(UnphasedOptions &options, double &loglikelihood,
                             vector<double> &gradient) {
    int phenoval;
    
    // empirical variance
    empiricalVariance.resize(gradient.size());
    for (int i = 0; i < empiricalVariance.size(); i++) {
        empiricalVariance[i].resize(gradient.size());
        for (int j = 0; j < empiricalVariance.size(); j++) {
            empiricalVariance[i][j] = 0;
        }
    }

    if (options.llhd) {
        cout << "frequency ";
        for (int i = 0; i < genoCode.size(); i++) {
            cout << frequency[i] << " ";
        }
        cout << "beta ";
        for (int i = 0; i < betasize; i++) {
            cout << beta[i] << " ";
        }
        cout << endl;
        cout << "betaparent ";
        for (int i = 0; i < betasize; i++) {
            cout << betaparent[i] << " ";
        }
        cout << endl;
        cout << "alpha ";
        for (int i = 0; i < genoCode.size(); i++) {
            cout << alpha[i] << " ";
        }
        cout << endl;
    }

    int nhap = genoCode.size();

    loglikelihood = 0;

    // loop through the families
    // global denominators for trios and ASPs
    bool haveGlobal[2];
    haveGlobal[0] = false;
    haveGlobal[1] = false;
    double globalDenomFamily[2];
    vector<double> globalGradFamily[2];
    familyCount[0] = 0;
    familyCount[1] = 0;
    for (int nfamily = 0; nfamily < familylist.size(); nfamily++) {
        NuclearFamily *family = &familylist[nfamily];
        if (family->status == "") {

            // total gradient for the family
            //vector<double> familyGradient(nparam*nhap,0);
            vector<double> familyGradient(sizeOfGradient(options), 0);

            // set up trait values and sexes for the sibs
            vector<double> sibTrait;
            vector<int> sibSex;
            vector<vector<double> > sibCovariate;
            for (int i = 0; consistentHaps[nfamily][i] != -1; i++) {
                int sib = consistentHaps[nfamily][i];
                if (family->sibs[sib].usable) {
                    // trait value
                    if (typeOfPhenotype == "quant") {
                        sibTrait.push_back(family->sibs[sib].trait[currentphenotype]);
                    } else if (typeOfPhenotype == "polytomous") {
                    	if (family->sibs[sib].affection[currentphenotype] == AFFECTED) {
                    	 	if (family->sibs[sib].affection[currentphenotype2] == AFFECTED) phenoval = 3;
                    	 	else phenoval = 1;
                    	 	} else { if (family->sibs[sib].affection[currentphenotype2] == AFFECTED) phenoval = 2;
                    	 				else phenoval = 4;
                    	 				}
                    	sibTrait.push_back(phenoval); 		
                    	} else {                    	
                        sibTrait.push_back(family->sibs[sib].affection[currentphenotype] == AFFECTED);
                    }
                    // sex
                    sibSex.push_back(family->sibs[sib].sex);

                    // covariate values
                    vector<double> subjectCovariate;
                    for (int cov = 0; cov < covariateName.size(); cov++) {
                        string name = covariateName[cov];
                        if (factor[cov]) {
                            int level = 0;
                            double value = family->sibs[sib].trait[traithash[name]];
                            // the baseline value will be given level 0
                            // other values in the order in which they occur in factorLevels
                            for (set<double>::iterator ilevel = factorLevels[cov].begin(); ilevel != factorLevels[cov].find(value); ilevel++) if (*ilevel != baseline[cov]) {
                                    level++;
                                }
                            if (value == baseline[cov]) {
                                level = 0;
                            } else {
                                level++;
                            }
                            subjectCovariate.push_back(level);
                        } else {
                            if (name == "sibsex") {
                                subjectCovariate.push_back(family->sibs[sib].sex != baseline[cov]);
                            } else if (name == "parsex") {
                                subjectCovariate.push_back(0);
                            } else {
                                subjectCovariate.push_back(family->sibs[sib].trait[traithash[name]]);
                            }
                        }
                    }
                    sibCovariate.push_back(subjectCovariate);
                }
            }

            // score the family
            double prob = 0;
            scoreFamily(*family, nfamily, sibTrait, sibSex, sibCovariate,
                        options, prob, familyGradient,
                        haveGlobal, globalDenomFamily, globalGradFamily);
            if (options.llhd) {
                cout << "family " << family->name << " prob " << prob << " ";
            }
            loglikelihood += prob;
            if (options.llhd) {
                cout << loglikelihood << endl;
            }
            if (options.llhd) {
                for (int i = 0; i < gradient.size(); i++) {
                    cout << familyGradient[i] << " ";
                }
                cout << endl;
            }
            for (int i = 0; i < gradient.size(); i++) {
                gradient[i] += familyGradient[i];
            }

            // empirical variance
            for (int i = 0; i < gradient.size(); i++)
                for (int j = 0; j < gradient.size(); j++) {
                    empiricalVariance[i][j] += familyGradient[i] * familyGradient[j];
                }
        }
    }

    // loop through the unrelateds
    int nsubject = 0;
    haveGlobal[0] = false;
    haveGlobal[1] = false;
    double globalDenomUnrelated[2];
    vector<double> globalGradUnrelated[2];
    unrelatedCount[0] = 0;
    unrelatedCount[1] = 0;
    for (int nsubject = 0; nsubject < unrelateds.size(); nsubject++) {
        Subject *subject = &unrelateds[nsubject];
        if (subject->usable) {

            // total gradient for the subject
            //vector<double> subjectGradient(nparam*nhap,0);
            vector<double> subjectGradient(sizeOfGradient(options), 0);

            // trait value
            double subjectTrait;
            if (typeOfPhenotype == "quant") {
                subjectTrait = unrelateds[permuteUnrelated[nsubject]].trait[currentphenotype];
            } else if (typeOfPhenotype == "polytomous") {
            	if (unrelateds[permuteUnrelated[nsubject]].affection[currentphenotype] == AFFECTED) {
            		if (unrelateds[permuteUnrelated[nsubject]].affection[currentphenotype2] == AFFECTED) subjectTrait = 3;
            		else subjectTrait = 1;
            		} else { if (unrelateds[permuteUnrelated[nsubject]].affection[currentphenotype2] == AFFECTED) subjectTrait = 2;
            					else subjectTrait = 4;
            					}
            	}
            else {
                subjectTrait = (unrelateds[permuteUnrelated[nsubject]].affection[currentphenotype] == AFFECTED);
            }

            // covariate values
            vector<double> subjectCovariate(0);
            for (int cov = 0; cov < covariateName.size(); cov++) {
                string name = covariateName[cov];
                //if (subject->trait.size()<=traithash[name]) cout << "ALARM " << subject->trait.size() << " " << traithash[name] << endl;
                if (factor[cov]) {
                    int level = 0;
                    double value = subject->trait[traithash[name]];
                    // the baseline value will be given level 0
                    // other values in the order in which they occur in factorLevels
                    for (set<double>::iterator ii = factorLevels[cov].begin(); ii != factorLevels[cov].find(subject->trait[traithash[name]]); ii++) if (*ii != baseline[cov]) {
                            level++;
                        }
                    if (value == baseline[cov]) {
                        level = 0;
                    } else {
                        level++;
                    }
                    subjectCovariate.push_back(level);
                } else {
                    if (name == "sibsex") {
                        subjectCovariate.push_back(subject->sex != baseline[cov]);
                    } else {
                        if (name == "parsex") {
                            subjectCovariate.push_back(0);
                        } else {
                            subjectCovariate.push_back(subject->trait[traithash[name]]);
                        }
                    }
                }
            }

            // score the subject
            double prob = 0;
            if (consistentHaps[nsubject+familylist.size()].size() > 0) {
                scoreUnrelated(*subject, nsubject + familylist.size(), subjectTrait, subjectCovariate, options, prob, subjectGradient, haveGlobal, globalDenomUnrelated, globalGradUnrelated);
                if (options.llhd) {
                    cout << subject->id << " " << prob << " ";
                }
                if (prob) {
                    loglikelihood += prob;
                    if (options.llhd) {
                        cout << loglikelihood << endl;
                    }
                    for (int i = 0; i < gradient.size(); i++) {
                        gradient[i] += subjectGradient[i];
                    }
                }
            }
        }
    }

}



//  numericalHessian

void UnphasedAnalysis::numericalHessian(UnphasedOptions &options, vector<vector<double> > &v, vector<double> &g, double llhd) {
    const double epsilon = 1e-8;
    int nhap = genoCode.size();

    // working score vectors
    vector<double> gradient(g.size(), 0);

    int ix = 0;
    for (int i = 0; i < nhap; ix++, i++) if (!zero[i]) {
            beta[i] += epsilon;
            if (options.hhrr) {
                betaparent[i] = beta[i];
            }
            for (int j = 0; j < gradient.size(); j++) {
                gradient[j] = 0;
            }
            score(options, llhd, gradient);
            for (int j = 0; j < gradient.size(); j++) {
                v[ix][j] = (gradient[j] - g[j]) / epsilon;
            }
            beta[i] -= epsilon;
            if (options.hhrr) {
                betaparent[i] = beta[i];
            }
        }

    for (int i = 0; i < nhap; ix++, i++) if (!zero[i]) {
            frequency[i] += epsilon;
            for (int j = 0; j < gradient.size(); j++) {
                gradient[j] = 0;
            }
            score(options, llhd, gradient);
            for (int j = 0; j < gradient.size(); j++) {
                v[ix][j] = (gradient[j] - g[j]) / epsilon;
            }
            frequency[i] -= epsilon;
        }

    if (haveBetaParent(options)) {
        for (int i = 0; i < nhap; ix++, i++) if (!zero[i]) {
                betaparent[i] += epsilon;
                for (int j = 0; j < gradient.size(); j++) {
                    gradient[j] = 0;
                }
                score(options, llhd, gradient);
                for (int j = 0; j < gradient.size(); j++) {
                    v[ix][j] = (gradient[j] - g[j]) / epsilon;
                }
                betaparent[i] -= epsilon;
            }
        if (haveBetaParent0(options)) {
            betaparent0 += epsilon;
            for (int j = 0; j < gradient.size(); j++) {
                gradient[j] = 0;
            }
            score(options, llhd, gradient);
            for (int j = 0; j < gradient.size(); j++) {
                v[ix][j] = (gradient[j] - g[j]) / epsilon;
            }
            betaparent0 -= epsilon;
            ix++;
        }
    }

    if (haveAlpha(options)) {
        for (int i = 0; i < nhap; ix++, i++) if (!zero[i]) {
                alpha[i] += epsilon;
                for (int j = 0; j < gradient.size(); j++) {
                    gradient[j] = 0;
                }
                score(options, llhd, gradient);
                for (int j = 0; j < gradient.size(); j++) {
                    v[ix][j] = (gradient[j] - g[j]) / epsilon;
                }
                alpha[i] -= epsilon;
            }
        if (haveAlpha0(options)) {
            alpha0 += epsilon;
            for (int j = 0; j < gradient.size(); j++) {
                gradient[j] = 0;
            }
            score(options, llhd, gradient);
            for (int j = 0; j < gradient.size(); j++) {
                v[ix][j] = (gradient[j] - g[j]) / epsilon;
            }
            alpha0 -= epsilon;
            ix++;
        }
    }

    for (int i = 0; i < betaCovariate.size(); i++)
        for (int j = 0; j < betaCovariate[i].size(); j++) {
            for (int k = 0; k < nhap; ix++, k++) if (!zero[k]) {
                    betaCovariate[i][j][k] += epsilon;
                    if (options.hhrr) {
                        betaparentCovariate[i][j][k] = betaCovariate[i][j][k];
                    }
                    for (int l = 0; l < gradient.size(); l++) {
                        gradient[l] = 0;
                    }
                    score(options, llhd, gradient);
                    for (int l = 0; l < gradient.size(); l++) {
                        v[ix][l] = (gradient[l] - g[l]) / epsilon;
                    }
                    betaCovariate[i][j][k] -= epsilon;
                    if (options.hhrr) {
                        betaparentCovariate[i][j][k] = betaCovariate[i][j][k];
                    }
                }
            if (!confounder[i] && haveFamilies && !options.hhrr) {
                for (int k = 0; k < nhap; ix++, k++) if (!zero[k]) {
                        betaparentCovariate[i][j][k] += epsilon;
                        for (int l = 0; l < gradient.size(); l++) {
                            gradient[l] = 0;
                        }
                        score(options, llhd, gradient);
                        for (int l = 0; l < gradient.size(); l++) {
                            v[ix][l] = (gradient[l] - g[l]) / epsilon;
                        }
                        betaparentCovariate[i][j][k] -= epsilon;
                    }
            }
            if (typeOfPhenotype == "quant") {
                betaCovariate0[i][j] += epsilon;
                if (options.hhrr) {
                    betaparentCovariate0[i][j] = betaCovariate0[i][j];
                }
                for (int l = 0; l < gradient.size(); l++) {
                    gradient[l] = 0;
                }
                score(options, llhd, gradient);
                for (int l = 0; l < gradient.size(); l++) {
                    v[ix][l] = (gradient[l] - g[l]) / epsilon;
                }
                betaCovariate0[i][j] -= epsilon;
                if (options.hhrr) {
                    betaparentCovariate0[i][j] = betaCovariate0[i][j];
                }
                ix++;
                if (haveFamilies && !options.hhrr) {
                    betaparentCovariate0[i][j] += epsilon;
                    for (int l = 0; l < gradient.size(); l++) {
                        gradient[l] = 0;
                    }
                    score(options, llhd, gradient);
                    for (int l = 0; l < gradient.size(); l++) {
                        v[ix][l] = (gradient[l] - g[l]) / epsilon;
                    }
                    betaparentCovariate0[i][j] -= epsilon;
                    ix++;
                }
            }
        }
}


