/* UnphasedAnalysisModels.cpp - Association analysis in trios and unrelateds

   Extension of the original code enabling analysis of polytomous phenotypes
   
   Copyright (c) 2015 Alexandre Bureau
   Université Laval
   Département de médecine sociale et préventive
   1050 rue de la Médecine
   Québec (Québec), G1V 0A6, Canada
   alexandre.bureau@fmed.ulaval.ca
   
   Original code:
   
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

// local includes
#include "stl.h"
#include "UnphasedOptions.h"
#include "UnphasedAnalysis.h"
#include "pvalues.h"
#include "matinv.h"

//  weq Haplotype

static bool weq(const Haplotype &h1, const Haplotype &h2) {
    if (h1.size() != h2.size()) {
        return(false);
    }
    for (int i = 0; i < h1.size(); i++)
        if (h1[i] != h2[i] && h1[i] != 0 && h2[i] != 0) {
            return(false);
        }
    return(true);
}



//  fullNull

double UnphasedAnalysis::fullNull(UnphasedOptions &options,
                                  const string &which, int &df) {
    Haplotype condspecific;
    condspecific = options.condspecific;
    int conditionsize = options.condition.size() * (options.condgenotype ? 2 : 1);
    Haplotype specific;
    specific = options.specific;
    int testsize = options.window * (options.genotype ? 2 : 1);

    for (int i = 0; i < group.size(); i++)
        for (int j = 0; j < group[i].size(); j++) {
            group[i][j] = 0;
        }
    for (int i = 0; i < beta.size(); i++) {
        beta[i] = 0;
    }

    vector<Haplotype> partialCode;
    for (int i = 0; i < genoCode.size(); i++) if (!zero[i]) {
            Haplotype h = genoCode[i];
            Haplotype condhap;
            condhap.assign(h.begin(), h.begin() + conditionsize);
            // if rare, or not the conditioning haplotype, put in own group
            if (rare[i] ||
                    (condspecific.size() && !weq(condhap, condspecific))) {
                condhap.assign(h.begin(), h.begin() + conditionsize + testsize);
                condhap.enter(partialCode);
                group[0][i] = condhap.index(partialCode);
            }
            // otherwise put in conditioning group
            else {
                condhap.enter(partialCode);
                group[0][i] = condhap.index(partialCode);
            }
        }

    df = partialCode.size() - 1;
	if (typeOfPhenotype == "polytomous")  df *= K-1;

    return(getloglikelihood("Null", options, which, true));
}



//  fullAlternative

double UnphasedAnalysis::fullAlternative(UnphasedOptions &options,
        const string &which, int &df) {

    Haplotype condspecific;
    condspecific = options.condspecific;
    int conditionsize = options.condition.size() * (options.condgenotype ? 2 : 1);
    Haplotype specific;
    specific = options.specific;
    int testsize = options.window * (options.genotype ? 2 : 1);

    for (int i = 0; i < group.size(); i++)
        for (int j = 0; j < group[i].size(); j++) {
            group[i][j] = 0;
        }
    for (int i = 0; i < beta.size(); i++) {
        beta[i] = 0;
    }

    vector<Haplotype> partialCode;
    for (int i = 0; i < genoCode.size(); i++) if (!zero[i]) {
            Haplotype h = genoCode[i];
            Haplotype condhap;
            condhap.assign(h.begin(), h.begin() + conditionsize);
            bool condmatch = !rare[i] &&
                             (!condspecific.size() || weq(condhap, condspecific));
            Haplotype testhap;
            testhap.assign(h.begin() + conditionsize, h.begin() + conditionsize + testsize);
            bool testmatch = !specific.size() || weq(testhap, specific);
            // if rare, or not a conditioning haplotype, or a test haplotype
            // put in own group
            if (rare[i] ||
                    (condspecific.size() && !weq(condhap, condspecific)) ||
                    (!specific.size() || weq(testhap, specific))) {
                testhap.assign(h.begin(), h.begin() + conditionsize + testsize);
                testhap.enter(partialCode);
                group[0][i] = testhap.index(partialCode);
            }
            // otherwise put in conditioning group
            // ie has conditioning haplotype but not test haplotype
            else {
                condhap.enter(partialCode);
                group[0][i] = condhap.index(partialCode);
            }
        }

    df = partialCode.size() - 1;
	if (typeOfPhenotype == "polytomous")  df *= K-1;

    return(getloglikelihood("Alternative", options, which, false));
}



//  haploMainAlternative

double UnphasedAnalysis::haploMainAlternative(UnphasedOptions &options,
        const string &which,
        int &df) {

    Haplotype condspecific;
    condspecific = options.condspecific;
    int conditionsize = options.condition.size() * (options.condgenotype ? 2 : 1);
    Haplotype specific;
    specific = options.specific;
    int testsize = options.window * (options.genotype ? 2 : 1);

    for (int i = 0; i < group.size(); i++)
        for (int j = 0; j < group[i].size(); j++) {
            group[i][j] = 0;
        }
    for (int i = 0; i < beta.size(); i++) {
        beta[i] = 0;
    }

    // main effects of the conditioning haplotype
    vector<Haplotype> partialCode;
    vector<bool> owngroup(genoCode.size(), false);
    for (int i = 0; i < genoCode.size(); i++) if (!zero[i]) {
            Haplotype h = genoCode[i];
            Haplotype condhap;
            condhap.assign(h.begin(), h.begin() + conditionsize);
            // if rare, or not the conditioning haplotype, put in own group
            if (rare[i] ||
                    (condspecific.size() && !weq(condhap, condspecific))) {
                owngroup[i] = true;
                condhap.assign(h.begin(), h.begin() + conditionsize + testsize);
                condhap.enter(partialCode);
                group[0][i] = condhap.index(partialCode);
            }
            // otherwise put in conditioning group
            else {
                condhap.enter(partialCode);
                group[0][i] = condhap.index(partialCode);
            }
        }

    df = partialCode.size() - 1;

    // main effects of the test haplotype
    partialCode.clear();
    for (int i = 0; i < genoCode.size(); i++) if (!zero[i]) {
            Haplotype h = genoCode[i];
            Haplotype condhap;
            condhap.assign(h.begin(), h.begin() + conditionsize);
            Haplotype testhap;
            testhap.assign(h.begin() + conditionsize, h.begin() + conditionsize + testsize);
            // if not already in own group
            if (!owngroup[i]) {
                // group haplotypes having the specific test haplotype
                if (!specific.size() || weq(testhap, specific)) {
                    testhap.enter(partialCode);
                    group[1][i] = testhap.index(partialCode);
                }
                // otherwise group everything together
                else {
                    testhap.clear();
                    testhap.enter(partialCode);
                    group[1][i] = testhap.index(partialCode);
                }
            }
        }
    if (partialCode.size()) {
        df += partialCode.size() - 1;
    }
	if (typeOfPhenotype == "polytomous")  df *= K-1;

    return(getloglikelihood("Alternative", options, which, false));
}



//  alleleMainAlternative

double UnphasedAnalysis::alleleMainAlternative(UnphasedOptions &options,
        const string &which,
        int &df) {

    Haplotype condspecific;
    condspecific = options.condspecific;
    int conditionsize = options.condition.size() * (options.condgenotype ? 2 : 1);
    Haplotype specific;
    specific = options.specific;
    int testsize = options.window * (options.genotype ? 2 : 1);

    for (int i = 0; i < group.size(); i++)
        for (int j = 0; j < group[i].size(); j++) {
            group[i][j] = 0;
        }
    for (int i = 0; i < beta.size(); i++) {
        beta[i] = 0;
    }

    // main effects of the conditioning haplotype
    vector<Haplotype> partialCode;
    vector<bool> owngroup(genoCode.size(), false);
    df = 0;
    if (options.condition.size()) {
        for (int i = 0; i < genoCode.size(); i++) if (!zero[i]) {
                Haplotype h = genoCode[i];
                Haplotype condhap;
                condhap.assign(h.begin(), h.begin() + conditionsize);
                // if rare, or not the conditioning haplotype, put in own group
                if (rare[i] ||
                        (condspecific.size() && !weq(condhap, condspecific))) {
                    owngroup[i] = true;
                    condhap.assign(h.begin(), h.begin() + conditionsize + testsize);
                    condhap.enter(partialCode);
                    group[0][i] = condhap.index(partialCode);
                }
                // otherwise put in conditioning group
                else {
                    condhap.enter(partialCode);
                    group[0][i] = condhap.index(partialCode);
                }
            }
        df = partialCode.size() - 1;
    }

    // main effects of the test alleles
    for (int k = 0; k < options.window; k++) {
        partialCode.clear();
        for (int i = 0; i < genoCode.size(); i++) if (!zero[i]) {
                Haplotype h = genoCode[i];
                Haplotype condhap;
                condhap.assign(h.begin(), h.begin() + conditionsize);
                Haplotype testhap;
                testhap.resize(options.genotype ? 2 : 1);
                testhap[0] = h[conditionsize+k*(options.genotype?2:1)];
                if (options.genotype) {
                    testhap[1] = h[conditionsize+k*(options.genotype?2:1)+1];
                }
                // if not already in own group
                if (!owngroup[i]) {
                    // group haplotypes having the same allele
                    // the specific option is ignored
                    testhap.enter(partialCode);
                    group[k+(options.condition.size()>0)][i] = testhap.index(partialCode);
                }
            }
        if (partialCode.size() && !(options.model == "commonmain" && k > 0)) {
            df += partialCode.size() - 1;
        }
    }
    if (typeOfPhenotype == "polytomous")  df *= K-1;

    return(getloglikelihood("Alternative", options, which, false));
}



//  gxgNull

double UnphasedAnalysis::gxgNull(UnphasedOptions &options, const string &which,
                                 int &df) {

    Haplotype condspecific;
    condspecific = options.condspecific;
    int conditionsize = options.condition.size() * (options.condgenotype ? 2 : 1);
    Haplotype specific;
    specific = options.specific;
    int testsize = options.window * (options.genotype ? 2 : 1);

    for (int i = 0; i < group.size(); i++)
        for (int j = 0; j < group[i].size(); j++) {
            group[i][j] = 0;
        }
    for (int i = 0; i < beta.size(); i++) {
        beta[i] = 0;
    }

    // main effects of the conditioning haplotype
    vector<Haplotype> partialCode;
    for (int i = 0; i < genoCode.size(); i++) if (!zero[i]) {
            Haplotype h = genoCode[i];
            Haplotype condhap;
            condhap.assign(h.begin(), h.begin() + conditionsize);
            // if rare, put in own group
            if (rare[i]) {
                condhap.assign(h.begin(), h.begin() + conditionsize + testsize);
                condhap.enter(partialCode);
                group[0][i] = condhap.index(partialCode);
            }
            // otherwise put in conditioning group
            else {
                condhap.enter(partialCode);
                group[0][i] = condhap.index(partialCode);
            }
        }

    df = partialCode.size() - 1;

    // main effects of the test haplotype
    partialCode.clear();
    for (int i = 0; i < genoCode.size(); i++) if (!zero[i]) {
            Haplotype h = genoCode[i];
            Haplotype testhap;
            testhap.assign(h.begin() + conditionsize, h.begin() + conditionsize + testsize);
            // if not already in own group
            if (!rare[i]) {
                // group haplotypes having the specific test haplotype
                testhap.enter(partialCode);
                group[1][i] = testhap.index(partialCode);
            }
        }
    if (partialCode.size()) {
        df += partialCode.size() - 1;
    }
	if (typeOfPhenotype == "polytomous")  df *= K-1;

    return(getloglikelihood("Null", options, which, true));
}



//  gxgAlternative

double UnphasedAnalysis::gxgAlternative(UnphasedOptions &options,
                                        const string &which, int &df) {

    Haplotype condspecific;
    condspecific = options.condspecific;
    int conditionsize = options.condition.size() * (options.condgenotype ? 2 : 1);
    Haplotype specific;
    specific = options.specific;
    int testsize = options.window * (options.genotype ? 2 : 1);

    for (int i = 0; i < group.size(); i++)
        for (int j = 0; j < group[i].size(); j++) {
            group[i][j] = 0;
        }
    for (int i = 0; i < beta.size(); i++) {
        beta[i] = 0;
    }

    // main effects of the conditioning haplotype
    vector<Haplotype> partialCode;
    vector<bool> owngroup(genoCode.size(), false);
    for (int i = 0; i < genoCode.size(); i++) if (!zero[i]) {
            Haplotype h = genoCode[i];
            Haplotype condhap;
            condhap.assign(h.begin(), h.begin() + conditionsize);
            Haplotype testhap;
            testhap.assign(h.begin() + conditionsize, h.begin() + conditionsize + testsize);
            // if rare, or a haplotype of interest, put in own group
            if (rare[i] ||
                    (!condspecific.size() || weq(condspecific, condhap)) &&
                    (!specific.size() || weq(specific, testhap))) {
                owngroup[i] = true;
                condhap.assign(h.begin(), h.begin() + conditionsize + testsize);
                condhap.enter(partialCode);
                group[0][i] = condhap.index(partialCode);
            }
            // otherwise put in conditioning group
            else {
                condhap.enter(partialCode);
                group[0][i] = condhap.index(partialCode);
            }
        }

    df = partialCode.size() - 1;

    // main effects of the test haplotype
    partialCode.clear();
    for (int i = 0; i < genoCode.size(); i++) if (!zero[i]) {
            Haplotype h = genoCode[i];
            Haplotype testhap;
            testhap.assign(h.begin() + conditionsize, h.begin() + conditionsize + testsize);
            // if not already in own group
            if (!owngroup[i]) {
                // group haplotypes having the specific test haplotype
                testhap.enter(partialCode);
                group[1][i] = testhap.index(partialCode);
            }
        }
    if (partialCode.size()) {
        df += partialCode.size() - 1;
    }
	if (typeOfPhenotype == "polytomous")  df *= K-1;

    return(getloglikelihood("Alternative", options, which, false));
}



//  pairwiseNull

double UnphasedAnalysis::pairwiseNull(UnphasedOptions &options, const string &which,
                                      int &df) {
    Haplotype compare1;
    Haplotype compare2;
    compare1 = options.compare1;
    compare2 = options.compare2;
    int conditionsize = options.condition.size() * (options.condgenotype ? 2 : 1);
    int testsize = options.window * (options.genotype ? 2 : 1);

    for (int i = 0; i < group.size(); i++)
        for (int j = 0; j < group[i].size(); j++) {
            group[i][j] = 0;
        }
    for (int i = 0; i < beta.size(); i++) {
        beta[i] = 0;
    }

    // group all haplotypes matching either target
    vector<Haplotype> partialCode;
    for (int i = 0; i < genoCode.size(); i++) if (!zero[i]) {
            Haplotype h = genoCode[i];
            Haplotype testhap;
            testhap.assign(h.begin(), h.begin() + conditionsize + testsize);
            // if rare, or not matching either< target, put in own group
            if (rare[i] ||
                    !weq(testhap, compare1) && !weq(testhap, compare2)) {
                testhap.enter(partialCode);
                group[0][i] = testhap.index(partialCode);
            }
            // otherwise put in one group
            else {
                testhap.clear();
                testhap.enter(partialCode);
                group[0][i] = testhap.index(partialCode);
            }
        }

    df = partialCode.size() - 1;
	if (typeOfPhenotype == "polytomous")  df *= K-1;
	
    return(getloglikelihood("Null", options, which, true));

}



//  pairwiseAlternative

double UnphasedAnalysis::pairwiseAlternative(UnphasedOptions &options, const string &which, int &df) {
    Haplotype compare1;
    Haplotype compare2;
    compare1 = options.compare1;
    compare2 = options.compare2;
    int conditionsize = options.condition.size() * (options.condgenotype ? 2 : 1);
    int testsize = options.window * (options.genotype ? 2 : 1);

    for (int i = 0; i < group.size(); i++)
        for (int j = 0; j < group[i].size(); j++) {
            group[i][j] = 0;
        }
    for (int i = 0; i < beta.size(); i++) {
        beta[i] = 0;
    }

    // make two groups for the two targets
    // first group: all haplotypes matching first target but not the second
    // second group: all haplotypes matching the second target
    df = 0;
    vector<Haplotype> partialCode;
    for (int i = 0; i < genoCode.size(); i++) if (!zero[i]) {
            Haplotype h = genoCode[i];
            Haplotype testhap;
            testhap.assign(h.begin(), h.begin() + conditionsize + testsize);
            // if rare, or not matching either target, put in own group
            if (rare[i] ||
                    !weq(testhap, compare1) && !weq(testhap, compare2)) {
                testhap.enter(partialCode);
                group[0][i] = testhap.index(partialCode);
            } else {
                if (weq(testhap, compare1) && !weq(testhap, compare2)) {
                    compare1.enter(partialCode);
                    group[0][i] = compare1.index(partialCode);
                } else {
                    compare2.enter(partialCode);
                    group[0][i] = compare2.index(partialCode);
                }
            }
        }

    df = partialCode.size() - 1;
	if (typeOfPhenotype == "polytomous")  df *= K-1;

    return(getloglikelihood("Alternative", options, which, false));
}



//  individual tests

double UnphasedAnalysis::individualtests(UnphasedOptions &options, const string &which) {

    // individual score tests calculated at beta=0
    // allowing for estimation of the nuisance parameters

    // maximise likelihood under global null
    int nhap = genoCode.size();
    if (which == "asymptotic") {
        *outStream << "Individual tests: " << flush;
    }
    for (int i = 0; i < group.size(); i++)
        for (int j = 0; j < group[i].size(); j++) {
            group[i][j] = 0;
        }
    for (int i = 0; i < beta.size(); i++) {
        beta[i] = 0;
    }
    int ndf = 0;
    double null = fullNull(options, which, ndf); //getloglikelihood("Null",options,which,true);

    // get score and information at the maximum
    vector<double> gradient(sizeOfGradient(options), 0);
    double llhd = 0;
    score(options, llhd, gradient);

    vector<vector<double> >variance(gradient.size());
    for (int i = 0; i < variance.size(); i++) {
        variance[i].resize(variance.size(), 0);
    }
    numericalHessian(options, variance, gradient, llhd);

    // set up group to reflect tag markers
    vector<Haplotype> partialCode;
    int conditionsize = options.condition.size() * (options.condgenotype ? 2 : 1);
    int testsize = options.window * (options.genotype ? 2 : 1);
    for (int i = 0; i < nhap; i++) if (!zero[i]) {
            Haplotype h = genoCode[i];
            Haplotype testhap;
            testhap.assign(h.begin(), h.begin() + conditionsize + testsize);
            testhap.enter(partialCode);
            group[0][i] = testhap.index(partialCode);
        }

    // calculate score allowing for tag markers
    vector<double> workingGradient(gradient.size(), 0);
    for (int i = 0; i < nhap; i++) if (!zero[i]) {
            workingGradient[group[0][i]] += gradient[i];
        }
    for (int i = nhap; i < gradient.size(); i++) {
        workingGradient[i] = gradient[i];
    }

    // calculate variance allowing for tag markers
    vector<vector<double> > workingVariance(variance.size());
    for (int i = 0; i < variance.size(); i++) {
        workingVariance[i].resize(variance.size(), 0);
    }
    for (int i = 0; i < nhap; i++) if (!zero[i]) {
            for (int j = 0; j < nhap; j++) if (!zero[j]) {
                    workingVariance[group[0][i]][group[0][j]] += variance[i][j];
                }
            for (int j = nhap; j < variance.size(); j++) {
                workingVariance[group[0][i]][j] += variance[i][j];
                workingVariance[j][group[0][i]] += variance[j][i];
            }
        }
    for (int i = nhap; i < variance.size(); i++)
        for (int j = nhap; j < variance.size(); j++) {
            workingVariance[i][j] = variance[i][j];
        }

    // adjust variance of betas to reflect estimation of nuisance parameters
    vector<vector<double> > Vaa, Vab, Vbb;

    // variance-covariance matrix of the nuisance parameters only
    Vaa.resize(gradient.size() - nhap);
    for (int i = 0; i < Vaa.size(); i++) {
        Vaa[i].resize(Vaa.size(), 0);
        for (int j = 0; j < Vaa.size(); j++) {
            Vaa[i][j] = workingVariance[i+nhap][j+nhap];
        }
    }

    // covariance between betas and nuisance parameters
    Vab.resize(nhap);
    for (int i = 0; i < nhap; i++) {
        Vab[i].resize(Vaa.size(), 0);
        for (int j = 0; j < Vaa.size(); j++) {
            Vab[group[0][i]][j] = workingVariance[group[0][i]][j+nhap];
        }
    }

    // obtain variance-covariance matrix of betas
    svdinvert(Vaa);
    sandwich(Vab, Vaa, Vbb);

    // the score tests
    for (int i = 0; i < nhap; i++) {
        if (!zero[i] && !rare[i]) {
            chisq[i] = workingGradient[group[0][i]] * workingGradient[group[0][i]] /
                       (workingVariance[group[0][i]][group[0][i]] - Vbb[group[0][i]][group[0][i]]);
            pvalue[i] = pochisq(chisq[i], 1);
            bestpvalue = min(bestpvalue, pvalue[i]);
            multipleTests++;
        } else {
            chisq[i] = 0;
            pvalue[i] = 1;
        }
    }

    // return the null likelihood obtained in computing the score tests
    return(null);
}



