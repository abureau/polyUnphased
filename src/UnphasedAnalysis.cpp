/* UnphasedAnalysis.cpp - Association analysis in trios and unrelateds

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

// system includes
#include <iostream>
#include <strstream>
#include <fstream>
#include <cmath>
#include <string>
#include <sys/time.h>
#include <ctime>

// local includes
#include "stl.h"
#include "UnphasedOptions.h"
#include "UnphasedAnalysis.h"
#include "pvalues.h"
#include "asran.h"

static Random ran;

void debug(vector<double> t) {
    for (int i = 0; i < t.size(); i++) {
        cout << t[i] << " ";
    }
    cout << endl;
}

void debug(vector<int> t) {
    for (int i = 0; i < t.size(); i++) {
        cout << t[i] << " ";
    }
    cout << endl;
}

void debug(vector<vector<double> >&a) {
    for (int i = 0; i < a.size(); i++) {
        for (int j = 0; j < a.size(); j++) {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}

void brk(int i, int j) {
    cout << "breakpoint " << i << " " << j << endl;
}

UnphasedAnalysis::UnphasedAnalysis(ostream *os = &cout): LinkageData(os) {}

UnphasedAnalysis::~UnphasedAnalysis() {}

//  haveBetaParent

bool UnphasedAnalysis::haveBetaParent(UnphasedOptions &options) {
    return(haveFamilies && options.parentrisk && !options.hhrr);
}



//  haveBetaParent0

bool UnphasedAnalysis::haveBetaParent0(UnphasedOptions &options) {
    return(haveBetaParent(options) && typeOfPhenotype == "quant" && options.uncentred);
}



//  haveAlpha

bool UnphasedAnalysis::haveAlpha(UnphasedOptions &options) {
    return(typeOfPhenotype == "quant" && (haveFamilies || betaCovariate.size()));
}



//  haveAlpha0

bool UnphasedAnalysis::haveAlpha0(UnphasedOptions &options) {
    return(haveAlpha(options) && betaCovariate.size());
}



//  getLocus

Locus UnphasedAnalysis::getLocus(string &marker) {
    vector<Locus>::iterator j = locus.begin();
    while (j != locus.end() && (j->name != marker || j->type != 3)) {
        j++;
    }
    if (j == locus.end()) {
        *outStream << endl << "ERROR: marker " << marker << " not found\n";
        exit(-1);
    }
    return *j;
}

//  checkMarkerNames

void UnphasedAnalysis::checkMarkerNames(vector<string> &marker) {
    for (vector<string>::iterator i = marker.begin(); i != marker.end(); i++) {
        vector<Locus>::iterator j = locus.begin();
        while (j != locus.end() && (j->name != *i || j->type != 3)) {
            j++;
        }
        if (j == locus.end()) {
            *outStream << endl << "ERROR: marker " << *i << " not found\n";
            exit(-1);
        }
    }
}



//  checkTraitNames

// check phenotypes exist in the data file
void UnphasedAnalysis::checkTraitNames(vector<string> &trait) {
    for (int i = 0; i < trait.size(); i++) {
        if (trait[i] != "sibsex" && trait[i] != "parsex") {
            int j = 0;
            while (j < locus.size() && !(locus[j].type % 4 == 0 && locus[j].name == trait[i])) {
                j++;
            }

            if (j == locus.size()) {
                *outStream << endl << "ERROR: trait " << trait[i] << " not found\n";
                exit(-1);
            }
        }
    }
}



//  clear

// clear all haplotype data
void UnphasedAnalysis::clear1() {
    haploCode.resize(0);
    genoCode.resize(0);
    realisationProb.resize(0);
    homix.resize(0);
    chisq.resize(0);
    pvalue.resize(0);
    rare.resize(0);
    zero.resize(0);
    for (int i = 0; i < 8; i++)
    	familyCount[i].resize(0);
    unrelatedCount[0].resize(0);
    unrelatedCount[1].resize(0);
    frequency.resize(0);
    betaparent.resize(0);
    alpha.resize(0);
    beta.resize(0);
    stderror.resize(0);
    for (int i = 0; i < betaCovariate.size(); i++) {
        for (int j = 0; j < betaCovariate[i].size(); j++) {
            betaCovariate[i][j].resize(0);
            betaparentCovariate[i][j].resize(0);
            stderrorCovariate[i][j].resize(0);
        }
    }
    for (int i = 0; i < group.size(); i++) {
        group[i].resize(0);
    }
    reference.clear();
}



//  resizeArrays

void UnphasedAnalysis::resizeArrays() {
    int size = genoCode.size();
    for (int i = 0; i < 8; i++) {
        familyCount[i].resize(size);
        familyCount[i] = 0;
    }
    for (int i = 0; i < 2; i++) {
        unrelatedCount[i].resize(size);
        unrelatedCount[i] = 0;
    }
    frequency.resize(size);
    frequency = 0;
    betaparent.resize(size*(K-1));
    betaparent = 0;
    alpha.resize(size);
    alpha = 0;
    beta.resize(size*(K-1));
    beta = 0;
    stderror.resize(size*(K-1));
    stderror = 0;
    for (int i = 0; i < betaCovariate.size(); i++)
        for (int j = 0; j < betaCovariate[i].size(); j++) {
            betaCovariate[i][j].resize(size);
            betaCovariate[i][j] = 0;
            betaparentCovariate[i][j].resize(size);
            betaparentCovariate[i][j] = 0;
            stderrorCovariate[i][j].resize(size);
            stderrorCovariate[i][j] = 0;
        }

    chisq.resize(size);
    chisq = 0;
    pvalue.resize(size);
    pvalue = 0;
    rare.resize(size);
    rare = false;
    zero.resize(size);
    zero = false;
}



//  readalleles

void UnphasedAnalysis::readalleles(UnphasedOptions &options) {
    allele.resize(options.window + options.condition.size() + options.tag.size());
    for (int k = 0; k < allele.size(); k++) {
        allele[k].clear();
        // set of alleles at this marker
        set<int> currentallele;
        for (PedigreeSet::iterator i = pedigree.begin();
                i != pedigree.end(); i++) {
            for (vector<Subject>::iterator j = i->second.begin();
                    j != i->second.end(); j++) {
                //  currentallele.insert(j->marker[currentmarker[k]][0]);
                //  currentallele.insert(j->marker[currentmarker[k]][1]);
                currentallele.insert(j->marker[currentmarker[k]] / 16);
                currentallele.insert(j->marker[currentmarker[k]] % 16);
            }
        }
        currentallele.erase(0);
        for (set<int>::iterator i = currentallele.begin(); i != currentallele.end(); i++) {
            allele[k].push_back(*i);
        }
    }
}



//  readLevels

// read the levels of factorial covariates and resize elements of betaCovariate
void UnphasedAnalysis::readLevels(UnphasedOptions &options) {
    int nconfounder = options.confounder.size();
    int nmodifier = options.modifier.size();
    factorLevels.resize(factor.size());
    for (int i = 0; i < nconfounder + nmodifier; i++) {
        factorLevels[i].clear();
        if (factor[i]) {
            int traitIndex = traithash[covariateName[i]];
            for (PedigreeSet::iterator j = pedigree.begin();
                    j != pedigree.end(); j++) {
                for (vector<Subject>::iterator k = j->second.begin();
                        k != j->second.end(); k++) {
                    factorLevels[i].insert(k->trait[traitIndex]);
                }
            }
            betaCovariate[i].resize(factorLevels[i].size() - 1);
            betaparentCovariate[i].resize(factorLevels[i].size() - 1);
            betaCovariate0[i].resize(factorLevels[i].size() - 1);
            betaparentCovariate0[i].resize(factorLevels[i].size() - 1);
            stderrorCovariate[i].resize(factorLevels[i].size() - 1);
        }
        if (!factor[i]) {
            factorLevels[i].insert(0);
            betaCovariate[i].resize(1);
            betaparentCovariate[i].resize(1);
            betaCovariate0[i].resize(1);
            betaparentCovariate0[i].resize(1);
            stderrorCovariate[i].resize(1);
        }
    }
}



//  readBaselines

// set the baselines of factorial covariates
void UnphasedAnalysis::readBaselines(UnphasedOptions &options) {
    int nconfounder = options.confounder.size();
    int nmodifier = options.modifier.size();
    for (int i = 0; i < nconfounder + nmodifier; i++) {
        if (factor[i]) {
            // if no baseline given, default to first level
            if (i >= options.baseline.size() ||
                    factorLevels[i].find(options.baseline[i]) == factorLevels[i].end()) {
                baseline[i] = *factorLevels[i].begin();
            } else {
                baseline[i] = options.baseline[i];
            }
        } else {
            if (covariateName[i] == "parsex" || covariateName[i] == "sibsex") {
                if (i >= options.baseline.size()) {
                    baseline[i] = FEMALE;
                } else {
                    if (options.baseline[i] != 1 && options.baseline[i] != 2) {
                        *outStream << endl << "ERROR: baseline level for " << covariateName[i] << " must be 1 or 2" << endl;
                        exit(-1);
                    }
                    baseline[i] = options.baseline[i];
                }
            } else {
                baseline[i] = 0;
            }
        }
    }
}



//  permutedata

void UnphasedAnalysis::permutedata(const string &which) {
    permute.resize(familylist.size());
    for (int i = 0; i < permute.size(); i++) {
        permute[i] = (which == "asymptotic" || ran.rand() < 0.5);
    }
    permuteUnrelated.resize(unrelateds.size());
    for (int i = 0; i < permuteUnrelated.size(); i++) {
        permuteUnrelated[i] = i;
    }
    if (which != "asymptotic") {
        // random shuffle of the labels
        for (int i = 0; i < permuteUnrelated.size(); i++) {
            // swap with randomly chosen element
            int ix = (int)(ran.rand() * (permuteUnrelated.size() - i)) + i;
            int tmp = permuteUnrelated[ix];
            permuteUnrelated[ix] = permuteUnrelated[i];
            permuteUnrelated[i] = tmp;
        }
    }
}



//  listMarkerCombinations

void UnphasedAnalysis::listMarkerCombinations(vector<vector<int> > &markerCombination, vector<int> combination, int window, UnphasedOptions &options, int n) {
    if (n < window) {
        for (int i = (n == 0 ? 0 : (combination[n-1] + 1)); i < options.marker.size(); i++) {
            combination.push_back(i);
            listMarkerCombinations(markerCombination, combination, window, options, n + 1);
            combination.pop_back();
            if (n > 0 && !options.allcombinations) {
                i = options.marker.size();
            }
        }
    } else {
        markerCombination.push_back(combination);
    }
}



//  analysemarkers

void UnphasedAnalysis::analysemarkers(UnphasedOptions &options, const string &which, string &trait, ofstream &dumpfile) {

    string filename = options.listMarkerFilename;
    if (filename != "") {
        ifstream infile(filename.c_str());
        if (!infile) {
            *outStream << "Failed to open list marker file " << filename << endl;
            exit(-1);
        }
        *outStream << "Reading list marker file " << filename << "..." << endl;

        string line;
        while ((line = getline(infile)) != "") {
            options.condition.clear();
            options.marker.clear();
            options.tag.clear();
            options.window = 1;
            istrstream instr(line.c_str());
            vector<string> words;
            char buf[line.length()];
            words.push_back("markerlist");
            while (instr >> buf, !instr.fail()) {
                string word = buf;
                words.push_back(word);
            }
            char **words_c = new char *[words.size()];
            for (int i = 0; i < words.size(); i++) {
                words_c[i] = new char[words[i].size() + 1];
                words[i].copy(words_c[i], words[i].size());
                words_c[i][words[i].size()] = 0;
            }
            UnphasedOptions tmpOptions(words.size(), words_c);
            options.condition = tmpOptions.condition;
            options.marker = tmpOptions.marker;
            options.tag = tmpOptions.tag;
            options.window = tmpOptions.window;
            analysemarkers_main(options, which, trait, dumpfile);
        }
    } else {
        analysemarkers_main(options, which, trait, dumpfile);
    }
}


// main analysis routine
// loops through selected markers and analyse each one
void UnphasedAnalysis::analysemarkers_main(UnphasedOptions &options, const string &which, string &trait, ofstream &dumpfile) {

  // open Tabular output files for this trait
  if (options.tabularFilename != "") {
        string filename;
        filename = options.tabularFilename + "." + trait + ".unrelated.out";
        tabularUnrelateds.open(filename.c_str());
        if (!tabularUnrelateds.is_open()) {
            *outStream << "ERROR: could not open unrelated tabular file " << filename << endl;
            exit(-1);
        }
        filename = options.tabularFilename + "." + trait + ".family.out";
        tabularFamilies.open(filename.c_str());
        if (!tabularFamilies.is_open()) {
            *outStream << "ERROR: could not open family tabular file " << filename << endl;
            exit(-1);
        }
	outputTabularHeaders(options);
        tabularUnrelateds.precision(4);
        tabularFamilies.precision(4);
  }

    // loop through window sizes
    for (int window = (options.allwindows ? 1 : options.window);
            window <= (options.allwindows ? options.marker.size() : options.window);
            window++) {

        vector<vector<int> > markerCombination;
        vector<int> combination;
        listMarkerCombinations(markerCombination, combination, window, options, 0);
        options.window = window;

        // fix conditioning markers
        currentmarker.resize(options.condition.size() + options.window + options.tag.size());
        for (int i = 0; i < options.condition.size(); i++) {
            currentmarker[i] = markerhash[options.condition[i]];
        }

        // fix tag markers
        for (int j = 0; j < options.tag.size(); j++) {
            currentmarker[options.condition.size()+options.window+j] = markerhash[options.tag[j]];
        }

        // loop through test marker combinations
        for (int i = 0; i < markerCombination.size(); i++) {

            // store current markers
            for (int j = 0; j < options.window; j++) {
                currentmarker[options.condition.size()+j] = markerhash[options.marker[markerCombination[i][j]]];
            }

            clear1();
            readalleles(options);

            // exploratory pass to make list of consistent haplotypes
            // and identify rare haplotypes
            exploratory(options, which);

            // null and alternative likelihoods
            double null, alternative;
            int df = 0, ndf = 0;

            // individual haplotype tests
            if (options.individual) {
                null = individualtests(options, which);
            }

            // get the null and alternative likelihoods
            // according to the type of analysis performed
            if (options.compare1.size() && options.compare2.size()) {
                null = pairwiseNull(options, which, ndf);
                alternative = pairwiseAlternative(options, which, df);
            } else {
                if (options.model == "full" || options.model == "") {
                    if (!options.individual || options.condition.size() != 0 || options.rare != 0) {
                        null = fullNull(options, which, ndf);
                    }
                    alternative = fullAlternative(options, which, df);
                }
                if (options.model == "haplomain") {
                    if (!options.individual || options.condition.size() != 0 || options.rare != 0) {
                        null = fullNull(options, which, ndf);
                    }
                    alternative = haploMainAlternative(options, which, df);
                }
                if (options.model == "allelemain" ||
                        options.model == "commonmain") {
                    if (!options.individual || options.condition.size() != 0 || options.rare != 0) {
                        null = fullNull(options, which, ndf);
                    }
                    alternative = alleleMainAlternative(options, which, df);
                }
                if (options.model == "gxg") {
                    null = gxgNull(options, which, ndf);
                    alternative = gxgAlternative(options, which, df);
                }
                if (options.model == "null") {
                    if (!options.individual || options.rare != 0) {
                        null = fullNull(options, which, ndf);
                    }
                    alternative = null;
                    df = ndf;
                }
            }
            df -= ndf;
            if (!options.individual) {
                double thispvalue = pochisq(2 * (alternative - null), df);
                bestpvalue = min(bestpvalue, thispvalue);
                multipleTests++;
            }

            //output results
            if (which == "asymptotic" || options.outputPermutation) {
                outputResults(markerCombination[i], trait, null, alternative, df, options);
                if (options.LD) {
                    outputLD(options);
                }
                if (which == "asymptotic" && options.dumpFilename != "") {
                    dumpHaplotypes(options, dumpfile);
                }
            }
            if (options.tabularFilename != "") {
                outputTabular(options, markerCombination[i], null, alternative, df);
            }

        } // for i
    } // for window

    // close tabularfile
     if (tabularFamilies.is_open()) {
         tabularFamilies.close();
     }
     if (tabularUnrelateds.is_open()) {
         tabularUnrelateds.close();
     }

}



//  analysetraits

void UnphasedAnalysis::analysetraits(UnphasedOptions &options, const string &which, ofstream &dumpfile) {

	string pheno1, pheno2;
    // loop over diseases, joint affection statuses and traits
    for (int diseasecount = 0; diseasecount < options.disease.size(); diseasecount++) {
        currentphenotype = diseasehash[options.disease[diseasecount]];
        typeOfPhenotype = "binary";
        analysemarkers(options, which, options.disease[diseasecount], dumpfile);
    }
    for (int polycount = 0; polycount < options.polytomous.size(); polycount++) {
        currentphenotype = diseasehash[options.polytomous[polycount]];
        typeOfPhenotype = "polytomous";
        polypheno = true;
        K = Kvec[currentphenotype];
        analysemarkers(options, which, options.polytomous[polycount], dumpfile);
        }    
    for (int jointcount = 0; jointcount < options.joint.size(); jointcount++) {
        // Extraire les noms de statut d'affection de options.joint[jointcount]
        // Nom du premier statut d'affection, avant le ":"
	std::size_t seppos = options.joint[jointcount].find(":");
        pheno1 = options.joint[jointcount].substr(0,seppos);
        // Nom du deuxième statut d'affection, après le ":" 
        pheno2 = options.joint[jointcount].substr(seppos+1,options.joint[jointcount].size()-seppos);
			// Additional debugging code
        /*    if (options.llhd) {
            	cout << "pheno1 " << pheno1;
            	cout << " pheno2 " << pheno2 << endl;
                } */
        currentphenotype = diseasehash[pheno1];
        currentphenotype2 = diseasehash[pheno2];
        typeOfPhenotype = "polytomous";
        polypheno = false;
        K = 4;
        analysemarkers(options, which, options.joint[jointcount], dumpfile);
    }
    for (int traitcount = 0; traitcount < options.trait.size(); traitcount++) {
        currentphenotype = traithash[options.trait[traitcount]];
        typeOfPhenotype = "quant";
        analysemarkers(options, which, options.trait[traitcount], dumpfile);
    }
    if (options.disease.size() == 0 && options.trait.size() == 0 && options.joint.size() == 0 && options.polytomous.size() == 0) {
        // default analysis is the first disease
        // if none found, the first trait
        int i = 0;
        while (i < locus.size() && locus[i].type != 1) {
            i++;
        }
        if (i == locus.size()) {
            int j = 0;
            while (j < locus.size() && locus[j].type % 4 != 0) {
                j++;
            }
            if (j == locus.size()) {
                *outStream << endl << "ERROR: no affection status or quantitative trait found\n";
                exit(-1);
            }
            currentphenotype = traithash[locus[j].name];
            typeOfPhenotype = "quant";
            analysemarkers(options, which, locus[j].name, dumpfile);
        } else {
            currentphenotype = diseasehash[locus[i].name];
            typeOfPhenotype = "binary";
            analysemarkers(options, which, locus[i].name, dumpfile);
        }
    }
}



//  run

// run the analysis
void UnphasedAnalysis::run(UnphasedOptions &options) {

    // open the dumpfile
    ofstream dumpfile;
    if (options.dumpFilename != "") {
        dumpfile.open(options.dumpFilename.c_str());
        if (!dumpfile.is_open()) {
            *outStream << "ERROR: could not open dump file " << options.dumpFilename << endl;
            exit(-1);
        }
    }

    // timing the calculation
    struct timeval starttime, endtime;
    gettimeofday(&starttime, NULL);

    // loop through the pedigree files
    for (int file = 0; file < options.inputfile.size(); file++) {
        // read input files
        if (options.datafile != "") {
            readdatafile(options.datafile);
        } else {
            defaultdatafile(options.inputfile[file]);
        }
        if (options.mapfile != "") {
            readmapfile(options.mapfile);
        }
        if (options.bedfile == "") {
            readpedfile(options.inputfile[file]);
        } else {
            readpedfile(options.inputfile[file], options.bedfile);
        }
        if (options.phenofile != "") {
            readphenofile(options.phenofile);
        }
        makehashtables();

        checkMarkerNames(options.tag);
        checkMarkerNames(options.condition);
        checkMarkerNames(options.marker);
        checkTraitNames(options.trait);
        checkTraitNames(options.modifier);
        checkTraitNames(options.confounder);
        checkTraitNames(options.factor);

        // initialise random number seed
        ran.init(options.randomseed, 1, 1);

        // set up lists of covariates
        int nconfounder = options.confounder.size();
        int nmodifier = options.modifier.size();
        betaCovariate.resize(nconfounder + 1 * nmodifier); // was 2
        betaparentCovariate.resize(betaCovariate.size());
        betaCovariate0.resize(nconfounder + 1 * nmodifier);
        betaparentCovariate0.resize(betaCovariate.size());
        stderrorCovariate.resize(betaCovariate.size());
        confounder.resize(nconfounder + nmodifier, false);
        baseline.resize(nconfounder + nmodifier, 0);
        for (int i = 0; i < nconfounder; i++) {
            confounder[i] = true;
        }
        covariateName.clear();
        for (int i = 0; i < nconfounder; i++) {
            covariateName.push_back(options.confounder[i]);
        }
        for (int i = 0; i < nmodifier; i++) {
            covariateName.push_back(options.modifier[i]);
        }
        factor.resize(nconfounder + nmodifier, false);
        for (int i = 0; i < factor.size(); i++) {
            int j = 0;
            while (j < options.factor.size() && covariateName[i] != options.factor[j]) {
                j++;
            }
            factor[i] = (j < options.factor.size() && covariateName[i] != "sibsex" && covariateName[i] != "parsex");
            if (covariateName[i] == "sibsex" || covariateName[i] == "parsex") {
                baseline[i] = FEMALE;
            }
        }

        readLevels(options);
        readBaselines(options);

        // if no marker names given, do them all
        if (options.marker.size() == 0) {
            for (vector<Locus>::iterator i = locus.begin(); i != locus.end(); i++)
                if (i->type == 3) {
                    options.marker.push_back(i->name);
                }
        }

        // break pedigrees into nuclear families and report errors
        familylist.extract(pedigree);
        for (FamilyList::iterator j = familylist.begin(); j != familylist.end(); j++) {
            j->report();
        }
        // get the single unrelated subjects
        unrelateds.clear();
        for (PedigreeSet::iterator j = pedigree.begin(); j != pedigree.end(); j++) {
            if (j->second.size() == 1 && j->second[0].pid == "0" && j->second[0].mid == "0") {
                Subject subject = j->second[0];
                subject.id = j->first;
                if (!options.chrY || subject.sex == MALE) {
                    unrelateds.push_back(subject);
                }
            }
        }
        // on chrY, take first male sib from each nuclear family as an unrelated
        if (options.chrY) {
            for (FamilyList::iterator family = familylist.begin();
                    family != familylist.end(); family++)
                if (family->status == "") {
                    bool haveSib = false;
                    for (int i = 0; i < family->sibs.size() && !haveSib; i++)
                        if (family->sibs[i].sex == MALE) {
                            unrelateds.push_back(family->sibs[i]);
                            haveSib = true;
                        }
                }
            familylist.clear();
        }

        // note which parameters will be in the linear models
        haveFamilies = familylist.size() > 0;
        haveUnrelateds = unrelateds.size() > 0;

        // main analysis
        if (!options.briefOutput) {
            outputOptions(options);
        }

        permutedata("asymptotic");
        bestpvalue = 1;
        multipleTests = 0;
        empiricalDistribution.clear();
        analysetraits(options, "asymptotic", dumpfile);
        empiricalDistribution.push_back(bestpvalue);

        // global permutation test
        if (options.npermutation) {
            double bestp_orig = bestpvalue;
            // truncate bestp_orig to 8 decimal places
            //bestp_orig=(int)(bestp_orig*1e8)/1e8;
            //cout << "bestp_orig " << bestp_orig <<endl;
            globalpvalue = 1;
            *outStream << "Running " << options.npermutation << " permutations..." << endl;
            for (int i = 0; i < options.npermutation; i++) {
                //if (!options.outputPermutation) *outStream << i+1 << flush;
                permutedata("permutation");
                bestpvalue = 1;
                multipleTests = 0;
                if (options.outputPermutation) {
                    *outStream << "Permutation #" << i + 1 << endl;
                }
                analysetraits(options, "permutation", dumpfile);
                // truncate p-value to 8 decimal places
                //bestpvalue=(int)(bestpvalue*1e8)/1e8;
                if (bestpvalue <= bestp_orig) {
                    globalpvalue++;
                }
                empiricalDistribution.push_back(bestpvalue);
                if (!options.outputPermutation) {
                    for (int j = 1; j <= 10; j++)
                        if (i + 1 == options.npermutation * j / 10) {
                            if (j == 1) {
                                *outStream << "Done " << flush;
                            }
                            *outStream << j * 10 << "%" << flush;
                            if (j < 10) {
                                *outStream << "..." << flush;
                            }
                        }
                }
                //if (!options.outputPermutation) *outStream << "\r";
            }
            *outStream << "done" << endl;
            globalpvalue /= options.npermutation + 1;
            *outStream << "Best p-value from " << multipleTests << " tests: " << bestp_orig << endl;
            *outStream << "Adjusted p-value from permutation test: " << globalpvalue << "  standard error: " <<
                       sqrt(globalpvalue*(1.0 - globalpvalue) / options.npermutation) << endl;
            sort(empiricalDistribution.begin(), empiricalDistribution.end());
            *outStream << "Empirical " << options.quantile << "% quantile of the best p-value: " << empiricalDistribution[(int)(empiricalDistribution.size()*options.quantile/100.0)] << endl;
        }

    }

    // close dumpfile
    if (dumpfile.is_open()) {
        dumpfile.close();
    }

    // output timings
    gettimeofday(&endtime, NULL);
    long sex = endtime.tv_sec - starttime.tv_sec;
    long usex = endtime.tv_usec - starttime.tv_usec;
    if (options.time) {
        *outStream << endl << "Elapsed time " << sex + usex / 1000000.0 << " seconds" << endl;
    }

}


