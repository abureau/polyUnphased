/* UnphasedAnalysis.h - Association analysis in trios and unrelateds

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

#include <fstream>
#include "stl.h"
#include "NuclearFamily.h"
#include "Haplotype.h"

#ifndef __UNPHASEDANALYSIS__

class UnphasedAnalysis: public LinkageData {
public:
    UnphasedAnalysis();
    UnphasedAnalysis(ostream *);
    UnphasedAnalysis(string &, ostream *);
    UnphasedAnalysis(string &, string &, ostream *);
    UnphasedAnalysis(string &, string &, string &, ostream *);
    UnphasedAnalysis(string &, string &, string &, string &, ostream *);
    ~UnphasedAnalysis();
    void run(class UnphasedOptions &);
protected:
    void outputOptions(class UnphasedOptions &);
    void outputResults(vector<int> &, string &, double, double, int, class UnphasedOptions &);
    void outputTabularHeaders(class UnphasedOptions &);
    void outputTabular(class UnphasedOptions &, vector<int>&, double, double, int);
    void dumpHaplotypes(class UnphasedOptions &, ofstream &);
    Locus getLocus(string &);
    void checkMarkerNames(vector<string> &);
    void checkTraitNames(vector<string> &);
    void outputLD(class UnphasedOptions &);
    void score(class UnphasedOptions &, double &, bool);
    void score(class UnphasedOptions &, double &, vector<double> &);
    void scoreFamily(class NuclearFamily &, int, vector<double> &, vector<int> &,
                     vector<vector<double> > &,
                     class UnphasedOptions &, double &, vector<double> &,
                     bool *, double *, vector<double> *);
    void scoreUnrelated(class Subject &, int, double &, vector<double> &,
                        class UnphasedOptions &, double &, vector<double> &,
                        bool *, double *, vector<double> *);
    void numericalHessian(UnphasedOptions &, vector<vector<double> > &, vector<double> &, double);
    void clear1();
    void resizeArrays();
    void readalleles(UnphasedOptions &);
    void readLevels(UnphasedOptions &);
    void readBaselines(UnphasedOptions &);
    int informativesibs(bool [], NuclearFamily &, UnphasedOptions &);
    bool genotyped(Subject &, int);
    bool missingParentalGenotypes(NuclearFamily &, UnphasedOptions &, int);
    bool missingSubjectGenotypes(Subject &, int);
    bool ambiguous(class NuclearFamily &, int, class UnphasedOptions &);
    bool ambiguous(class Subject &, class UnphasedOptions &);
    void makepvalues();
    int parametersPerHaplotype(UnphasedOptions &);
    int sizeOfGradient(UnphasedOptions &);
    int freeParameters(class UnphasedOptions &);
    void gradientFreeParameters(vector<double> &, vector<double> &,
                                class UnphasedOptions &);
    double amotry(vector<vector<double> > &, vector<double> &, vector<double> &,
                  int, double, class UnphasedOptions &, int &);
    double amoeba(class UnphasedOptions &, bool, int &);
    double evaluate(vector<double> &, vector<double> &,
                    class UnphasedOptions &, int &);
    double lnsrch(vector<double> &, double, vector<double> &,
                  vector<double> &, vector<double> &,
                  double, int &, UnphasedOptions &, int &);
    double mydfpmin(UnphasedOptions &, bool, int &);
    double getloglikelihood(const string &, class UnphasedOptions &, const string &, bool);
    void findUsableSubjects(class UnphasedOptions &);
    void listConsistent(class UnphasedOptions &);
    void countFamilyPhenotypes();
    void listHaplotypesFamily(vector<int> &,
                              Haplotype &, Haplotype &, Haplotype &, Haplotype &,
                              NuclearFamily &, int, UnphasedOptions &, bool &, int);
    void listHaplotypesUnrelated(vector<int> &, Haplotype &, Haplotype &, Subject &,
                                 UnphasedOptions &, int);
    void listConsistentAllSibsFamily(Haplotype &, Haplotype &, Haplotype &, Haplotype &,
                                     vector<vector<int> > &, vector<bool> &,
                                     vector<int> &,
                                     class UnphasedOptions &, int, int);
    void listConsistentAllSibsSibship(Haplotype &, Haplotype &,
                                      vector<vector<int> > &, vector<bool> &,
                                      vector<int> &,
                                      class UnphasedOptions &, int);
    void registerHap(class Haplotype &);
    void registerGeno(class Haplotype &);
    bool consistentGenotypes(Haplotype &, Haplotype &, Haplotype &, Haplotype &,
                             vector<bool> &, vector<bool> &, vector<bool> &,
                             UnphasedOptions &, NuclearFamily &);
    bool consistentGenotypes(Haplotype &, Haplotype &, class UnphasedOptions &, class Subject &);
    void exploratory(class UnphasedOptions &, const string &);

    double individualtests(class UnphasedOptions &, const string &);

    double fullNull(UnphasedOptions &, const string &, int &);
    double fullAlternative(UnphasedOptions &, const string &, int &);
    double haploMainAlternative(UnphasedOptions &, const string &, int &);
    double alleleMainAlternative(UnphasedOptions &, const string &, int &);
    double gxgNull(UnphasedOptions &, const string &, int &);
    double gxgAlternative(UnphasedOptions &, const string &, int &);
    double pairwiseNull(UnphasedOptions &, const string &, int &);
    double pairwiseAlternative(UnphasedOptions &, const string &, int &);

    void permutedata(const string &);
    void listMarkerCombinations(vector<vector<int> > &, vector<int>, int, UnphasedOptions &, int);
    void analysetraits(class UnphasedOptions &, const string &, ofstream &);
    void analysemarkers(class UnphasedOptions &, const string &, string &, ofstream &);
    void analysemarkers_main(class UnphasedOptions &, const string &, string &, ofstream &);

    bool haveBetaParent(class UnphasedOptions &);
    bool haveBetaParent0(class UnphasedOptions &);
    bool haveAlpha(class UnphasedOptions &);
    bool haveAlpha0(class UnphasedOptions &);

    //  output files
    ofstream tabularUnrelateds;
    ofstream tabularFamilies;

    // pedigree file data
    FamilyList familylist;
    vector<Subject> unrelateds;
    vector<bool> permute;
    vector<int> permuteUnrelated;
    vector<Haplotype> haploCode, genoCode;
    vector<vector<int> > consistentHaps;
    vector<vector<double> > realisationProb;
    vector<bool> homix;

    // stores whether a covariate is a confounder or an effect modifier
    vector<bool> confounder;
    // stores whether a covariate is factorial or not
    vector<bool> factor;
    // stores the levels of each factorial covariate
    vector<set<double> > factorLevels;
    // stores the baseline level of each factorial covariate
    vector<double> baseline;

    // haplotype data
    valarray<double> familyCount[2], unrelatedCount[2],
             frequency, betaparent, alpha, beta,
             chisq, pvalue, stderror;
    double betaparent0, alpha0; // intercepts
    valarray<bool> rare, zero;

    // covariate data
    // set of effects for each level of each variable
    vector<vector<valarray<double> > > betaCovariate, betaparentCovariate;
    vector<vector<double> > betaCovariate0, betaparentCovariate0; // intercepts
    // standard errors of covariate effects relative to beta
    vector<vector<valarray<double> > > stderrorCovariate;
    // covariate names
    vector<string> covariateName;

    // design matrix (in a non standard format)
    vector<valarray<int> > group;

    // bookkeeping data
    vector<int> currentmarker;
    int currentphenotype, currentphenotype2;
    vector<vector<int> > allele;
    double globalpvalue;
    double bestpvalue;
    int multipleTests;
    vector<double> empiricalDistribution;
    string typeOfPhenotype;
    Haplotype reference;

    bool haveFamilies;
    bool haveUnrelateds;

    vector<vector<double> > empiricalVariance;

};

#endif
#define __UNPHASEDANALYSIS__
