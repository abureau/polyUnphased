/* UnphasedAnalysisExploratory.cpp - Association analysis in trios and unrelateds

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

// local includes
#include "stl.h"
#include "UnphasedOptions.h"
#include "UnphasedAnalysis.h"

inline bool seq(int x, int y) {
    return(x != 0 && y != 0 && x == y);
}
inline bool weq(int x, int y) {
    return(x == 0 || y == 0 || x == y);
}

//  weq Haplotype

static bool weq(const Haplotype &h1, const Haplotype &h2) {
    if (h1.size() != h2.size()) {
        return(false);
    }
    for (int i = 0; i < h1.size(); i++)
        if (!weq(h1[i], h2[i])) {
            return(false);
        }
    return(true);
}



//  genotyped

bool UnphasedAnalysis::genotyped(Subject &subject, int size) {
    bool untyped = true;
    for (int i = 0; i < size; i++)
        //     if (subject.marker[currentmarker[i]][0]!=0 ||
        //  subject.marker[currentmarker[i]][1]!=0) untyped=false;
        if (subject.marker[currentmarker[i]] / 16 != 0 ||
                subject.marker[currentmarker[i]] % 16 != 0) {
            untyped = false;
        }
    return(!untyped);
}



//  possiblegenotype

bool possiblegenotype(Subject &subject, int i, int a1, int a2) {
    //   if (weq(subject.marker[i][0],a1) && weq(subject.marker[i][1],a2) ||
    //       weq(subject.marker[i][0],a2) && weq(subject.marker[i][1],a1))
    if (weq(subject.marker[i] / 16, a1) && weq(subject.marker[i] % 16, a2) ||
            weq(subject.marker[i] / 16, a2) && weq(subject.marker[i] % 16, a1)) {
        return(true);
    }
    return(false);
}



//  determine transmissions in families

static void determinetransmission(NuclearFamily &family, UnphasedOptions &options,
                                  int sib, int marker,
                                  int &fa1, int &fa2, int &ma1, int &ma2,
                                  bool &intercross, bool &homix) {
    int f[2], m[2], s[2];
    //   f[0]=min(family.father.marker[marker][0],family.father.marker[marker][1]);
    //   f[1]=max(family.father.marker[marker][0],family.father.marker[marker][1]);
    //   m[0]=min(family.mother.marker[marker][0],family.mother.marker[marker][1]);
    //   m[1]=max(family.mother.marker[marker][0],family.mother.marker[marker][1]);
    //   s[0]=min(family.sibs[sib].marker[marker][0],family.sibs[sib].marker[marker][1]);
    //   s[1]=max(family.sibs[sib].marker[marker][0],family.sibs[sib].marker[marker][1]);
    f[0] = min(family.father.marker[marker] / 16, family.father.marker[marker] % 16);
    f[1] = max(family.father.marker[marker] / 16, family.father.marker[marker] % 16);
    m[0] = min(family.mother.marker[marker] / 16, family.mother.marker[marker] % 16);
    m[1] = max(family.mother.marker[marker] / 16, family.mother.marker[marker] % 16);
    s[0] = min(family.sibs[sib].marker[marker] / 16, family.sibs[sib].marker[marker] % 16);
    s[1] = max(family.sibs[sib].marker[marker] / 16, family.sibs[sib].marker[marker] % 16);
    // on chrX, treat father as homozygous
    if (options.chrX) {
        f[0] = f[1];
        if (family.sibs[sib].sex == MALE) {
            s[0] = f[1];
        }
    }

    fa1 = 0;
    fa2 = 0;
    ma1 = 0;
    ma2 = 0;

    // intercross
    // 1/2 x 1/2 -> 1/2 or 1/1
    if (seq(f[0], m[0]) && seq(f[1], m[1]) & !weq(f[0], f[1])) {
        // het
        if (seq(s[0], f[0]) && seq(s[1], f[1])) {
            fa1 = f[0];
            fa2 = f[1];
            ma1 = f[1];
            ma2 = f[0];
            intercross = true;
            return;
        }
        // hom
        if (seq(s[0], s[1])) {
            fa1 = s[0];
            fa2 = (fa1 == f[0]) ? f[1] : f[0];
            ma1 = s[1];
            ma2 = (ma1 == m[0]) ? m[1] : m[0];
            homix = true;
            return;
        }
    }

    // eight phases
    // should only be one unique solution
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                if (seq(f[i], s[k]) && seq(m[j], s[1-k])) {
                    fa1 = f[i];
                    fa2 = f[1-i];
                    ma1 = m[j];
                    ma2 = m[1-j];
                }

    // on chrX, force father to be homozygous
    if (options.chrX) {
        fa2 = fa1;
    }
}

#ifdef COMMENT
// assigns 'value' into 'a11' assuming it is transmitted from 'parent1'
static void assign(Subject &parent1, Subject &parent2, Subject &sib, int i,
                   int &a11, int &a12, int &a21, int &a22,
                   int value) {
    // transmitted allele
    asn(a11, value);

    int p11 = parent1.marker[i][0];
    int p12 = parent1.marker[i][1];
    int p21 = parent2.marker[i][0];
    int p22 = parent2.marker[i][1];
    int s1 = sib.marker[i][0];
    int s2 = sib.marker[i][1];

    // nontransmitted allele
    if (seq(p11, value)) {
        asn(a12, p12);
    }
    if (seq(p12, value)) {
        asn(a12, p11);
    }

    // other parent
    if (!weq(s1, value)) {
        asn(a21, s1);
        if (seq(p21, s1)) {
            asn(a22, p22);
        }
        if (seq(p22, s1)) {
            asn(a22, p21);
        }
    }
    if (!weq(s2, value)) {
        asn(a21, s2);
        if (seq(p21, s2)) {
            asn(a22, p22);
        }
        if (seq(p22, s2)) {
            asn(a22, p21);
        }
    }
}

static void determinetransmission(NuclearFamily &family, UnphasedOptions &options,
                                  int k, int i,
                                  int &fa1, int &fa2, int &ma1, int &ma2,
                                  bool &intercross, bool &homix) {
    int f1 = family.father.marker[i][0];
    int f2 = family.father.marker[i][1];
    int m1 = family.mother.marker[i][0];
    int m2 = family.mother.marker[i][1];
    int s1 = family.sibs[k].marker[i][0];
    int s2 = family.sibs[k].marker[i][1];
    //    fa1=s1;
    //    fa2=(f1==fa1?f2:f1);
    //    ma1=s2;
    //    ma2=(m1==ma1?m2:m1);
    //    return;

    // homozygous sib
    // only assign if parents fully typed; avoids Curtis & Sham bias
    if (seq(s1, s2) && f1 && (f2 || options.chrX) && m1 && m2) {
        assign(family.father, family.mother, family.sibs[k], i, fa1, fa2, ma1, ma2, s1);
        assign(family.mother, family.father, family.sibs[k], i, ma1, ma2, fa1, fa2, s1);
    }

    // exclude "semi-type" bias
    // sib must be fully typed, or parents must have no typed allele in common
    if (s1 && (s2 || options.chrX && family.sibs[k].sex == MALE) ||
            !seq(f1, m1) && !seq(f1, m2) && (options.chrX || !seq(f2, m1) &&
                    !seq(f2, m2))) {
        // homozygous parent
        // chrX: pretend father is homozygous
        if ((seq(f1, f2) || options.chrX)) {
            assign(family.father, family.mother, family.sibs[k], i, fa1, fa2, ma1, ma2, f1);
        }
        if (seq(m1, m2)) {
            assign(family.mother, family.father, family.sibs[k], i, ma1, ma2, fa1, fa2, m1);
        }

        // sib allele absent in parent
        if (!weq(s1, m1) && !weq(s1, m2)) {
            assign(family.father, family.mother, family.sibs[k], i, fa1, fa2, ma1, ma2, s1);
        }
        if (!weq(s1, f1) && (!weq(s1, f2) || options.chrX)) {
            assign(family.mother, family.father, family.sibs[k], i, ma1, ma2, fa1, fa2, s1);
        }
        if (!weq(s2, m1) && !weq(s2, m2)) {
            assign(family.father, family.mother, family.sibs[k], i, fa1, fa2, ma1, ma2, s2);
        }
        if (!weq(s2, f1) && (!weq(s2, f2) || options.chrX)) {
            assign(family.mother, family.father, family.sibs[k], i, ma1, ma2, fa1, fa2, s2);
        }
        // chrX: male sib, allele must come from mother
        if (options.chrX && family.sibs[k].sex == MALE) {
            assign(family.mother, family.father, family.sibs[k], i, ma1, ma2, fa1, fa2, s1);
        }
    }

    // intercross cases
    // if parents have the same genotype
    if (seq(min(f1, f2), min(m1, m2)) && seq(max(f1, f2), max(m1, m2)) &&
            // which is heterozygous
            !weq(f1, f2)
            // and we do not distinguish parental sex
            // I don't think we need this any more
            //&& !options.parsex
       ) {
        // if sib has same (het) genotype
        // intercross is set for just this marker
        if (seq(min(f1, f2), min(s1, s2)) && seq(max(f1, f2), max(s1, s2))) {
            fa1 = f1;
            fa2 = f2;
            ma1 = fa2;
            ma2 = fa1;
            intercross = true;
        }
        // if sib is homozygous
        // homix is set across all markers
        if (seq(s1, s2)) {
            homix = true;
        }
    }

    // if chrX, one paternal chromosome is missing
    if (options.chrX && family.sibs[k].sex == MALE) {
        fa2 = 0;
    }

}
#endif



//  findUsableSubjects

void UnphasedAnalysis::findUsableSubjects(UnphasedOptions &options) {
    // families
    for (FamilyList::iterator family = familylist.begin(); family != familylist.end();
            family++)
        if (family->status == "") {
            for (int sib = 0; sib < family->sibs.size(); sib++) {
                bool usesib = (genotyped(family->sibs[sib], options.condition.size() + options.window + options.tag.size()));
                if (typeOfPhenotype == "quant") {
                    usesib = usesib && !family->sibs[sib].missingtrait[currentphenotype];
                } else {
                    usesib = usesib && family->sibs[sib].affection[currentphenotype] != UNSURE;
                    // On considÃ¨re currentphenotype2 seulement si l'option joint est spÃ©cifiÃ©e
                    if (typeOfPhenotype == "polytomous" && currentphenotype2>=0) {
	                    usesib = usesib && family->sibs[sib].affection[currentphenotype2] != UNSURE; 
	                }                       
                }
                for (int conf = 0; conf < options.confounder.size(); conf++)
                    if (options.confounder[conf] != "sibsex" && options.confounder[conf] != "parsex") {
                        usesib = usesib && family->sibs[sib].trait.size() > traithash[options.confounder[conf]] && !family->sibs[sib].missingtrait[traithash[options.confounder[conf]]];
                    }
                for (int modif = 0; modif < options.modifier.size(); modif++)
                    if (options.modifier[modif] != "sibsex" && options.modifier[modif] != "parsex") {
                        usesib = usesib && family->sibs[sib].trait.size() > traithash[options.modifier[modif]] && !family->sibs[sib].missingtrait[traithash[options.modifier[modif]]];
                    }
                family->sibs[sib].usable = usesib;
            }
        }

    // unrelateds
    for (vector<Subject>::iterator subject = unrelateds.begin(); subject != unrelateds.end(); subject++) {
        bool usesubject = genotyped(*subject, options.window + options.condition.size() + options.tag.size());
        if (typeOfPhenotype == "quant") {
            usesubject = usesubject && !unrelateds[permuteUnrelated[subject-unrelateds.begin()]].missingtrait[currentphenotype];
        } else {
            usesubject = usesubject && unrelateds[permuteUnrelated[subject-unrelateds.begin()]].affection[currentphenotype] != UNSURE;
        }
        for (int conf = 0; conf < options.confounder.size(); conf++)
            if (options.confounder[conf] != "sibsex" && options.confounder[conf] != "parsex") {
                usesubject = usesubject && subject->trait.size() > traithash[options.confounder[conf]] && !subject->missingtrait[traithash[options.confounder[conf]]];
            }
        for (int modif = 0; modif < options.modifier.size(); modif++)
            if (options.modifier[modif] != "sibsex" && options.modifier[modif] != "parsex") {
                usesubject = usesubject && subject->trait.size() > traithash[options.modifier[modif]] && !subject->missingtrait[traithash[options.modifier[modif]]];
            }
        subject->usable = usesubject;

    }
}



//  listConsistent

// store a list of haplotypes consistent with each subject
void UnphasedAnalysis::listConsistent(UnphasedOptions &options) {
    consistentHaps.clear();
    haploCode.clear();
    genoCode.clear();

    // loop through families
    for (FamilyList::iterator family = familylist.begin();
            family != familylist.end(); family++)
        if (family->status == "") {
            int size = options.window + options.condition.size() + options.tag.size();
            Haplotype ft(size);
            Haplotype fnt(size);
            Haplotype mt(size);
            Haplotype mnt(size);
            bool homoIntercross = false;

            // the consistent haplotypes for this sibship
            // one row for each sib listing all consistent solutions
            vector<vector<int> > consistentHapsSibship;

            // use sibship model if parents are missing and sibship option used
            // but only informative if there is more than one phenotype in the sibs
            family->sibship = false;
            if (options.sibship && !genotyped(family->father, size) && !genotyped(family->mother, size)) {
                bool virgin = true;
                double firstTrait = 0;
                for (int i = 0; i < family->sibs.size(); i++)
                    if (family->sibs[i].usable) {
                        double trait = (typeOfPhenotype == "quant") ?
                                       family->sibs[i].trait[currentphenotype] :
                                       (double)family->sibs[i].affection[currentphenotype];
                        if (virgin) {
                            firstTrait = trait;
                        } else if (trait != firstTrait) {
                            family->sibship = true;
                        }
                        virgin = false;
                    }
            }

            // loop through the informative sibs
            for (int i = 0; i < family->sibs.size(); i++) {
                vector<int> haps;
                if (family->sibs[i].usable)
                    if (!family->sibship) {
                        listHaplotypesFamily(haps, ft, fnt, mt, mnt, *family, i, options, homoIntercross, 0);
                    } else {
                        listHaplotypesUnrelated(haps, ft, mt, family->sibs[i], options, 0);
                    }
                consistentHapsSibship.push_back(haps);
                if (options.show) {
                    *outStream << family->name << " " << i << " ";
                    for (int j = 0; j < haps.size(); j++) {
                        *outStream << haps[j] << " ";
                    }
                    *outStream << endl;
                }
            }

            // merge consistentHapsSibship into a single list for this sibship
            vector<int> sibshipHaps;
            // first store the indices of the sibs with any consistent haplotypes
            int nconsistentSib = 0;
            for (int i = 0; i < family->sibs.size(); i++)
                if (consistentHapsSibship[i].size()) {
                    sibshipHaps.push_back(i);
                }
            // signal the end of sib indices
            sibshipHaps.push_back(-1);

            // now list them out
            ft.clear();
            fnt.clear();
            mt.clear();
            mnt.clear();

            vector<bool> MchrX;
            MchrX.clear();
            for (int i = 0; i < family->sibs.size(); i++)
                if (consistentHapsSibship[i].size()) {
                    MchrX.push_back(options.chrX && family->sibs[i].sex == MALE);
                }
            // recursively check for consistency with all sibs
            if (family->sibship) {
                listConsistentAllSibsSibship(ft, mt, consistentHapsSibship, MchrX, sibshipHaps, options, 0);
            } else {
                listConsistentAllSibsFamily(ft, fnt, mt, mnt, consistentHapsSibship, MchrX, sibshipHaps, options, 1, 0);
            }
            // store consistent haplotypes in the master list
            if (options.show) {
                *outStream << family->name << " ";
                for (int j = 0; j < sibshipHaps.size(); j++) {
                    *outStream << sibshipHaps[j] << " ";
                }
                *outStream << endl;
            }

            consistentHaps.push_back(sibshipHaps);

            // note whether we have a homozygous intercross
            homix.push_back(homoIntercross);
        }

    // loop through unrelateds
    for (vector<Subject>::iterator subject = unrelateds.begin();
            subject != unrelateds.end(); subject++) {
        int size = options.window + options.condition.size() + options.tag.size();
        Haplotype hap1(size);
        Haplotype hap2(size);
        vector<int> haps;
        if (subject->usable) {
            listHaplotypesUnrelated(haps, hap1, hap2, *subject, options, 0);
        }
        if (options.show) {
            *outStream << subject->id << " ";
            for (int j = 0; j < haps.size(); j++) {
                *outStream << genoCode[haps[j]].str(options.condgenotype * options.condition.size(), options.genotype, ACGT) << " ";
            }
            *outStream << endl;
        }
        consistentHaps.push_back(haps);
        hap1.clear();
        hap2.clear();
    }
    realisationProb.resize(consistentHaps.size());
}

// list haplotypes in a family
void UnphasedAnalysis::listHaplotypesFamily(vector<int> &haps,
        Haplotype &ft, Haplotype &fnt, Haplotype &mt,
        Haplotype &mnt, NuclearFamily &family, int sib,
        UnphasedOptions &options,
        bool &homoIntercross, int i) {
    // construct the haplotypes
    if (i < ft.size()) {
        int marker = currentmarker[i];

        // transmitted/untransmitted alleles at this marker
        int fa1 = 0, fa2 = 0, ma1 = 0, ma2 = 0;

        // intercross applies to just this marker
        // homoIntercross applies ("OR") across all markers and sibs
        bool intercross = false;

        // determine transmitted alleles from the family data
        determinetransmission(family, options, sib, marker, fa1, fa2, ma1, ma2,
                              intercross, homoIntercross);

        // number of loci that could introduce phase uncertainty
        int phaseLoci = options.condition.size() * (options.condgenotype ? 0 : 1) +
                        (options.window + options.tag.size()) * (options.genotype ? 0 : 1);

        // return if only using certain haplotypes
        if (options.certain && intercross && phaseLoci > 1) {
            return;
        }

        // ignore homoIntercross if estimating missing data
        if (!options.certain) {
            homoIntercross = false;
        }

        // ignore homoIntercross if just one marker
        if (phaseLoci <= 1) {
            homoIntercross = false;
        }

        // loop through the possible alleles when they are missing
        bool MchrX = options.chrX && family.sibs[sib].sex == MALE;
        bool FchrX = options.chrX && family.sibs[sib].sex == FEMALE;
        int nallele = allele[i].size();
        bool done1 = false;
        for (int i1 = 0; !done1 && i1 < nallele; i1++) {
            if (!ma1 && options.missing) {
                mt[i] = allele[i][i1];
            } else {
                mt[i] = ma1;
                done1 = true;
            }
            bool done2 = false;
            for (int i2 = 0; !done2 && i2 < nallele; i2++) {
                if (!ma2 && options.missing) {
                    mnt[i] = allele[i][i2];
                } else {
                    mnt[i] = ma2;
                    done2 = true;
                }
                bool done3 = false;
                for (int i3 = 0; !done3 && i3 < nallele; i3++) {
                    if (!fa1 && options.missing) {
                        ft[i] = allele[i][i3];
                    } else {
                        ft[i] = fa1;
                        done3 = true;
                    }
                    bool done4 = false;
                    for (int i4 = 0; !done4 && i4 < nallele; i4++) {
                        if (!fa2 && options.missing) {
                            // force homozygote father on chrX
                            if (!options.chrX || allele[i][i4] == ft[i]) {
                                fnt[i] = allele[i][i4];
                            }
                        } else {
                            fnt[i] = fa2;
                            done4 = true;
                        }
                        if (possiblegenotype(family.father, marker, ft[i], fnt[i]) &&
                                possiblegenotype(family.mother, marker, mt[i], mnt[i]) &&
                                possiblegenotype(family.sibs[sib], marker, MchrX ? mt[i] : ft[i], mt[i])) {
                            // recursively move to the next marker
                            listHaplotypesFamily(haps, ft, fnt, mt, mnt, family, sib, options,
                                                 homoIntercross, i + 1);
                            // if intercross, do other phase
                            // relevant for genotype tests also (twice the probability...)
                            if (intercross) {
                                swap(ft[i], fnt[i]);
                                swap(mt[i], mnt[i]);
                                listHaplotypesFamily(haps, ft, fnt, mt, mnt, family, sib, options,
                                                     homoIntercross, i + 1);
                            }
                        }
                        fnt[i] = 0;
                    }
                    ft[i] = 0;
                }
                mnt[i] = 0;
            }
            mt[i] = 0;
        }

    }

    // complete haplotype - store the possibilities
    else {
        if (ft.nozero() && fnt.nozero() &&
                mt.nozero() && mnt.nozero()) {
            ft.enter(haploCode);
            fnt.enter(haploCode);
            mt.enter(haploCode);
            mnt.enter(haploCode);
            haps.push_back(ft.index(haploCode));
            haps.push_back(fnt.index(haploCode));
            haps.push_back(mt.index(haploCode));
            haps.push_back(mnt.index(haploCode));
            // store the genotypes for posterity
            Haplotype genotype[8];
            int ncondition = options.condgenotype * options.condition.size();
            genotype[0] = ft.formgeno(mt, ncondition, options.genotype);
            genotype[1] = ft.formgeno(mnt, ncondition, options.genotype);
            genotype[2] = fnt.formgeno(mt, ncondition, options.genotype);
            genotype[3] = fnt.formgeno(mnt, ncondition, options.genotype);
            genotype[4] = mt.formgeno(ft, ncondition, options.genotype);
            genotype[5] = mnt.formgeno(ft, ncondition, options.genotype);
            genotype[6] = mt.formgeno(fnt, ncondition, options.genotype);
            genotype[7] = mnt.formgeno(fnt, ncondition, options.genotype);
            if (options.chrX && family.sibs[sib].sex == MALE) {
                // on chrX, make male sibs homozygous for maternal haplotypes
                genotype[0] = mt.formgeno(mt, ncondition, options.genotype);
                genotype[1] = genotype[0];
                genotype[2] = mnt.formgeno(mnt, ncondition, options.genotype);
                genotype[3] = genotype[2];
                for (int j = 4; j < 8; j++) {
                    genotype[j] = genotype[j-4];
                }
            }
            for (int j = 0; j < 8; j++) {
                genotype[j].enter(genoCode);
            }
        }
    }
}

// find haplotype solutions that are consistent with all sibs
void UnphasedAnalysis::listConsistentAllSibsFamily(Haplotype &ft,
        Haplotype &fnt,
        Haplotype &mt,
        Haplotype &mnt,
        vector<vector<int> > &consistentHapsSibship,
        vector<bool> &MchrX,
        vector<int> &sibshipHaps,
        UnphasedOptions &options,
        int nphase,
        int i) {
    if (i < consistentHapsSibship.size()) {
        // if nothing for this sib, move on
        if (!consistentHapsSibship[i].size()) {
            listConsistentAllSibsFamily(ft, fnt, mt, mnt, consistentHapsSibship, MchrX, sibshipHaps, options, nphase, i + 1);
            return;
        }
        // recursively form combinations of trans/untrans haplotypes
        for (int j = 0; j < consistentHapsSibship[i].size(); j += 4) {
            int ft1 = consistentHapsSibship[i][j];
            int fnt1 = consistentHapsSibship[i][j+1];
            int mt1 = consistentHapsSibship[i][j+2];
            int mnt1 = consistentHapsSibship[i][j+3];
            // only allow one paternal and one maternal genotype
            if (ft.size() == 0 ||
                    min(ft1, fnt1) == min(ft[0], fnt[0]) &&
                    max(ft1, fnt1) == max(ft[0], fnt[0]) &&
                    min(mt1, mnt1) == min(mt[0], mnt[0]) &&
                    max(mt1, mnt1) == max(mt[0], mnt[0])) {
                ft.push_back(ft1);
                fnt.push_back(fnt1);
                mt.push_back(mt1);
                mnt.push_back(mnt1);
                int thisphase = 1;
                // count alternative phases from homozygous parents
                if (ft.size() > 1) {
                    thisphase += (ft1 == fnt1) + (mt1 == mnt1) + (ft1 == fnt1 && mt1 == mnt1);
                }
                listConsistentAllSibsFamily(ft, fnt, mt, mnt, consistentHapsSibship, MchrX, sibshipHaps, options, nphase * thisphase, i + 1);
                ft.pop_back();
                fnt.pop_back();
                mt.pop_back();
                mnt.pop_back();
            }
        }
    } else {
        // save all offspring haplotypes
        int ncondition = options.condgenotype * options.condition.size();
        sibshipHaps.push_back(nphase);
        sibshipHaps.push_back(1);
        sibshipHaps.push_back(1);
        for (int j = 0; j < ft.size(); j++) {
            Haplotype ftj = haploCode[ft[j]];
            Haplotype fntj = haploCode[fnt[j]];
            Haplotype mtj = haploCode[mt[j]];
            Haplotype mntj = haploCode[mnt[j]];
            Haplotype genotype[8];
            genotype[0] = ftj.formgeno(mtj, ncondition, options.genotype);
            genotype[1] = ftj.formgeno(mntj, ncondition, options.genotype);
            genotype[2] = fntj.formgeno(mtj, ncondition, options.genotype);
            genotype[3] = fntj.formgeno(mntj, ncondition, options.genotype);
            genotype[4] = mtj.formgeno(ftj, ncondition, options.genotype);
            genotype[5] = mntj.formgeno(ftj, ncondition, options.genotype);
            genotype[6] = mtj.formgeno(fntj, ncondition, options.genotype);
            genotype[7] = mntj.formgeno(fntj, ncondition, options.genotype);
            if (MchrX[j]) {
                // on chrX, make male sibs homozygous for maternal haplotypes
                genotype[0] = mtj.formgeno(mtj, ncondition, options.genotype);
                genotype[1] = genotype[0];
                genotype[2] = mntj.formgeno(mntj, ncondition, options.genotype);
                genotype[3] = genotype[2];
                for (int k = 4; k < 8; k++) {
                    genotype[k] = genotype[k-4];
                }
            }
            for (int k = 0; k < 8; k++) {
                sibshipHaps.push_back(genotype[k].index(genoCode));
            }
        }
    }
}

// find haplotype solutions that are consistent with all sibs
void UnphasedAnalysis::listConsistentAllSibsSibship(Haplotype &ft,
        Haplotype &mt,
        vector<vector<int> > &consistentHapsSibship,
        vector<bool> &MchrX,
        vector<int> &sibshipHaps,
        UnphasedOptions &options,
        int i) {
    if (i < consistentHapsSibship.size()) {
        // if nothing for this sib, move on
        if (!consistentHapsSibship[i].size()) {
            listConsistentAllSibsSibship(ft, mt, consistentHapsSibship, MchrX, sibshipHaps, options, i + 1);
            return;
        }
        // see if each of the phases of this sib are consistent with the others
        bool consistent = false;
        for (int j = 0; j < consistentHapsSibship[i].size(); j += 2) {
            int ft1 = consistentHapsSibship[i][j];
            int mt1 = consistentHapsSibship[i][j+1];
            ft.push_back(ft1);
            mt.push_back(mt1);
            listConsistentAllSibsSibship(ft, mt, consistentHapsSibship, MchrX, sibshipHaps, options, i + 1);
            ft.pop_back();
            mt.pop_back();
        }
    } else {

        // see if more than 2 paternal or maternal haplotypes
        vector<int> fhaps, mhaps;
        for (int j = 0; j < ft.size(); j++) {
            if (fhaps.size() == 2 && ft[j] != fhaps[0] && ft[j] != fhaps[1]) {
                return;
            }
            if (mhaps.size() == 2 && mt[j] != mhaps[0] && mt[j] != mhaps[1]) {
                return;
            }
            if (fhaps.size() == 0 || fhaps.size() == 1 && ft[j] != fhaps[0]) {
                fhaps.push_back(ft[j]);
            }
            if (mhaps.size() == 0 || mhaps.size() == 1 && mt[j] != mhaps[0]) {
                mhaps.push_back(mt[j]);
            }
        }

        // just one phase, as we can't deduce homozygous parents
        sibshipHaps.push_back(1);
        bool knowFather = (fhaps.size() == 2);
        bool knowMother = (mhaps.size() == 2);
        sibshipHaps.push_back(knowFather);
        sibshipHaps.push_back(knowMother);

        if (fhaps.size() == 1) {
            fhaps.push_back(fhaps[0]);
        }
        if (mhaps.size() == 1) {
            mhaps.push_back(mhaps[0]);
        }

        // save all offspring haplotypes
        int ncondition = options.condgenotype * options.condition.size();
        for (int j = 0; j < ft.size(); j++) {
            Haplotype ftj = haploCode[ft[j]];
            Haplotype fntj = haploCode[fhaps[ft[j] == fhaps[0]]];
            Haplotype mtj = haploCode[mt[j]];
            Haplotype mntj = haploCode[mhaps[mt[j] == mhaps[0]]];

            Haplotype genotype[8];
            genotype[0] = ftj.formgeno(mtj, ncondition, options.genotype);
            genotype[1] = ftj.formgeno(mntj, ncondition, options.genotype);
            genotype[2] = fntj.formgeno(mtj, ncondition, options.genotype);
            genotype[3] = fntj.formgeno(mntj, ncondition, options.genotype);
            genotype[4] = mtj.formgeno(ftj, ncondition, options.genotype);
            genotype[5] = mntj.formgeno(ftj, ncondition, options.genotype);
            genotype[6] = mtj.formgeno(fntj, ncondition, options.genotype);
            genotype[7] = mntj.formgeno(fntj, ncondition, options.genotype);
            if (MchrX[j]) {
                // on chrX, make male sibs homozygous for maternal haplotypes
                genotype[0] = mtj.formgeno(mtj, ncondition, options.genotype);
                genotype[1] = genotype[0];
                genotype[2] = mntj.formgeno(mntj, ncondition, options.genotype);
                genotype[3] = genotype[2];
                for (int k = 4; k < 8; k++) {
                    genotype[k] = genotype[k-4];
                }
            }
            for (int k = 0; k < 8; k++) {
                sibshipHaps.push_back(genotype[k].index(genoCode));
            }
        }

    }
}


// list haplotypes for single unrelated
void UnphasedAnalysis::listHaplotypesUnrelated(vector<int> &haps,
        Haplotype &hap1, Haplotype &hap2,
        Subject &subject, UnphasedOptions &options, int i) {
    // construct the haplotypes
    if (i < hap1.size()) {
        int marker = currentmarker[i];
        // transmitted/untransmitted alleles at this marker
        //     int a1=subject.marker[marker][0];
        //     int a2=subject.marker[marker][1];
        int a1 = subject.marker[marker] / 16;
        int a2 = subject.marker[marker] % 16;

        // force a homozygote for male chrX
        bool chrX = options.chrX && subject.sex == MALE || options.chrY;
        if (chrX) {
            a1 = max(a1, a2);
            a2 = a1;
        }

        // number of loci that could introduce phase uncertainty
        int phaseLoci = options.condition.size() * (options.condgenotype ? 0 : 1) +
                        (options.window + options.tag.size()) * (options.genotype ? 0 : 1);

        // for certain haplotypes only,
        // reject if heterozgyous and looking at more than one phased locus
        if (options.certain && !chrX &&
                (a1 != a2 && phaseLoci > 1 &&
                 !(options.condgenotype && i < options.condition.size()) &&
                 !options.genotype || !a1 || !a2)) {
            return;
        }

        // loop through the possible alleles when they are missing
        int nallele = allele[i].size();
        bool done1 = false;
        for (int i1 = 0; !done1 && i1 < nallele; i1++) {
            if (!a1 && options.missing) {
                hap1[i] = allele[i][i1];
            } else {
                hap1[i] = a1;
                done1 = true;
            }
            bool done2 = false;
            for (int i2 = 0; !done2 && i2 < nallele; i2++) if (!chrX || i2 == i1) {
                    if (!a2 && options.missing) {
                        hap2[i] = allele[i][i2];
                    } else {
                        hap2[i] = a2;
                        done2 = true;
                    }
                    if (hap1[i] == 0 || hap2[i] == 0) {
                        return;
                    }
                    if (possiblegenotype(subject, marker, hap1[i], hap2[i])) {
                        // recursively move to the next marker
                        listHaplotypesUnrelated(haps, hap1, hap2, subject, options, i + 1);
                        // if heterozygous, do other phase
                        if (hap1[i] != hap2[i] &&
                                // other phase of missing genotype will come separately
                                (a1 || a2)) {
                            swap(hap1[i], hap2[i]);
                            listHaplotypesUnrelated(haps, hap1, hap2, subject, options, i + 1);
                            swap(hap1[i], hap2[i]);
                        }
                    }
                    hap2[i] = 0;
                }
            hap1[i] = 0;
        }
    }

    // complete haplotype - store the possibilities
    else {
        if (hap1.nozero() && hap2.nozero()) {
            // make second haplotype zero for male chrX
            int ncondition = options.condgenotype * options.condition.size();
            Haplotype geno1 = hap1.formgeno(hap2, ncondition, options.genotype);
            Haplotype geno2 = hap2.formgeno(hap1, ncondition, options.genotype);
            geno1.enter(genoCode);
            geno2.enter(genoCode);
            hap1.enter(haploCode);
            hap2.enter(haploCode);
            haps.push_back(geno1.index(genoCode));
            haps.push_back(geno2.index(genoCode));
        }
    }
}



//  exploratory

void UnphasedAnalysis::exploratory(UnphasedOptions &options, const string &which) {
    // make exploratory pass through the data
    // storing a list of compatible haplotypes in each family and subject
    if (which == "asymptotic") {
        *outStream << "Exploratory pass..." << flush;
    }
    findUsableSubjects(options);
    listConsistent(options);

    if (which == "asymptotic") {
        *outStream << "done" << endl;
    }

    if (options.show) {
        for (int i = 0; i < genoCode.size(); i++) {
            Haplotype hap = genoCode[i];
            *outStream << hap.str(options.condition.size()*options.condgenotype, options.genotype, ACGT) << " " << i << endl;
        }
        for (int i = 0; i < haploCode.size(); i++) {
            Haplotype hap = haploCode[i];
            *outStream << hap.str(ACGT) << " " << i << endl;
        }
    }

    resizeArrays();
    if (options.model == "allelemain" || options.model == "commonmain") {
        group.resize((options.condition.size() > 0) + options.window);
    } else {
        if (options.model == "haplomain" || options.model == "gxg") {
            group.resize(2);
        } else {
            group.resize(1);
        }
    }
    for (int i = 0; i < group.size(); i++) {
        group[i].resize(genoCode.size(), 0);
    }

    // if a reference haplotype is specified, force it into the list of genotypes
    if (options.reference.size()) {
        reference = options.reference;
        reference.enter(genoCode);
    } else if (genoCode.size()) {
        reference = genoCode[0];
    }

    // estimate frequencies and set up "rare" and "zero"
    // estimate under the full alternative model, the one with max likelihood
    if (options.rare || options.zero) {
        int df = 0;
        if (which == "asymptotic") {
            *outStream << "Identifying rare haplotypes: " << flush;
        }
        fullAlternative(options, which, df);
        double totalCount[2];
        for (int i = 0; i < 2; i++) {
            totalCount[i] = 0;
            for (int j = 0; j < genoCode.size(); j++) {
                totalCount[i] += familyCount[i][j] + unrelatedCount[i][j];
            }
        }
        for (int i = 0; i < genoCode.size(); i++) {
            double x0 = familyCount[0][i] + unrelatedCount[0][i];
            double x1 = familyCount[1][i] + unrelatedCount[1][i];
            if (!options.cellcount) {
                x0 /= totalCount[0];
                x1 /= totalCount[1];
            }
            if (options.userare == "both" && typeOfPhenotype != "quant") {
                rare[i] = x0 < options.rare && x1 < options.rare;
            }
            if (options.userare == "either" && typeOfPhenotype != "quant") {
                rare[i] = x0 < options.rare || x1 < options.rare;
            }
            if (options.userare == "case" && typeOfPhenotype != "quant") {
                rare[i] = x1 < options.rare;
            }
            if (options.userare == "control" || typeOfPhenotype == "quant") {
                rare[i] = x0 < options.rare;
            }

            // zero haplotypes must be zero in both
            zero[i] = ((familyCount[0][i] + unrelatedCount[0][i]) / totalCount[0] < options.zero) &&
                      (typeOfPhenotype == "quant" || (familyCount[1][i] + unrelatedCount[1][i]) / totalCount[1] < options.zero);

            // if reference haplotype specified, don't let it be zero
            zero[i] = zero[i] && genoCode[i] != options.reference;
        }
    }

    // if no reference haplotype specified, make it the first non-zero one
    if (!options.reference.size()) {
        reference.clear();
        vector<Haplotype> sortedHaps;
        sortedHaps = genoCode;
        sort(sortedHaps.begin(), sortedHaps.end());
        int conditionsize = options.condition.size() * (options.condgenotype ? 2 : 1);
        int testsize = options.window * (options.genotype ? 2 : 1);
        Haplotype specific;
        specific = options.specific;
        for (int i = 0; i < sortedHaps.size() && !reference.size(); i++) {
            int j = sortedHaps[i].index(genoCode);
            Haplotype testhap;
            testhap.assign(sortedHaps[i].begin() + conditionsize, sortedHaps[i].begin() + conditionsize + testsize);
            if (!zero[j] && !(options.specific.size() && weq(testhap, specific))) {
                reference = sortedHaps[i];
            }
        }
    }
}



