/* NuclearFamily.cpp - managing nuclear families

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
#include <strstream>
#include <stdio.h>
#include "NuclearFamily.h"

NuclearFamily::NuclearFamily() {}

// get family members and check family structure
NuclearFamily::NuclearFamily(string s, vector<Subject> &subject) {
    name = s;
    // identify the parents, may be missing
    hasfather = false;
    hasmother = false;
    status = "";
    for (vector<Subject>::iterator i = subject.begin();
            i != subject.end(); i++) {
        bool isfather = false;
        bool ismother = false;
        vector<Subject>::iterator j = subject.begin();
        while (j != subject.end() && (j->pid == "0" || j->pid != i->id)) {
            j++;
        }
        if (j != subject.end()) {
            if (i->sex == FEMALE) {
                status = "parental sex mismatch";
            }
            if (i == j) {
                status = "father ID same as subject ID";
            }
            if (hasfather) {
                status = "two fathers";
            }
            if (status == "") {
                father = *i;
                hasfather = true;
                isfather = true;
            }
        }
        j = subject.begin();
        while (j != subject.end() && (j->mid == "0" || j->mid != i->id)) {
            j++;
        }
        if (j != subject.end()) {
            if (i->sex == MALE) {
                status = "parental sex mismatch";
            }
            if (i == j) {
                status = "mother ID same as subject ID";
            }
            if (hasmother) {
                status = "two mothers";
            }
            if (status == "") {
                mother = *i;
                hasmother = true;
                ismother = true;
            }
        }
        if (!isfather && !ismother) {
            sibs.push_back(*i);
        }
    }
    if (sibs.size() == 0) {
        status = "no offspring";
    } else {
        // create dummy parents if required
        //     int maxid=0;
        //     for(int i=0;i<sibs.size();i++)
        //       if (sibs[i].id>maxid) maxid=sibs[i].id;
        //    vector<int> nullgenotype(2,0);
        char nullgenotype = 0;
        if (!hasfather) {
            father.id = "__auto_father__"; //maxid+1;
            father.pid = "0";
            father.mid = "0";
            father.sex = MALE;
            father.affection.clear();
            for (int i = 0; i < sibs[0].affection.size(); i++) {
                father.affection.push_back(UNSURE);
            }
            father.marker.clear();
            for (int i = 0; i < sibs[0].marker.size(); i++) {
                father.marker.push_back(nullgenotype);
            }
            father.trait.clear();
            for (int i = 0; i < sibs[0].trait.size(); i++) {
                father.trait.push_back(0);
                father.missingtrait.push_back(true);
            }
        }
        if (!hasmother) {
            mother.id = "__auto_mother__"; //maxid+2;
            mother.pid = "0";
            mother.mid = "0";
            mother.sex = FEMALE;
            mother.affection.clear();
            for (int i = 0; i < sibs[0].affection.size(); i++) {
                mother.affection.push_back(UNSURE);
            }
            mother.marker.clear();
            for (int i = 0; i < sibs[0].marker.size(); i++) {
                mother.marker.push_back(nullgenotype);
            }
            mother.trait.clear();
            for (int i = 0; i < sibs[0].trait.size(); i++) {
                mother.trait.push_back(0);
                mother.missingtrait.push_back(true);
            }
        }
        sibship = (!hasfather && !hasmother);
    }
}

// report errors
void NuclearFamily::report() {
    if (status != "") {
        cout << "Warning: " << status << " in pedigree "
             << name << ": skipped" << endl;
        for (vector<Subject>::iterator i = sibs.begin(); i < sibs.end(); i++) {
            printf("%s %s %s %s %d\n", name.c_str(), i->id.c_str(), i->pid.c_str(), i->mid.c_str(), i->sex);
        }
    }
}

NuclearFamily::~NuclearFamily() {}

FamilyList::FamilyList() {}

inline bool seq(string x, string y) {
    return(x != "0" && y != "0" && x == y);
}

// extract nuclear families from a set of pedigrees
void FamilyList::extract(PedigreeSet &pedigree) {
    for (map<string, vector<Subject> >::iterator i = pedigree.begin();
            i != pedigree.end(); i++) {
        // loop through each subject in turn
        for (vector<Subject>::iterator j = i->second.begin();
                j != i->second.end(); j++) {
            vector<Subject> subjects;
            subjects.clear();
            // check whether this subject is the first with these parents
            bool firstSib = true;
            for (vector<Subject>::iterator k = i->second.begin(); k < j; k++)
                if (seq(k->pid, j->pid) && seq(k->mid, j->mid)) {
                    firstSib = false;
                }
            // if true then extract all sibs with the same parents
            // and the parents themselves
            if (firstSib) {
                for (vector<Subject>::iterator k = i->second.begin();
                        k != i->second.end(); k++) {
                    if (seq(k->pid, j->pid) && seq(k->mid, j->mid) ||
                            seq(k->id, j->pid) || seq(k->id, j->mid) || k == j) {
                        subjects.push_back(*k);
                    }
                }
            }
            // add this to the list of nuclear families
            if (subjects.size() > 1 ||
                    subjects.size() == 1 && (subjects[0].pid != "0" || subjects[0].mid != "0")) {
                NuclearFamily family(i->first, subjects);
                push_back(family);
            }
        }
    }
}

FamilyList::~FamilyList() {}
