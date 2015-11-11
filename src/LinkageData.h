/* LinkageData.h - manage Linkage format files

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

#include "stl.h"
#include <fstream>
#include <string>

using namespace std;

#ifndef __LINKAGEDATA__

// Number of levels of polytomous outcome
const int K = 4;
const int AFFECTED = 2;
const int UNAFFECTED = 1;
const int UNSURE = 0;
const int MALE = 1;
const int FEMALE = 2;
const string MISSING = "-";

class Subject {
public:
    string id;
    string pid;
    string mid;
    int sex;
    vector<int> affection;
    //  vector<vector<int> > marker;
    vector<unsigned char> marker; // more efficient for GWAS data
    vector<double> trait;
    vector<bool> missingtrait;
    bool usable;
};

class Locus {
public:
    string name;
    int type;
    int nallele;
    int chromosome;
    int position;
    vector<double> frequency;
    int nliability;
    vector<vector<double> > penetrance;
    vector<vector<double> > malepenetrance;
    vector<vector<double> > variance;
    float hetmultiplier;
    vector<int> allele;
    void readLocus(ifstream &, ostream * &, int);
};

typedef map<string, vector<Subject> > PedigreeSet;

class LinkageData {
public:
    LinkageData();
    LinkageData(ostream *);
    bool ACGT;
protected:
    void readpedfile(string &);
    void readpedfile(string &, string &);
    void readdatafile(string &);
    void readphenofile(string &);
    void defaultdatafile(string &);
    void readmapfile(string &);
    void writedatafile();
    void makehashtables();
    PedigreeSet pedigree;
    int nlocus;
    int risklocus;
    int sexlinked;
    int program;
    int mutlocus;
    double malemut;
    double femalemut;
    int hapfreq;
    vector<Locus> locus;
    vector<int> locuscount;
    map<string, int> diseasehash;
    map<string, int> markerhash;
    map<string, int> traithash;
    vector<int> order;
    int sexdiff;
    int interference;
    vector<double> recombination;
    ostream *outStream;
};

#endif
#define __LINKAGEDATA__
