/* UnphasedOptions.h - options for UNPHASED

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


#include <fstream>
#include "Options.h"
#include "getline.h"

#ifndef __UNPHASEDOPTIONS__

class UnphasedOptions: public Options {
public:
    vector<string> marker;
    vector<string> condition;
    vector<string> tag;
    vector<string> trait;
    vector<string> disease;
    vector<string> polytomous;
    vector<string> joint;
    vector<string> confounder;
    vector<string> modifier;
    vector<string> factor;
    vector<int> baseline;
    vector<int> condspecific;
    vector<int> reference;
    vector<int> compare1;
    vector<int> compare2;
    vector<int> specific;
    string datafile;
    string mapfile;
    string bedfile;
    string phenofile;
    string outputFilename;
    string dumpFilename;
    string tabularFilename;
    string optionFilename;
    string listMarkerFilename;
    string model;
    string userare;
    int npermutation;
    int quantile;
    int window;
    int restarts;
    int randomseed;
    double rare;
    double zero;
    double epsilon;
    double variance;
    double covariance;
    bool certain;
    bool missing;
    bool chrX;
    bool chrY;
    bool LD;
    bool neldermead;
    bool slow;
    bool cellcount;
    bool outputPermutation;
    bool briefOutput;
    bool mostlikely;
    bool condgenotype;
    bool genotype;
    bool individual;
    bool testconfounders;
    bool testmodifiers;
    bool nolinkage;
    bool sibship;
    bool llhd;
    bool show;
    bool parentrisk;
    bool hhrr;
    bool normal;
    bool uncentred;
    bool follow;
    bool time;
    bool allwindows;
    bool allcombinations;
    bool onefbc;
    bool imputeDump;
    bool checkgradient;

    // constructor
    UnphasedOptions(int, char **);
    UnphasedOptions(vector<string>);

    // destructor
    ~UnphasedOptions();

protected:
    void parse();
    void readOptionFile();
void unphasedOption(string,vector<int> &);
};

#endif
#define __UNPHASEDOPTIONS__
