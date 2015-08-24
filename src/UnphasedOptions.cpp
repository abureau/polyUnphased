/* UnphasedOptions.cpp - options for UNPHASED

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


#include "UnphasedOptions.h"

using namespace std;

// constructor
UnphasedOptions::UnphasedOptions(int argc, char **argv): Options(argc, argv, false) {
    option("-optionfile", optionFilename, "");
    if (optionFilename != "") {
        readOptionFile();
    }
    option("-pedfile", inputfile);
    option("-marker", marker);
    option("-condition", condition);
    option("-tag", tag);
    option("-trait", trait);
    option("-disease", disease);
    option("-joint", joint);
    option("-confounder", confounder);
    option("-modifier", modifier);
    option("-factor", factor);
    option("-baseline", baseline);
    unphasedOption("-condspecific", condspecific);
    unphasedOption("-reference", reference);
    unphasedOption("-compare", compare1);
    unphasedOption("-with", compare2);
    unphasedOption("-specific", specific);
    option("-datafile", datafile, "");
    option("-mapfile", mapfile, "");
    option("-bedfile", bedfile, "");
    option("-phenofile", phenofile, "");
    option("-output", outputFilename, "");
    option("-dumpfile", dumpFilename, "");
    option("-tabularfile", tabularFilename, "");
    option("-listmarkerfile", listMarkerFilename, "");
    option("-model", model, "full");
    option("-userare", userare, "both");
    option("-permutation", npermutation, 0);
    option("-quantile", quantile, 5);
    option("-window", window, 1);
    option("-restarts", restarts, 0);
    option("-randomseed", randomseed, 1);
    option("-rare", rare, 0);
    option("-zero", zero, 1e-8);
    option("-epsilon", epsilon, 1e-8);
    option("-variance", variance, 1.0);
    option("-covariance", covariance, 0.0);
    option("-certain", certain, false);
    option("-missing", missing, false);
    option("-chrX", chrX, false);
    option("-chrY", chrY, false);
    option("-LD", LD, false);
    option("-neldermead", neldermead, false);
    option("-slow", slow, false);
    option("-cellcount", cellcount, false);
    option("-permoutput", outputPermutation, false);
    option("-brief", briefOutput, false);
    option("-mostlikely", mostlikely, false);
    option("-condgenotype", condgenotype, false);
    option("-genotype", genotype, false);
    option("-individual", individual, false);
    option("-testconfounders", testconfounders, false);
    option("-testmodifiers", testmodifiers, false);
    option("-nolinkage", nolinkage, false);
    option("-sibship", sibship, false);
    option("-llhd", llhd, false);
    option("-show", show, false);
    option("-parentrisk", parentrisk, false);
    option("-hhrr", hhrr, false);
    option("-normal", normal, false);
    option("-uncentred", uncentred, false);
    option("-follow", follow, false);
    option("-time", time, false);
    option("-allwindows", allwindows, false);
    option("-allcombinations", allcombinations, false);
    option("-onefbc", onefbc, false);
    option("-imputeDump", imputeDump, false);
    option("-checkgradient", checkgradient, false);
    parse();
}

// destructor
UnphasedOptions::~UnphasedOptions() {}

void UnphasedOptions::parse() {
    Options::parse();
    if (missing) {
        certain = false;
    }
    if (certain) {
        restarts = 0;
    }
    if (rare >= 1.0 && !cellcount) {
        rare /= 100.0;
    }
    if (genotype) {
        condgenotype = true;
    }
    if (hhrr) {
        parentrisk = false;
    }
    if (sibship) {
        nolinkage = false;
        if (genotype || condgenotype) {
            cout << "WARNING: -genotype not yet implemented with -sibship" << endl;
            genotype = false;
        }
        //    parentrisk=true;
    }
    if (model != "" &&
            model != "full" &&
            model != "haplomain" &&
            model != "allelemain" &&
            model != "commonmain" &&
            model != "gxg" &&
            model != "null") {
        cout << "ERROR: unrecognised argument to -model: " << model << endl;
        exit(-1);
    }
}


void UnphasedOptions::readOptionFile() {
    ifstream infile(optionFilename.c_str());
    if (!infile) {
        cout << "Failed to open option file " << optionFilename << endl;
        exit(-1);
    }
    cout << "Reading option file " << optionFilename << "..." << flush;
    string line;
    while ((line = getline(infile)) != "") {
        istrstream instr(line.c_str());
        bool virgin = true;
        while (!instr.eof()) {
            string arg;
            instr >> arg;
            if (virgin && arg[0] != '-') {
                arg = "-" + arg;
            }
            argc++;
            argv.push_back(arg);
            virgin = false;
        }
    }
    cout << "done" << endl;
}

// several int arguments which could be ACGT on the command line
void UnphasedOptions::unphasedOption(string s, vector<int> &arg) {
    OptionInfo o(-1);
    arg.clear();
    for (int i = 1; i < argc; i++) {
        if (s.find(argv[i]) == 0) {
            for (i++; i < argc && argv[i][0] != '-'; i++) {
                int val;
                const string letters = "NnAaCcGgTt";
                int index = letters.find(argv[i][0]) / 2;
                if (index >= 0 && index < 5 && argv[i][1]==0) {
                    val=index;
                }
                else {    
                    istrstream istr(argv[i].c_str());
                    istr >> val;
                }
                arg.push_back(val);
            }
            o.nargs = arg.size();
        }
    }
    optionList.insert(make_pair(s, o));
}

