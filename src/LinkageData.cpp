/* LinkageData.cpp - managing Linkage format files

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
#include <fstream>
#include <strstream>
#include <cstdlib>
#include <string>
#include <cstring>
#include "stl.h"
#include "LinkageData.h"
#include "getline.h"
#include "UnphasedAnalysis.h"

LinkageData::LinkageData(ostream *os = &cout) {
    outStream = os;
}

//  readpedfile */

// read a LINKAGE format pedfile

void LinkageData::readpedfile(string &filename) {
    ifstream infile(filename.c_str());

    if (!infile) {
        *outStream << "Failed to open pedigree file " << filename << endl;
        exit(-1);
    }
    *outStream << "Reading pedigree file " << filename << "..." << flush;

	Kvec.resize(nlocus, 2);
    string line;
    while ((line = getline(infile)) != "") {
        Subject subject;
        istrstream instr(line.c_str());
        char buf[256];
        buf[0] = 0;
        // pedigree name
        if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
            *outStream << "\nPedigree file error: premature end of line: expecting pedigree name\n";
            exit(-1);
        }
        string name = buf;

        // subject id
        if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
            *outStream << "\nPedigree file error: premature end of line: expecting subject ID\n";
            exit(-1);
        }
        subject.id = buf;

        // paternal id
        if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
            *outStream << "\nPedigree file error: premature end of line: expecting paternal ID\n";
            exit(-1);
        }
        subject.pid = buf;

        // maternal id
        if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
            *outStream << "\nPedigree file error: premature end of line: expecting maternal ID\n";
            exit(-1);
        }
        subject.mid = buf;

        // sex
        if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
            *outStream << "\nPedigree file error: premature end of line: expecting sex\n";
            exit(-1);
        }
        subject.sex = atoi(buf);

        // locus info
        for (int i = 0; i < nlocus; i++) {
            switch (locus[i].type) {
            case 4:
            case 0: {
                // quantitative trait or covariate
                if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
                    *outStream << "\nPedigree file error: premature end of line: expecting value for trait: " << locus[i].name.c_str() << endl;
                    exit(-1);
                }
                string value = buf;
                if (value == MISSING) {
                    subject.trait.push_back(0);
                    subject.missingtrait.push_back(true);
                } else {
                    subject.trait.push_back(atof(value.c_str()));
                    subject.missingtrait.push_back(false);
                }
                break;
            }
            case 1: {
                // affection status
                if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
                    *outStream << "\nPedigree file error: premature end of line: expecting status for affection: " << locus[i].name.c_str() << endl;
                    exit(-1);
                }
                //subject.affection.insert(make_pair(locus[i].name,atoi(buf)));
                // Determine number of levels of the phenotype (assuming all levels are represented)
                if ((int) atoi(buf) > 4) {
                    *outStream << "\nPedigree file error: value > 4 for an affection locus: " << locus[i].name.c_str() << endl;
                    exit(-1);
                }
                if ((int) atoi(buf) > 4) {
                    *outStream << "\nPedigree file error: negative value for an affection locus: " << locus[i].name.c_str() << endl;
                    exit(-1);
                }
                // Déterminer la valeur maximale du phénotype dans Kvec
                if ((int) atoi(buf) > Kvec[i]) Kvec[i] = (int) atoi(buf);
                subject.affection.push_back(atoi(buf));
                break;
            }
            case 3:
                // marker alleles
                vector<int> genotype(2, 0);
                for (int j = 0; j < 2; j++) {
                    if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
                        *outStream << "\nPedigree file error: premature end of line: expecting alleles for marker: " << locus[i].name.c_str() << endl;
                        exit(-1);
                    }
                    const string letters = "NnAaCcGgTt";
                    int index = letters.find(buf[0]) / 2;
                    if (index >= 0 && index < 5) {
                        ACGT = true;
                        genotype[j] = index;
                    } else {
                        genotype[j] = atoi(buf);
                        if (genotype[j] != 0) {
                            ACGT = false;
                        }
                    }
                }
                //if (genotype[0]==0 || genotype[1]>0 && genotype[0]>genotype[1])
                //swap(genotype[0],genotype[1]);
                //subject.marker.insert(make_pair(locus[i].name,genotype));

                // previous version for vector<int>
                //  subject.marker.push_back(genotype);

                // new version for char
                subject.marker.push_back((unsigned char)(genotype[0] * 16 + genotype[1]));
                break;
            } // switch (locus[i].type)
        } // for locus

        subject.usable = true;
        pedigree[name].push_back(subject);

    } // while line
//    *outStream << "done" << endl;
    *outStream << "done " << Kvec[0] << " " << Kvec[1] << endl;
}



//  readpedfile binary */

// read a PLINK binary format pedfile

void LinkageData::readpedfile(string &filename, string &bedfilename) {
    ifstream infile(filename.c_str());
    if (!infile) {
        *outStream << "Failed to open pedigree file " << filename << endl;
        exit(-1);
    }
    ifstream bedfile(bedfilename.c_str());
    if (!infile) {
        *outStream << "Failed to open binary pedigree file " << filename << endl;
        exit(-1);
    }

    *outStream << "Reading binary pedigree files " << filename << ", " << bedfilename << "..." << flush;

    int magic1 = bedfile.get();
    int magic2 = bedfile.get();
    if (magic1 != 108 || magic2 != 27) {
        *outStream << endl << "ERROR: " << filename << " not reconised as a binary pedigree file" << endl;
        exit(-1);
    }

    // list of pedigree names and ids in the order of the pedfile
    vector<string> nameList;
    vector<int> idList;

    string line;
    while ((line = getline(infile)) != "") {
        Subject subject;
        istrstream instr(line.c_str());
        char buf[256];
        buf[0] = 0;
        // pedigree name
        if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
            *outStream << "\nPedigree file error: premature end of line: expecting pedigree name\n";
            exit(-1);
        }
        string name = buf;

        // subject id
        if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
            *outStream << "\nPedigree file error: premature end of line: expecting subject ID\n";
            exit(-1);
        }
        subject.id = buf;

        // paternal id
        if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
            *outStream << "\nPedigree file error: premature end of line: expecting paternal ID\n";
            exit(-1);
        }
        subject.pid = buf;

        // maternal id
        if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
            *outStream << "\nPedigree file error: premature end of line: expecting maternal ID\n";
            exit(-1);
        }
        subject.mid = buf;

        // sex
        if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
            *outStream << "\nPedigree file error: premature end of line: expecting sex\n";
            exit(-1);
        }
        subject.sex = atoi(buf);

        // locus info
        for (int i = 0; i < 1; i++) {
            switch (locus[i].type) {
            case 4:
            case 0: {
                // quantitative trait or covariate
                if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
                    *outStream << "\nPedigree file error: premature end of line: expecting value for trait: " << locus[i].name.c_str() << endl;
                    exit(-1);
                }
                string value = buf;
                if (value == MISSING) {
                    subject.trait.push_back(0);
                    subject.missingtrait.push_back(true);
                } else {
                    subject.trait.push_back(atof(value.c_str()));
                    subject.missingtrait.push_back(false);
                }
                break;
            }
            case 1: {
                // affection status
                if ((instr >> buf, instr.fail()) || buf[0] == '<' && buf[1] == '<') {
                    *outStream << "\nPedigree file error: premature end of line: expecting status for affection: " << locus[i].name.c_str() << endl;
                    exit(-1);
                }
                //subject.affection.insert(make_pair(locus[i].name,atoi(buf)));
                subject.affection.push_back(atoi(buf));
                break;
            }
            }
        }
        subject.usable = true;
        nameList.push_back(name);
        idList.push_back(pedigree[name].size());
        pedigree[name].push_back(subject);
    } // while line

    // detect SNP-major or individual-major mode
    if (bedfile.get()) {
        for (int i = 1; i < nlocus; i++) {
            for (int subject = 0; subject < nameList.size(); subject += 4) {
                int byte = bedfile.get();
                for (int j = 0; subject + j < nameList.size() && j < 4; j++) {
                    // marker alleles
                    int allele1 = locus[i].allele[0];
                    int allele2 = locus[i].allele[1];
                    int genotype;
                    switch (byte >> (j * 2) & 3) {
                    case 0: {
                        genotype = allele1 * 16 + allele1;    // 1/1
                        break;
                    }
                    case 1: {
                        genotype = 0;    // 0/0
                        break;
                    }
                    case 2: {
                        genotype = allele1 * 16 + allele2;    // 1/2
                        break;
                    }
                    case 3: {
                        genotype = allele2 * 16 + allele2;    // 2/2
                        break;
                    }
                    }
                    pedigree[nameList[subject+j]][idList[subject+j]].marker.push_back((char)genotype);
                } // for j
            } // for subject
        } // for locus
    } // SNP-major
    else {
        for (int subject = 0; subject < nameList.size(); subject++) {
            for (int i = 1; i < nlocus; i += 4) {
                int byte = bedfile.get();
                for (int j = 0; i + j < nlocus && j < 4; j++) {
                    // marker alleles
                    int allele1 = locus[i+j].allele[0];
                    int allele2 = locus[i+j].allele[1];
                    int genotype;
                    switch (byte >> (j * 2) & 3) {
                    case 0: {
                        genotype = allele1 * 16 + allele1;    // 1/1
                        break;
                    }
                    case 1: {
                        genotype = 0;    // 0/0
                        break;
                    }
                    case 2: {
                        genotype = allele1 * 16 + allele2;    // 1/2
                        break;
                    }
                    case 3: {
                        genotype = allele2 * 16 + allele2;    // 2/2
                        break;
                    }
                    }
                    pedigree[nameList[subject]][idList[subject]].marker.push_back((char)genotype);
                } // for j
            } // for locus
        } // for subject

    } // ind-major

    *outStream << "done" << endl;
}



//  readdatafile */

// read a LINKAGE or QTDT format data file

void LinkageData::readdatafile(string &filename) {
    if (filename == "") {
        return;
    }

    ifstream infile(filename.c_str());
    if (!infile) {
        *outStream << "Failed to open data file " << filename << endl;
        exit(-1);
    }
    *outStream << "Reading data file " << filename << flush;

    // first string in data file
    char firstChar = infile.peek();

    // QTDT format
    if (firstChar == 'M' || firstChar == 'T' || firstChar == 'C' || firstChar == 'A' ||
            firstChar == 'm' || firstChar == 't' || firstChar == 'c' || firstChar == 'a') {
        *outStream << " (QTDT format)..." << flush;
        nlocus = 0;
        locuscount.resize(4, 0);
        string line;
        while ((line = getline(infile)) != "") {
            istrstream instr(line.c_str());
            char buf[256];
            buf[0] = 0;
            // read locus type
            if (instr >> buf, instr.fail()) {
                *outStream << "\nData file error: premature end of line: expecting column code\n";
                exit(-1);
            }
            string code = buf;
            if (code != "M" && code != "T" && code != "C" && code != "A" &&
                    code != "m" && code != "t" && code != "c" && code != "a") {
                *outStream << endl << "ERROR: unrecognised code in data file: " << code << endl;
                exit(-1);
            }

            Locus newLocus;
            newLocus.chromosome = 0;
            newLocus.position = nlocus;
            // assign locus info
            if (code == "M" || code == "m") {
                newLocus.type = 3;
                newLocus.nallele = 2;
                newLocus.frequency.resize(2, 0.5);
            }
            if (code == "T" || code == "t") {
                newLocus.type = 0;
                newLocus.nallele = 2;
                newLocus.frequency.resize(2, 0.5);
                newLocus.nliability = 1;
                newLocus.penetrance.resize(1);
                newLocus.penetrance[0].resize(3, 0);
                newLocus.variance.resize(1);
                newLocus.variance[0].resize(1, 1);
                newLocus.hetmultiplier = 1;
            }
            if (code == "C" || code == "c") {
                newLocus.type = 4;
            }
            if (code == "A" || code == "a") {
                newLocus.type = 1;
                newLocus.nallele = 2;
                newLocus.frequency.resize(2, 0.5);
                newLocus.nliability = 1;
                newLocus.penetrance.resize(1);
                newLocus.penetrance[0].resize(3);
                newLocus.penetrance[0][0] = 0; // dominant model
                newLocus.penetrance[0][1] = 1;
                newLocus.penetrance[0][2] = 1;
            }

            // read locus name
            if (instr >> buf, instr.fail()) {
                ostrstream ostr;
                ostr << ++locuscount[newLocus.type%4] << ends;
                newLocus.name = ostr.str();
            } else {
                newLocus.name = buf;
            }
            locus.push_back(newLocus);
            nlocus++;
        }
    }

    // linkage format
    else {

        *outStream << " (Linkage format)..." << flush;

        // first line of data file
        infile >> nlocus >> risklocus >> sexlinked >> program;
        while (infile.peek() != EOF && infile.get() != '\n') {}

        // second line of data file
        infile >> mutlocus >> malemut >> femalemut >> hapfreq;
        while (infile.peek() != EOF && infile.get() != '\n') {}

        // third line: order of loci, skip this line as not sure how GH uses it
        //    order.resize(nlocus);
        //    for(int i=0;i<nlocus;i++) {
        //      infile >> order[i];
        //    }
        order.resize(1);
        infile >> order[0];

        while (infile.peek() != EOF && infile.get() != '\n') {}

        // locus information
        locus.resize(nlocus);
        locuscount.resize(4, 0);
        for (int i = 0; i < nlocus; i++) {
            locus[i].readLocus(infile, outStream, sexlinked);
            if (locus[i].name == "") {
                ostrstream ostr;
                ostr << ++locuscount[locus[i].type % 4] << ends;
                locus[i].name = ostr.str();
            }
            locus[i].position = i;
            locus[i].chromosome = 0;
        }

        // sex difference, interference
        if (infile >> sexdiff >> interference, infile.fail()) {
            *outStream << "\nData file error: premature end of file, expecting sex-difference\n";
            exit(-1);
        }
        while (infile.peek() != EOF && infile.get() != '\n') {}

        // recombination fractions
        recombination.resize(nlocus - 1);
        for (int i = 0; i < nlocus - 1; i++) {
            double x;
            if ((infile >> x, infile.fail()) && i < nlocus - 2) {
                *outStream << "\nData file error: premature end of file, expecting recombination fraction\n";
                exit(-1);
            }
            recombination[i] = x;
        }
        while (infile.peek() != EOF && infile.get() != '\n') {}
    }

    // skip the rest of the data file
    *outStream << "done" << endl;
}

//  readmapfile */

// read a PLINK map file

void LinkageData::readmapfile(string &filename) {
    if (filename == "") {
        return;
    }

    ifstream infile(filename.c_str());
    if (!infile) {
        *outStream << "Failed to open map file " << filename << endl;
        exit(-1);
    }
    *outStream << "Reading map file " << filename << "..." << flush;

    int locusIndex = 0;
    int markercount = 0;
    string line;
    while ((line = getline(infile)) != "") {
        while (locusIndex < locus.size() && locus[locusIndex].type != 3) {
            locusIndex++;
        }
        istrstream instr(line.c_str());
        char buf[256];
        buf[0] = 0;
        char chr[256];
        chr[0] = 0;
        // read chr
        if (instr >> chr, instr.fail()) {
            *outStream << "\nMap file error: premature end of line: expecting chromosome\n";
            exit(-1);
        }
        // read name
        if (instr >> buf, instr.fail()) {
            *outStream << "\nMap file error: premature end of line: expecting name for marker " << markercount << endl;
            exit(-1);
        } else {
            if (locusIndex < locus.size()) {
                locus[locusIndex].name = buf;
                locus[locusIndex].chromosome = atoi(chr);
                // ignore distance
                instr >> buf;
                //  store position
                if (instr >> buf, instr.fail()) {
                    *outStream << "\nMap file error: premature end of line: expecting position1\n";
                    exit(-1);
                } else {
                    locus[locusIndex].position = atoi(buf);
                }
            } else {
                // insert a SNP
                locus.resize(locus.size() + 1);
                nlocus++;
                locus[locusIndex].name = buf;
                locus[locusIndex].type = 3;
                locus[locusIndex].nallele = 2;
                locus[locusIndex].chromosome = atoi(chr);
                // equal frequencies for all alleles
                locus[locusIndex].frequency.resize(locus[locusIndex].nallele);
                for (int j = 0; j < locus[locusIndex].nallele; j++) {
                    locus[locusIndex].frequency[j] = 1.0 / locus[locusIndex].nallele;
                }
                locuscount[3]++;
                // look for extended map file format
                // ignore distance
                instr >> buf;
                //  store position
                if (instr >> buf, instr.fail()) {
                    *outStream << "\nMap file error: premature end of line: expecting position2\n";
                    exit(-1);
                } else {
                    locus[locusIndex].position = atoi(buf);
                }
            }
            //  store allele data
            locus[locusIndex].allele.resize(2);
            locus[locusIndex].allele[0] = 1;
            locus[locusIndex].allele[1] = 2;
            for (int i = 0; i < 2; i++) {
                if (instr >> buf, !instr.fail()) {
                    const string letters = "NnAaCcGgTt";
                    int index = letters.find(buf[0]) / 2;
                    if (index >= 0 && index < 5) {
                        ACGT = true;
                        locus[locusIndex].allele[i] = index;
                    } else {
                        locus[locusIndex].allele[i] = atoi(buf);
                        if (locus[locusIndex].allele[i] != 0) {
                            ACGT = false;
                        }
                    }
                } // if instr
            } // for i
        }

        locusIndex++;
        markercount++;
    }

    // skip the rest of the data file
    *outStream << "done" << endl;
}



//  readphenofile */

// read a PLINK format phenotype file

void LinkageData::readphenofile(string &filename) {
    if (filename == "") {
        return;
    }

    ifstream infile(filename.c_str());
    if (!infile) {
        *outStream << "Failed to open phenotype file " << filename << endl;
        exit(-1);
    }
    *outStream << "Reading phenotype file " << filename << flush;

    *outStream << "..." << flush;
    string line;
    bool firstline = true;
    int npheno = 0;
    while ((line = getline(infile)) != "") {
        istrstream instr(line.c_str());
        char buf[256];
        buf[0] = 0;
        // read first two words
        if (instr >> buf, instr.fail()) {
            *outStream << "\nPhenotype file error: premature end of line: expecting pedigree ID\n";
            exit(-1);
        }
        string pid = buf;
        if (instr >> buf, instr.fail()) {
            *outStream << "\nPhenotype file error: premature end of line: expecting subject ID\n";
            exit(-1);
        }
        string sid = buf;
        if (firstline && pid == "FID" && sid == "IID") {
            // header line
            npheno = 0;
            while (instr >> buf, !instr.fail()) {
                Locus newLocus;
                newLocus.type = 4;
                newLocus.name = buf;
                newLocus.chromosome = 0;
                npheno++;
                newLocus.position = npheno;
                locus.push_back(newLocus);
                locuscount[0]++;
            }
        } else {
            bool haveSubject = true;
            if (pedigree.find(pid) == pedigree.end()) {
                //  *outStream << "\nPhenotype file error: pedigree not found: " << pid << endl;
                //exit(-1);
                haveSubject = false;
            }
            int subjectIndex = -1;
            for (int i = 0; i < pedigree[pid].size(); i++)
                if (pedigree[pid][i].id == sid) {
                    subjectIndex = i;
                }
            if (subjectIndex == -1) {
                // *outStream << "\nPhenotype file error: subject not found: pedigree " << pid << " subject " << sid << endl;
                //exit(-1);
                haveSubject = false;
            }
            if (haveSubject) {
                int count = 0;
                while (instr >> buf, !instr.fail()) {
                    string value = buf;
                    if (value == MISSING) {
                        pedigree[pid][subjectIndex].trait.push_back(0);
                        pedigree[pid][subjectIndex].missingtrait.push_back(true);
                    } else {
                        pedigree[pid][subjectIndex].trait.push_back(atof(value.c_str()));
                        pedigree[pid][subjectIndex].missingtrait.push_back(false);
                    }
                    if (firstline) {
                        // insert new locus
                        Locus newLocus;
                        newLocus.type = 4;
                        ostrstream ostr;
                        ostr << ++locuscount[0] << ends;
                        newLocus.name = ostr.str();
                        newLocus.chromosome = 0;
                        newLocus.position = locuscount[0];
                        locus.push_back(newLocus);
                    }
                    count++;
                }
                if (firstline) {
                    npheno = count;
                } else {
                    if (count < npheno) {
                        *outStream << "\nPhenotype file error: too few phenotypes: pedigree " << pid << " subject " << sid << endl;
                        exit(-1);
                    }
                    if (count > npheno) {
                        *outStream << "\nPhenotype file error: too many phenotypes: pedigree " << pid << " subject " << sid << endl;
                        exit(-1);
                    }
                }
            }
        }
        firstline = false;
    }

    *outStream << "done" << endl;
}



//  readLocus */

// read information for one locus
void Locus::readLocus(ifstream &infile, ostream *&outStream, int sexlinked) {
    if (infile >> type >> nallele, infile.fail()) {
        *outStream << "\nData file error: premature end of file, expecting locus type\n";
        exit(-1);
    }
    if (type != 0 && type != 1 && type != 3 && type != 4) {
        *outStream << "\nData file error: locus type read as " << type << endl;
        exit(-1);
    }

    // read the optional locus name
    char c = 0;
    while (infile.peek() != EOF && c != '#' && c != '\n') {
        c = infile.get();
    }
    if (c == '#') {
        char buf[256];
        infile >> buf;
        name = buf;
        while (infile.peek() != EOF && infile.get() != '\n') {}
    }

    // no more info for covariates
    if (type == 4) {
        return;
    }

    // allele frequencies
    frequency.resize(nallele);
    for (int i = 0; i < nallele; i++) {
        float x;
        if (infile >> x, infile.fail()) {
            *outStream << "\nData file error: premature end of file, expecting allele frequencies\n";
            exit(-1);
        }
        if (x < 0 || x > 1) {
            *outStream << "\nData file error: frequency out of range: " << type << nallele << x << endl;
            exit(-1);
        }
        frequency[i] = x;
    }
    while (infile.peek() != EOF && infile.get() != '\n') {}

    switch (type) {
    case 1: {
        // affection status
        if (infile >> nliability, infile.fail()) {
            *outStream << "\nData file error: premature end of file, expecting number of liability classes\n";
            exit(-1);
        }
        while (infile.peek() != EOF && infile.get() != '\n') {}
        penetrance.resize(nliability);
        if (sexlinked) {
            malepenetrance.resize(nliability);
        }
        for (int i = 0; i < nliability; i++) {
            penetrance[i].resize(nallele*(nallele + 1) / 2);
            for (int j = 0; j < penetrance[i].size(); j++) {
                float x;
                if (infile >> x, infile.fail()) {
                    *outStream << "\nData file error: premature end of file, expecting penetrances\n";
                    exit(-1);
                }
                if (x < 0 || x > 1) {
                    *outStream << "\nData file error: penetrance out of range: " << x << endl;
                    exit(-1);
                }
                penetrance[i][j] = x;
            }
            while (infile.peek() != EOF && infile.get() != '\n') {}
            if (sexlinked) {
                malepenetrance[i].resize(nallele);
                for (int j = 0; j < nallele; j++) {
                    float x;
                    if (infile >> x, infile.fail()) {
                        *outStream << "\nData file error: premature end of file, expecting male penetrances\n";
                        exit(-1);
                    }
                    if (x < 0 || x > 1) {
                        *outStream << "\nData file error: male penetrance out of range: " <<
                                   x << endl;
                    }
                    malepenetrance[i][j] = x;
                }
                while (infile.peek() != EOF && infile.get() != '\n') {}
            }
        }
        break;
    }
    case 0: {
        // quantitative trait
        if (infile >> nliability, infile.fail()) {
            *outStream << "\nData file error: premature end of file, expecting number of quantitative variables\n";
            exit(-1);
        }
        while (infile.peek() != EOF && infile.get() != '\n') {}
        penetrance.resize(nliability); // actually the genotypic means
        if (sexlinked) {
            malepenetrance.resize(nliability);
        }
        for (int i = 0; i < nliability; i++) {
            penetrance[i].resize(nallele*(nallele + 1) / 2);
            for (int j = 0; j < penetrance[i].size(); j++) {
                float x;
                if (infile >> x, infile.fail()) {
                    *outStream << "\nData file error: premature end of file, expecting genotypic means\n";
                    exit(-1);
                }
                penetrance[i][j] = x;
            }
            while (infile.peek() != EOF && infile.get() != '\n') {}
            if (sexlinked) {
                malepenetrance[i].resize(nallele);
                for (int j = 0; j < nallele; j++) {
                    float x;
                    if (infile >> x, infile.fail()) {
                        *outStream << "\nData file error: premature end of file, expecting male allelic means\n";
                        exit(-1);
                    }
                    malepenetrance[i][j] = x;
                }
                while (infile.peek() != EOF && infile.get() != '\n') {}
            }
        }
        // variance-covariance matrix
        variance.resize(nliability);
        for (int i = 0; i < variance.size(); i++) {
            variance[i].resize(nliability);
            for (int j = i; j < variance.size(); j++) {
                float x;
                if (infile >> x, infile.fail()) {
                    *outStream << "\nData file error: premature end of file, expecting quantitative trait covariances\n";
                    exit(-1);
                }
                variance[i][j] = x;
                variance[j][i] = x;
            }
        }
        while (infile.peek() != EOF && infile.get() != '\n') {}
        if (infile >> hetmultiplier, infile.fail()) {
            *outStream << "\nData file error: premature end of file, expecting heterozygote multiplier\n";
            exit(-1);
        }
        while (infile.peek() != EOF && infile.get() != '\n') {}
    }

    case 3:
    {}
    default:
    {}
    }
}



//  defaultdatafile */

void LinkageData::defaultdatafile(string &filename) {

    ifstream infile(filename.c_str());
    if (!infile) {
        *outStream << "Failed to open pedigree file " << filename << endl;
        exit(-1);
    }
    string line;
    int ncol = 0;
    bool qt = false;
    vector<set<int> > allele;
    while ((line = getline(infile)) != "") {
        // count the number of columns up to any comments
        istrstream istr(line.c_str());
        char buf[256];
        buf[0] = 0;
        int j = 0;
        while ((istr >> buf, !istr.fail()) && buf[0] != '<' && buf[1] != '<') {
            // count the number of alleles at each marker
            if (j > 5) {
                if (allele.size() < j / 2 - 2) {
                    allele.resize(j / 2 - 2);
                }
                allele[j/2-3].insert(atoi(buf));
            }
            // check for quantitative trait in column 6
            if (j == 5) {
                for (int i = 0; i < strlen(buf); i++)
                    if (buf[i] != '0' + UNSURE && buf[i] != '0' + UNAFFECTED &&
                            buf[i] != '0' + AFFECTED) {
                        qt = true;
                    }
            }
            j++;
        }
        if (ncol && j != ncol) {
            *outStream << "\nPedigree file error: this line has " << (j < ncol ? "fewer" : "more") << " columns than the previous line\n" << line.c_str() << endl;
            exit(-1);
        } else {
            ncol = j;
        }
    }
    if (ncol % 2) {
        *outStream << "\nPedigree file error: could not determine number of markers\n" << "There are " << ncol << " columns, should be an even number" << endl;
        exit(-1);
    }

    nlocus = (ncol - 6) / 2 + 1; // assume one trait locus followed by markers
    risklocus = 0;
    sexlinked = 0;
    program = 5;
    mutlocus = 0;
    malemut = 0;
    femalemut = 0;
    hapfreq = 0;
    order.push_back(1);
    sexdiff = 0;
    interference = 0;
    locus.resize(nlocus);
    locuscount.resize(4, 0);
    locus[0].chromosome = 0;
    locus[0].position = 1;
    // column 6 is either a quantitative trait...
    if (qt) {
        locus[0].name = "Trait";
        locus[0].type = 0;
        locus[0].nallele = 0;
    }
    // ...or an affection status
    else {
        locus[0].name = "Disease";
        locus[0].type = 1;
        locus[0].nallele = 2;
        locus[0].frequency.push_back(0.999);
        locus[0].frequency.push_back(0.001);
        locus[0].nliability = 1;
        locus[0].penetrance.resize(1);
        locus[0].penetrance[0].push_back(0);
        locus[0].penetrance[0].push_back(1);
        locus[0].penetrance[0].push_back(1);
    }
    locuscount[locus[0].type]++;
    for (int i = 1; i < nlocus; i++) {
        allele[i-1].erase(0);
        ostrstream ostr;
        ostr << i << ends;
        locus[i].name = ostr.str();
        locus[i].type = 3;
        locus[i].position = i;
        locus[i].chromosome = 0;
        // not using the number of observed alleles
        //locus[i].nallele=allele[i-1].size();
        // using the maximum allele value
        locus[i].nallele = *(max_element(allele[i-1].begin(), allele[i-1].end()));
        // equal frequencies for all alleles
        locus[i].frequency.resize(locus[i].nallele);
        for (int j = 0; j < locus[i].nallele; j++) {
            locus[i].frequency[j] = 1.0 / locus[i].nallele;
        }
        locuscount[3]++;
    }
    recombination.resize(nlocus - 1, 0);
}



//  writedatafile */

// write a Linkage format data file

void LinkageData::writedatafile() {
    // first line of data file
    *outStream << nlocus << " " <<  risklocus << " " <<  sexlinked << " " <<  program << endl;

    // second line of data file
    *outStream << mutlocus << " " <<  malemut << " " <<  femalemut << " " << hapfreq << endl;

    // third line: order of loci
    for (int i = 1; i <= nlocus; i++) {
        *outStream << i << " ";
    }

    *outStream << endl;

    // locus information
    locus.resize(nlocus);
    for (int i = 0; i < nlocus; i++) {
        *outStream << locus[i].type << " " << locus[i].nallele << endl;
        for (int j = 0; j < locus[i].frequency.size(); j++) {
            *outStream << locus[i].frequency[j] << " ";
        }
        *outStream << endl;
        switch (locus[i].type) {
        case 0: {
            for (int j = 0; j < locus[i].penetrance.size(); j++)
                for (int k = 0; k < locus[i].penetrance[j].size(); k++) {
                    *outStream << locus[i].penetrance[j][k] << " ";
                }
            *outStream << endl;
            for (int j = 0; j < locus[i].variance.size(); j++)
                for (int k = 0; k < locus[i].variance[j].size(); k++) {
                    *outStream << locus[i].variance[j][k] << " ";
                }
            *outStream << endl;
            *outStream << locus[i].hetmultiplier << endl;
            break;
        }
        case 1: {
            *outStream << locus[i].nliability << endl;
            for (int j = 0; j < locus[i].nliability; j++) {
                for (int k = 0; k < locus[i].penetrance[j].size(); k++) {
                    *outStream << locus[i].penetrance[j][k] << " ";
                }
                *outStream << endl;
                if (sexlinked) {
                    for (int k = 0; k < locus[i].malepenetrance[j].size(); k++) {
                        *outStream << locus[i].malepenetrance[j][k] << " ";
                    }
                    *outStream << endl;
                }
            }
            break;
        }
        }
    }

    // sex difference, interference
    *outStream << sexdiff << " " << interference << endl;

    // recombination fractions
    for (int i = 0; i < recombination.size(); i++) {
        *outStream << recombination[i] << " ";
    }
    *outStream << endl;

}



//  makehashtables */

void LinkageData::makehashtables() {
    int ntrait = 0;
    int naffect = 0;
    int nmarker = 0;
    for (vector<Locus>::iterator i = locus.begin(); i != locus.end(); i++) {
        if (i->type % 4 == 0) {
            traithash[i->name] = ntrait++;
        }
        if (i->type == 1) {
            diseasehash[i->name] = naffect++;
        }
        if (i->type == 3) {
            markerhash[i->name] = nmarker++;
        }
    }
}


