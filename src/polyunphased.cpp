/* polyunphased.cpp - main program for POLYUNPHASED

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

#include <iostream>
#include "UnphasedOptions.h"
#include "UnphasedAnalysis.h"

using namespace std;

int main(int argc, char **argv) {

    const string VERSION = "1.2";

    UnphasedOptions options(argc, argv);
    ostream *outStream;
    if (options.outputFilename != "") {
        outStream = (ostream *)(new ofstream(options.outputFilename.c_str()));
        if (!((ofstream *)outStream)->is_open()) {
            cout << endl << "ERROR: could not open output file " << options.outputFilename << endl;
            exit(-1);
        }
    } else {
        outStream = &cout;
    }
    *outStream << endl << "***START OF POLYUNPHASED " << VERSION << "***" << endl;
    UnphasedAnalysis analysis(outStream);
    analysis.run(options);
    *outStream << "***END OF POLYUNPHASED " << VERSION << "***" << endl << endl;;
    if (outStream != &cout) {
        ((ofstream *)outStream)->close();
    }
    return(0);
}
