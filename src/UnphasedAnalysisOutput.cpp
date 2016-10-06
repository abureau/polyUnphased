/* UnphasedAnalysisOutput.cpp - Association analysis in trios and unrelateds

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
#include <iomanip>
#include <strstream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

// local includes
#include "stl.h"
#include "UnphasedOptions.h"
#include "UnphasedAnalysis.h"
#include "pvalues.h"
#define COMMENT

void UnphasedAnalysis::outputTabularHeaders(UnphasedOptions &options) {
#ifdef COMMENT
  // headers for SNPs
  // not much use in a file with a mixture of SNPs and multiallelic markers
  // Code de dÃ©buggage
  //cout << "window " << options.window << " genoCode.size " << genoCode.size() << endl;
  //if (1 == options.window && (genoCode.size() + ! options.genotype == 3)) {
  // contournement d'un bug d'Unphased
  if (1 == options.window && options.condition.size() == 0 && ! options.genotype) {
	  // binary traits
	  if ("binary" == typeOfPhenotype) { // binary trait headers
            if (haveFamilies) {
	      if (! options.genotype) {
	      tabularFamilies << right << setw(4) << "CHR" << setw(12) << "SNP" << setw(13) << "BP" << setw(4) << "A1" << setw(4) << "A2" << setw(7) <<"T" << setw(7) << "U" << setw(13) << "OR" << setw(13) << "L95" << setw(13) << "U95" << setw(13) << "CHISQ" << setw(13) << "P" << left << endl;
	      } else {
		tabularFamilies << right << setw(4) << "CHR" << setw(12) << "SNP" << setw(4) << "A1" << setw(4) << "A2" << setw(9) << "TEST" << setw(15) << "AFF" << setw(15) << "UNAFF" << setw(13) << "CHISQ" << setw(5) << "DF" << setw(13) << "P" << left << endl;
	      }
	    }
            if (haveUnrelateds) {
	      if (! options.genotype) {
		tabularUnrelateds << right << setw (4) << "CHR" << setw(12) << "SNP" << setw(13) << "BP" << setw(4) << "A1" << setw(13) << "F_A" << setw(13) << "F_U" << setw(4) << "A2" << setw(13) << "CHISQ" << setw(13) << "P" << setw(13) << "OR" << setw(13) << "SE" << setw(13) << "L95" << setw(13) << "U95" << left << endl;
	      } else {
		tabularUnrelateds << right << setw(4) << "CHR" << setw(12) << "SNP" << setw(4) << "A1" << setw(4) << "A2" << setw(9) << "TEST" << setw(15) << "AFF" << setw(15) << "UNAFF" << setw(13) << "CHISQ" << setw(5) << "DF" << setw(13) << "P" << endl;
	      }
	    }
	  } else {	
	  if (typeOfPhenotype == "polytomous")  
	  	{ 
	  	// polytomous trait headers
            if (haveFamilies) {
	      if (! options.genotype) {
	      tabularFamilies << right << setw(4) << "CHR" << setw(12) << "SNP" << setw(13) << "BP"  << setw(4) << "A1" << setw(4) << "A2" ;
          for (int k = 1; k <= K; k++)	
	      	tabularFamilies	<< setw(6) << "T" << k << setw(6) << "U" << k ;
          for (int k = 1; k <= K-1; k++) {	
	      	tabularFamilies << setw(12) << "OR" << k << setw(12) << "L95_" << k << setw(12) << "U95_" << k ; 
			if (options.individual) tabularFamilies << setw(12) << "CHISQ" << k << setw(12) << "P" << k;
			}
	      tabularFamilies << setw(13) << "CHISQ" << setw(13) << "P" << left << endl;
	      }
	    }
	    } else {	  	 
	    // quantitative trait headers
            if (haveFamilies) {
	      if (! options.genotype) {
	      tabularFamilies << right << setw(4) << "CHR" << setw(12) << "SNP" << setw(13) << "BP" << setw(4) << "A1" << setw(4) << "A2" << setw(7) << "COUNT" << setw(13) << "BETA" << setw(13) << "SE" << setw(13) << "L95" << setw(13) << "U95" << setw(13) << "CHISQ" << setw(13) << "P" << left << endl;
	      } else {
	      tabularFamilies << right << setw(4) << "CHR" << setw(12) << "SNP" << setw(4) << "A1" << setw(4) << "A2" << setw(9) << "TEST" << setw(15) << "COUNT" << setw(13) << "CHISQ" << setw(5) << "DF" << setw(13) << "P" << left << endl;
	      }
	    }
	    if (haveUnrelateds) {
	      if (! options.genotype) {
		tabularUnrelateds << right << setw(4) << "CHR" << setw(12) << "SNP" << setw(13) << "BP" << setw(4) << "A1" << setw(13) << "F_ALL" << setw(4) << "A2" << setw(13) << "BETA" << setw(13) << "SE" << setw(13) << "L95" << setw(13) << "U95" << setw(13) << "CHISQ" << setw(13) << "P" << left << endl;
	      } else {
	      tabularUnrelateds << right << setw(4) << "CHR" << setw(12) << "SNP" << setw(4) << "A1" << setw(4) << "A2" << setw(9) << "TEST" << setw(15) << "COUNT" << setw(13) << "CHISQ" << setw(5) << "DF" << setw(13) << "P" << endl;
	      }
	    }
	  }
  }
  } else {
#endif
 // headers for haplotypes or multiallelic markers - and also SNPs
    bool alleles = (options.window + options.condition.size() + options.tag.size() == 1);
    string testedObject=(options.genotype?"GENOTYPE":(alleles?"ALLELE":"HAPLOTYPE")); 
    string oneMarker=(options.window==1?"MARKER":"MARKERS");

    if ("binary" == typeOfPhenotype) {
      // binary trait headers
      if (haveFamilies) {
	tabularFamilies << right << setw(4) << "CHR" << setw(9) << oneMarker << setw(13) << "BP_START" << setw(13) << "BP_END" << setw(13) << testedObject << setw(13) << "F_A" << setw(13) << "F_U" << setw(13) << "CHISQ" << setw(5) << "DF" << setw(13) << "P" << setw(13) << "OR" << setw(13) << "SE" << setw(13) << "L95" << setw(13) << "U95" << left << endl;
      }
      if (haveUnrelateds) {
	tabularUnrelateds << right << setw(4) << "CHR" << setw(9) << oneMarker << setw(13) << "BP_START" << setw(13) << "BP_END" << setw(13) << testedObject << setw(13) << "F_A" << setw(13) << "F_U" << setw(13) << "CHISQ" << setw(5) << "DF" << setw(13) << "P" << setw(13) << "OR" << setw(13) << "SE" << setw(13) << "L95" << setw(13) << "U95" << left << endl;
      }
    } else {
	  if (typeOfPhenotype == "polytomous")  
	  	{ 
	  	// polytomous trait headers
            if (haveFamilies) {
	      if (! options.genotype) {
	      tabularFamilies << right << setw(4) << "CHR" << setw(9) << oneMarker << setw(13) << "BP_START" << setw(13) << "BP_END" << setw(13) << testedObject ;
	      for (int k = 1; k <= K; k++)	
	      	tabularFamilies << setw(10) << "F" << k << "_A" << setw(10) << "F" << k << "_U";
	      for (int k = 1; k < K; k++)	
	      	tabularFamilies  << setw(12) << "CHISQ" << k << setw(5) << "DF" << setw(12) << "P" << k << setw(12) << "OR" << k << setw(12) << "SE" << k << setw(12) << "L95_" << k << setw(12) << "U95_" << k; 
	      tabularFamilies << left << endl;
	      }
	    }
	    } else {
      // quantitative trait headers
      if (haveFamilies) {
	tabularFamilies << right << setw(4) << "CHR" << setw(9) << oneMarker << setw(13) << "BP_START" << setw(13) << "BP_END" << setw(13) << testedObject << setw(13) << "F_ALL" << setw(13) << "CHISQ" << setw(5) << "DF" << setw(13) << "P" << setw(13) << "BETA" << setw(13) << "SE" << setw(13) << "L95" << setw(13) << "U95" << left << endl;
      }
      if (haveUnrelateds) {
	tabularUnrelateds << right << setw(4) << "CHR" << setw(9) << oneMarker << setw(13) << "BP_START" << setw(13) << "BP_END" << setw(13) << testedObject << setw(13) << "F_ALL" << setw(13) << "CHISQ" << setw(5) << "DF" << setw(13) << "P" << setw(13) << "BETA" << setw(13) << "SE" << setw(13) << "L95" << setw(13) << "U95" << left << endl;
      }
    }
}
}    
}

void UnphasedAnalysis::outputTabular(UnphasedOptions &options, vector<int>& combination, double null, double alternative, int df) {
    int     i, j, liIdx, liRef;
    double  laTotalFamilyCount[8] = {0, 0, 0, 0, 0, 0, 0, 0}, laTotalUnrelatedCount[2] = {0, 0}, lrFU, lrFA, lrExpected, lrChi2[3] = {0,0,0}, lrP[3]= {1,1,1}, lrOdds[3] = {1,1,1}, lrBeta, lrT, lrR2;
    string  lsHap1, lsHap2;
    stringstream lsMarker, lsRatios;
    Locus   loLocus1, loLocus2;

    // sort haplotypes into alphabetical order
    vector<Haplotype> sortedHaps;
    sortedHaps = genoCode;
    sort(sortedHaps.begin(), sortedHaps.end());

    //  Calculate the total frequencies.
    int nCount = 2;
    if (typeOfPhenotype == "polytomous") nCount = 2*K;
    for (j = 0; j < genoCode.size(); j++) {
            if (!zero[j]) {
    		for (i = 0; i < nCount; i++) 
                laTotalFamilyCount[i]    += familyCount[i][j];
    		for (i = 0; i < 2; i++) 
                laTotalUnrelatedCount[i] += unrelatedCount[i][j];
            }
        }

#ifdef COMMENT
    //  Do the output for SNPs
    // Not much use in a file with a mixutre of SNPs and multiallelic markers
    if (1 == options.window && sortedHaps.size() + ! options.genotype == 3) {
 
            //  Make the marker code.
            lsMarker << options.marker[combination[0]];
            loLocus1 = getLocus(options.marker[combination[0]]);

            switch (genoCode.size()) {
            case 0:
                return;
            case 1:
                /*
                            lsHap1 = genoCode[0].str(options.condition.size() * options.condgenotype, options.genotype, ACGT);
                            lsHap2 = "NA";
                            lrFA   = -1;
                            lrFU   = -1;
                            out << right;
                            out << "   1" << setw(12) << lsMarker.str() << setw(11) << siPosition;
                            siPosition++;
                            out << setw(5) << lsHap1 << setw(9) << lrFA;
                            if (lrFU > -0.5) {
                                out << setw(9) << lrFU << setw(5) << lsHap2;
                                out << setw(13) << lrChi2 << setw(13) << lrP << setw(13) << lrOdds;
                            } else {
                                Out << "       NA   NA           NA           NA           NA";
                            }
                            out << left << endl;
                */
                break;
            default:
	      // allele tests
	      if (! options.genotype) {
	      // one line for each non-reference allele
	      for (i = 0; i < sortedHaps.size(); i++) {
                    liIdx  = sortedHaps[i].index(genoCode);
		    liRef  = reference.index(genoCode);
		    if (liIdx != liRef) {
		    lsHap1 = sortedHaps[i].str(options.condition.size() * options.condgenotype, options.genotype, ACGT);
		    lsHap2 = reference.str(options.condition.size() * options.condgenotype, options.genotype, ACGT);
                    if (options.individual) {
                        lrChi2[0] = chisq[liIdx];
                        lrP[0]    = pvalue[liIdx];
                    	if ("polytomous" == typeOfPhenotype) {
                        	for (int k = 1; k < K-1; k++) {
                        		lrChi2[k] = chisq[liIdx + k*genoCode.size()];
                        		lrP[k]    = pvalue[liIdx + k*genoCode.size()];
                        		}
                        	}
                    } else if (sortedHaps.size() == 2) {
		      // for a SNP, global test = individual test
                        lrChi2[0] = 2 * (alternative - null);
                        lrP[0]    = pochisq(2 * (alternative - null), df);
                    }
                    
                    lrOdds[0] = exp(beta[liIdx] - beta[liRef]);
                    if ("polytomous" == typeOfPhenotype) {
                        for (int k = 1; k < K-1; k++)
                    		lrOdds[k] = exp(beta[liIdx + k*genoCode.size()] - beta[liRef + k*genoCode.size()]);
                    	}
		    lrBeta = beta[liIdx] - beta[liRef];
		    // allele tests in families
                    if (haveFamilies) {
                        tabularFamilies << right;
                        tabularFamilies << setw(4) << loLocus1.chromosome << setw(12) << lsMarker.str() << setw(13) << loLocus1.position << setw(4) << lsHap1 << setw(4) << lsHap2;
			if ("binary" == typeOfPhenotype) {
			  tabularFamilies << setw(7) << familyCount[1][liIdx] << setw(7) << familyCount[0][liIdx];
			  tabularFamilies << setw(13) << lrOdds[0] << setw(13) << exp(log(lrOdds[0])-1.96*stderror[liIdx]) << setw(13) << exp(log(lrOdds[0])+1.96*stderror[liIdx]) << setw(13) << lrChi2[0] << setw(13) << lrP[0];
			} else if ("polytomous" == typeOfPhenotype) {
			  for (int k = 0; k < K; k++)
			  	tabularFamilies << setw(7) << familyCount[k][liIdx] << setw(7) << familyCount[k+K][liIdx] ;
              for (int k = 0; k < K-1; k++)	{
			  	tabularFamilies << setw(13) << lrOdds[k] << setw(13) << exp(log(lrOdds[k])-1.96*stderror[liIdx + k*genoCode.size()]) << setw(13) << exp(log(lrOdds[k])+1.96*stderror[liIdx + k*genoCode.size()]);
				if (options.individual)
					tabularFamilies << setw(13) << lrChi2[k] << setw(13) << lrP[k];	
				}
			if (sortedHaps.size() == 2) tabularFamilies << setw(13) << lrChi2[0] << setw(13) << lrP[0];
			} else {
			  tabularFamilies << setw(7) << familyCount[0][liIdx];
			  tabularFamilies << setw(13) << lrBeta << setw(13) << stderror[liIdx] << setw(13) << lrBeta-1.96*stderror[liIdx] << setw(13) << lrBeta+1.96*stderror[liIdx] << setw(13) << lrChi2[0] << setw(13) << lrP[0];
			}
                        tabularFamilies << left << endl;
                    }
		    // allele tests in unrelateds
                    if (haveUnrelateds) {
                        tabularUnrelateds << right;
                        tabularUnrelateds << setw(4) << loLocus1.chromosome << setw(12) << lsMarker.str() << setw(13) << loLocus1.position;
                        tabularUnrelateds << setw(4) << lsHap1;
			if ("binary" == typeOfPhenotype) {
			  tabularUnrelateds << setw(13) << (unrelatedCount[1][liIdx] / laTotalUnrelatedCount[1]);
			  tabularUnrelateds << setw(13) << (unrelatedCount[0][liIdx] / laTotalUnrelatedCount[0]);
			  tabularUnrelateds << setw(4) << lsHap2;
			  tabularUnrelateds << setw(13) << lrChi2 << setw(13) << lrP << setw(13) << lrOdds[0] << setw(13) << stderror[liIdx] << setw(13) << exp(log(lrOdds[0])-1.96*stderror[liIdx]) << setw(13) << exp(log(lrOdds[0])+1.96*stderror[liIdx]);
			} else {
			  tabularUnrelateds << setw(13) << (unrelatedCount[0][liIdx] / laTotalUnrelatedCount[0]);
			  tabularUnrelateds << setw(4) << lsHap2;
			  tabularUnrelateds << setw(13) << lrBeta << setw(13) << stderror[liIdx] << setw(13) << lrBeta-1.96*stderror[liIdx] << setw(13) << lrBeta+1.96*stderror[liIdx] << setw(13) << lrChi2 << setw(13) << lrP;
			}
                        }
		    tabularUnrelateds << left << endl;
                    }
	      }
	      } else { // genotype tests for SNPs
		// get alleles from the het genotype
		lsHap1 = sortedHaps[1].str(options.condition.size() * options.condgenotype, options.genotype, ACGT)[2];
		lsHap2 = sortedHaps[1].str(options.condition.size() * options.condgenotype, options.genotype, ACGT)[0];
		// get indices of the three genotypes
		int liHom1 = sortedHaps[2].index(genoCode);
		int liHet = sortedHaps[1].index(genoCode);
		int liHom2 = sortedHaps[0].index(genoCode);
		// genotype tests in families
		if (haveFamilies) {
		      // omnibus genotype test
                        tabularFamilies << right;
                        tabularFamilies << setw(4) << loLocus1.chromosome << setw(12) << lsMarker.str() << setw(4) << lsHap1 << setw(4) << lsHap2 << setw(9) << "GENO";
			if ("binary" == typeOfPhenotype) {
			  lsRatios.str("");
			  lsRatios << familyCount[1][liHom1] << "/" << familyCount[1][liHet] << "/" << familyCount[1][liHom2];
			  tabularFamilies << setw(15) << lsRatios.str();
			}
                        lsRatios.str("");
                        lsRatios << familyCount[0][liHom1] << "/" << familyCount[0][liHet] << "/" << familyCount[0][liHom2];
                        tabularFamilies << setw(15) << lsRatios.str() << setw(13) << 2 *(alternative - null) << setw(5) << df << setw(13) << pochisq(2 *(alternative - null), df) << endl;
			if (options.individual) {
			// dominant model
			tabularFamilies << setw(4) << loLocus1.chromosome << setw(12) << lsMarker.str() << setw(4) << lsHap1 << setw(4) << lsHap2 << setw(9) << "DOM";
			if ("binary" == typeOfPhenotype) {
			  lsRatios.str("");
			  lsRatios << familyCount[1][liHom1] + familyCount[1][liHet] << "/" << familyCount[1][liHom2];
			  tabularFamilies << setw(15) << lsRatios.str();
			}
                        lsRatios.str("");
                        lsRatios << familyCount[0][liHom1] + familyCount[0][liHet] << "/" << familyCount[0][liHom2];
                        tabularFamilies << setw(15) << lsRatios.str() << setw(13) << chisq[liHom2] << setw(5) << 1 << setw(13) << pvalue[liHom2] << endl;
			// recessive model
                        tabularFamilies << setw(4) << loLocus1.chromosome << setw(12) << lsMarker.str() << setw(4) << lsHap1 << setw(4) << lsHap2 << setw(9) << "REC";
			if ("binary" == typeOfPhenotype) {
			  lsRatios.str("");
			  lsRatios << familyCount[1][liHom1] << "/" <<  familyCount[1][liHet] + familyCount[1][liHom2];
			  tabularFamilies << setw(15) << lsRatios.str();
			}
                        lsRatios.str("");
                        lsRatios << familyCount[0][liHom1] << "/" <<  familyCount[0][liHet] + familyCount[0][liHom2];
                        tabularFamilies << setw(15) << lsRatios.str() << setw(13) << chisq[liHom1] << setw(5) << 1 << setw(13) << pvalue[liHom1] << endl;
			}
		}
		// genotype tests in unrelateds
                    if (haveUnrelateds) {
		      // omnibus genotype tests
                        tabularUnrelateds << right;
                        tabularUnrelateds << setw(4) << loLocus1.chromosome << setw(12) << lsMarker.str() << setw(4) << lsHap1 << setw(4) << lsHap2 << setw(9) << "GENO";
			if ("binary" == typeOfPhenotype) {
			  lsRatios.str("");
			  lsRatios << unrelatedCount[1][liHom1] << "/" << unrelatedCount[1][liHet] << "/" << unrelatedCount[1][liHom2];
			  tabularUnrelateds << setw(15) << lsRatios.str();
			}
			lsRatios.str("");
                        lsRatios << unrelatedCount[0][liHom1] << "/" << unrelatedCount[0][liHet] << "/" << unrelatedCount[0][liHom2];
                        tabularUnrelateds << setw(15) << lsRatios.str() << setw(13) << 2 *(alternative - null) << setw(5) << df << setw(13) << pochisq(2 *(alternative - null), df) << endl;
			if (options.individual) {
			// dominant model
                        tabularUnrelateds << setw(4) << loLocus1.chromosome << setw(12) << lsMarker.str() << setw(4) << lsHap1 << setw(4) << lsHap2 << setw(9) << "DOM";
			if ("binary" == typeOfPhenotype) {
			  lsRatios.str("");
			  lsRatios << unrelatedCount[1][liHom1] + unrelatedCount[1][liHet] << "/" << unrelatedCount[1][liHom2];
			  tabularUnrelateds << setw(15) << lsRatios.str();
			}
                        lsRatios.str("");
                        lsRatios << unrelatedCount[0][liHom1] + unrelatedCount[0][liHet] << "/" << unrelatedCount[0][liHom2];
                        tabularUnrelateds << setw(15) << lsRatios.str() << setw(13) << chisq[liHom2] << setw(5) << 1 << setw(13) << pvalue[liHom2] << endl;
			// recessive model
                        tabularUnrelateds << setw(4) << loLocus1.chromosome << setw(12) << lsMarker.str() << setw(4) << lsHap1 << setw(4) << lsHap2 << setw(9) << "REC";
			if ("binary" == typeOfPhenotype) {
			  lsRatios.str("");
			  lsRatios << unrelatedCount[1][liHom1] << "/" <<  unrelatedCount[1][liHet] + unrelatedCount[1][liHom2];
			  tabularUnrelateds << setw(15) << lsRatios.str();
			}
                        lsRatios.str("");
                        lsRatios << unrelatedCount[0][liHom1] << "/" <<  unrelatedCount[0][liHet] + unrelatedCount[0][liHom2];
                        tabularUnrelateds << setw(15) << lsRatios.str() << setw(13) << chisq[liHom1] << setw(5) << 1 << setw(13) << pvalue[liHom1] << endl;
			}
		    }
	      }
	      }
        } else {
#endif
      // Do the output for haplotypes or multi-allelic markers - and SNPs
      static int siFamilyPosition = 1;
      static int siUnrelatedsPosition = 1;
            //  Make the marker code.
            for (i = 0; i < options.window; i++) {
                lsMarker << options.marker[combination[i]];
            }
            loLocus1 = getLocus(options.marker[combination[0]]);
            loLocus2 = getLocus(options.marker[combination[options.window-1]]);

            switch (genoCode.size()) {
            case 0:
                return;
            case 1:
                /*
                            lsHap1 = genoCode[0].str(options.condition.size() * options.condgenotype, options.genotype, ACGT);
                            lsHap2 = "NA";
                            lrFA   = -1;
                            lrFU   = -1;
                            out << right;
                            out << "   1" << setw(12) << lsMarker.str() << setw(11) << siPosition;
                            siPosition++;
                            out << setw(5) << lsHap1 << setw(9) << lrFA;
                            if (lrFU > -0.5) {
                                out << setw(9) << lrFU << setw(5) << lsHap2;
                                out << setw(13) << lrChi2 << setw(13) << lrP << setw(13) << lrOdds;
                            } else {
                                out << "       NA   NA           NA           NA           NA";
                            }
                            out << left << endl;
                */
                break;
            default:
                for (i = 0; i < sortedHaps.size(); i++) {
		  string markerString=" ";
		  for (j = 0; j < options.window; j++) {
		    markerString += options.marker[combination[j]];
		    if (j < options.window - 1) {
		      markerString += "|";
		    }
		  }
                    liIdx  = sortedHaps[i].index(genoCode);
                    liRef  = reference.index(genoCode);
                    lsHap1 = sortedHaps[i].str(options.condition.size() * options.condgenotype, options.genotype, ACGT);
                    if (options.individual) {
                        lrChi2[0] = chisq[liIdx];
                        lrP[0]    = pvalue[liIdx];
		                if ("polytomous" == typeOfPhenotype) {
                        	for (int k = 1; k < K-1; k++) {
                        		lrChi2[k] = chisq[liIdx + k*genoCode.size()];
                        		lrP[k]    = pvalue[liIdx + k*genoCode.size()];
                        		}
                        	}                        	
                    } 
               		lrOdds[0] = exp(beta[liIdx] - beta[liRef]);     
               		if ("polytomous" == typeOfPhenotype) {
                        for (int k = 1; k < K-1; k++)
                    		lrOdds[k] = exp(beta[liIdx + k*genoCode.size()] - beta[liRef + k*genoCode.size()]);
                    	}
			    lrBeta = beta[liIdx] - beta[liRef];
		    if (haveFamilies) {
		      if (0 == i) {
			tabularFamilies << right << setw(4) << loLocus1.chromosome << setw(9) << markerString << setw(13) << loLocus1.position << setw(13) << loLocus2.position << setw(13) << "OMNIBUS" << setw(13) << "NA";
			if ("binary" == typeOfPhenotype || "polytomous" == typeOfPhenotype) {
			  tabularFamilies << setw(13) << "NA";
			}
			tabularFamilies << setw(13) << 2 *(alternative - null) << setw(5) << df << setw(13) << pochisq(2 *(alternative - null), df)
					  << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA";
			tabularFamilies << left << endl;
		      }
		      tabularFamilies << right << setw(4) << loLocus1.chromosome << setw(9) << markerString << setw(13) << loLocus1.position << setw(13) << loLocus2.position <<  setw(13) << lsHap1;
		      if ("binary" == typeOfPhenotype) {
			tabularFamilies << setw(13) << familyCount[1][liIdx]/laTotalFamilyCount[1] << setw(13) << familyCount[0][liIdx]/laTotalFamilyCount[0];
			if (options.individual) {
			  tabularFamilies << setw(13) << lrChi2[0] << setw(5) << 1 << setw(13) << lrP[0] << setw(13) << lrOdds[0] << setw(13) << stderror[liIdx] << setw(13) << exp(log(lrOdds[0])-1.96*stderror[liIdx]) << setw(13) << exp(log(lrOdds[0])+1.96*stderror[liIdx]);
			} else {
			  tabularFamilies << setw(13) << "NA" << setw(5) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA";
			}
		      } else if ("polytomous" == typeOfPhenotype) {
			  for (int k = 0; k < K; k++)
			  	tabularFamilies << setw(13) << familyCount[k][liIdx]/laTotalFamilyCount[k] << setw(13) << familyCount[k+K][liIdx]/laTotalFamilyCount[k+K] ;
				if (options.individual) {
			  		for (int k = 0; k < K-1; k++)
			  			tabularFamilies << setw(13) << lrChi2[k] << setw(5) << 1 << setw(13) << lrP[k] << setw(13) << lrOdds[k] << setw(13) << stderror[liIdx + k*genoCode.size()] << setw(13) << exp(log(lrOdds[k])-1.96*stderror[liIdx + k*genoCode.size()]) << setw(13) << exp(log(lrOdds[k])+1.96*stderror[liIdx + k*genoCode.size()]);				
				}
				else {
			  		for (int k = 0; k < K-1; k++)
						tabularFamilies << setw(13) << "NA" << setw(5) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA" ;
					}
		    	}
		    else {
			tabularFamilies << setw(13) << familyCount[0][liIdx]/laTotalFamilyCount[0];
			if (options.individual) {
			  tabularFamilies << setw(13) << lrChi2[0] << setw(5) << 1 << setw(13) << lrP[0] << setw(13) << lrBeta << setw(13) << stderror[liIdx] << setw(13) << lrBeta-1.96*stderror[liIdx] << setw(13) << lrBeta+1.96*stderror[liIdx];
			} else {
			  tabularFamilies << setw(13) << "NA" << setw(5) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA";
			}
		      }
		      tabularFamilies << left << endl;
                    }
                    if (haveUnrelateds) {
                        if (0 == i) {
			  tabularUnrelateds << right << setw(4) << loLocus1.chromosome << setw(9) << markerString << setw(13) << loLocus1.position << setw(13) << loLocus2.position << setw(13) << "OMNIBUS" << setw(13) << "NA";
			  if ("binary" == typeOfPhenotype) {
			    tabularUnrelateds << setw(13) << "NA";
			  }
			  tabularUnrelateds << setw(13) << 2 *(alternative - null) << setw(5) << df << setw(13) << pochisq(2 *(alternative - null), df)
						<< setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA";
			  tabularUnrelateds << left << endl;
                        }
			tabularUnrelateds << right << setw(4) << loLocus1.chromosome << setw(9) << markerString << setw(13) << loLocus1.position << setw(13) << loLocus2.position << setw(13) << lsHap1;
			if ("binary" == typeOfPhenotype) {
			  tabularUnrelateds << setw(13) << (unrelatedCount[1][liIdx] / laTotalUnrelatedCount[1])
					      << setw(13) << (unrelatedCount[0][liIdx] / laTotalUnrelatedCount[0]);
			  if (options.individual) {
			    tabularUnrelateds << setw(13) << lrChi2 << setw(5) << "1" << setw(13) << lrP << setw(13) << lrOdds[0] << setw(13) << stderror[liIdx] << setw(13) << exp(log(lrOdds[0])-1.96*stderror[liIdx]) << setw(13) << exp(log(lrOdds[0])+1.96*stderror[liIdx]);
			  } else {
			    tabularUnrelateds << setw(13) << "NA" << setw(5) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA";
			  }
			} else {
			  tabularUnrelateds << setw(13) << (unrelatedCount[0][liIdx] / laTotalUnrelatedCount[0]);
			  if (options.individual) {
			    tabularUnrelateds << setw(13) << lrChi2 << setw(5) << "1" << setw(13) << lrP << setw(13) << lrBeta << setw(13) << stderror[liIdx] << setw(13) << lrBeta-1.96*stderror[liIdx] << setw(13) << lrBeta+1.96*stderror[liIdx];
			  } else {
			    tabularUnrelateds << setw(13) << "NA" << setw(5) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA" << setw(13) << "NA";
			  }
			}
                        tabularUnrelateds << left << endl;
                    }
                }
            }
	        }
}

inline double round(double x, double y) {
    return(fabs(x) < y ? 0 : x);
}

//  output options

void UnphasedAnalysis::outputOptions(UnphasedOptions &options) {
    *outStream << endl;
    *outStream << "ANALYSIS OPTIONS" << endl;
    outStream->setf(ios::left);
    *outStream << setw(50) << "Number of nuclear families: " << familylist.size() << endl;
    *outStream << setw(50) << "Number of unrelated subjects: " << unrelateds.size() << endl;
    *outStream << "Rare haplotype maximum " << (options.cellcount ? "count" : "frequency") << setw(18) << " (-rare): " << options.rare << endl;
    *outStream << setw(50) << "Zero haplotype maximum frequency (-zero): " << options.zero << endl;
    *outStream << setw(50) << "Quantitative trait variance (-variance): " << options.variance << endl;
    *outStream << setw(50) << "Residual sibling covariance (-covariance): " << options.covariance << endl;
    *outStream << setw(50) << "Model normal distribution (-normal): " << (options.normal ? "Y" : "N") << endl;
    *outStream << setw(50) << "Uncentred quantitative traits (-uncentred): " << (options.uncentred ? "Y" : "N") << endl;
    *outStream << setw(50) << "Restrict to certain haplotypes (-certain): " << (options.certain ? "Y" : "N") << endl;
    *outStream << setw(50) << "Estimate missing genotypes (-missing): " << (options.missing ? "Y" : "N") << endl;
    *outStream << setw(50) << "Test genotypes (-genotype): " << (options.genotype ? "Y" : "N") << endl;
    *outStream << setw(50) << "Condition on genotype effects (-condgenotype): " << (options.condgenotype ? "Y" : "N") << endl;
    *outStream << setw(50) << "Assume no prior linkage (-nolinkage): " << (options.nolinkage ? "Y" : "N") << endl;
    *outStream << setw(50) << "Model odds ratio in parents (-parentrisk): " << (options.parentrisk ? "Y" : "N") << endl;
    *outStream << setw(50) << "Use only one family based control (-onefbc): " << (options.onefbc ? "Y" : "N") << endl;
    *outStream << setw(50) << "Unmatched family based controls (-hhrr): " << (options.hhrr ? "Y" : "N") << endl;
    *outStream << setw(50) << "Faster analysis of sibships (-sibship): " << (options.sibship ? "Y" : "N") << endl;
    *outStream << setw(50) << "Chromosome X (-chrX): " << (options.chrX ? "Y" : "N") << endl;
    *outStream << setw(50) << "Chromosome Y (-chrY): " << (options.chrY ? "Y" : "N") << endl;
    *outStream << setw(50) << "Convergence threshold (-epsilon): " << options.epsilon << endl;
    *outStream << setw(50) << "Numerical maximisation (-neldermead): " << (options.neldermead ? "N-M" : "D-F-P") << endl;
    *outStream << setw(50) << "Random restarts (-restarts): " << options.restarts << endl;
    *outStream << setw(50) << "Random number seed (-randomseed): " << options.randomseed << endl;
    *outStream << setw(50) << "Number of permutations (-permutation): " << options.npermutation << endl;
    *outStream << setw(50) << "Permutation distribution quantile (-quantile): " << options.quantile << endl;

    *outStream << endl;

    *outStream << "OUTPUT OPTIONS" << endl;
    *outStream << setw(50) << "Output file (-output): " << options.outputFilename << endl;
    *outStream << setw(50) << "Brief output (-brief): " << (options.briefOutput ? "Y" : "N") << endl;
    *outStream << setw(50) << "Output LD statistics (-LD): " << (options.LD ? "Y" : "N") << endl;
    *outStream << setw(50) << "Output permutation analyses (-permoutput): " << (options.outputPermutation ? "Y" : "N") << endl;
    *outStream << setw(50) << "Haplotype dump file (-dumpfile): " << options.dumpFilename << endl;
    *outStream << setw(50) << "Tabular output files (-tabularfile): " << options.tabularFilename << endl;
    *outStream << setw(50) << "Just the most likely haplotypes (-mostlikely): " << (options.mostlikely ? "Y" : "N") << endl;

    *outStream << endl;
    if (options.confounder.size() + options.modifier.size()) {
        *outStream << "COVARIATE OPTIONS" << endl;
        if (options.confounder.size()) {
            *outStream << "Confounder" << (options.confounder.size() == 1 ? ": " : "s: ");
            for (int i = 0; i < options.confounder.size(); i++) {
                *outStream << options.confounder[i] << " ";
                if (options.confounder[i] == "parsex" || options.confounder[i] == "sibsex") {
                    *outStream << "(baseline " << ((baseline[i] == MALE) ? "Male" : "Female") << ") ";
                }
            }
            *outStream << endl;
        }
        if (options.modifier.size()) {
            *outStream << "Modifier" << (options.modifier.size() == 1 ? ": " : "s: ");
            for (int i = 0; i < options.modifier.size(); i++) {
                *outStream << options.modifier[i] << " ";
                if (options.modifier[i] == "parsex" || options.modifier[i] == "sibsex") {
                    *outStream << "(baseline " << ((baseline[i+options.confounder.size()] == MALE) ? "Male" : "Female") << ") ";
                }
            }
            *outStream << endl;
        }
        if (factor.size()) {
            int nfactor = 0;
            for (int i = 0; i < factor.size(); i++) {
                nfactor += factor[i];
            }
            if (nfactor) {
                *outStream << "Factor" << (nfactor == 1 ? " " : "s ");
                *outStream << "(baseline" << (nfactor == 1 ? "): " : "s): ");
                for (int i = 0; i < factor.size(); i++) {
                    if (factor[i]) {
                        if (i < options.confounder.size()) {
                            *outStream << options.confounder[i] << " ";
                        } else {
                            *outStream << options.modifier[i-options.confounder.size()] << " ";
                        }
                        *outStream << "(" << baseline[i] << ") ";
                    }
                }
                *outStream << endl;
            }
        }
    }
    *outStream << "---------" << endl;
}



//  output results

void UnphasedAnalysis::outputResults(vector<int> &combination, string &trait,
                                     double null, double alternative, int df,
                                     UnphasedOptions &options) {

    int nhap = genoCode.size(); // number of haplotypes

    if (typeOfPhenotype == "binary") {
        *outStream << endl << "ANALYSIS OF BINARY TRAIT: " << trait << endl;
    } else {
    	if (typeOfPhenotype == "polytomous")
    		*outStream << endl << "ANALYSIS OF POLYTOMOUS TRAIT: " << trait << endl;
    	else
        *outStream << endl << "ANALYSIS OF QUANTITATIVE TRAIT: " << trait << endl;
        }
    *outStream << endl << "MARKER OPTIONS" << endl;
    *outStream << "Test marker" << (options.window == 1 ? ": " : "s: ");
    for (int i = 0; i < options.window; i++) {
        *outStream << options.marker[combination[i]] << " " ;
    }
    *outStream << endl;
    if (options.specific.size()) {
        *outStream << "Specific test " << (options.genotype ? "genotype: " : ((options.specific.size() == 1) ? "allele: " : "haplotype: "));
        Haplotype hap;
        hap = options.specific;
        *outStream << hap.str(0, options.genotype, ACGT) << endl;
    }
    if (options.condition.size()) {
        *outStream << "Conditioning marker" << (options.condition.size() == 1 ? ": " : "s: ");
        for (int i = 0; i < options.condition.size(); i++) {
            *outStream << options.condition[i] << " ";
        }
        *outStream << endl;
        if (options.condspecific.size()) {
            *outStream << "Conditioning " << (options.genotype ? "genotype: " : ((options.condition.size() == 1) ? "allele: " : "haplotype: "));
            Haplotype hap;
            hap = options.condspecific;
            *outStream << hap.str(options.condition.size()*options.condgenotype, options.genotype, ACGT) << endl;
        }
    }
    if (options.tag.size()) {
        *outStream << "Tag marker" << (options.tag.size() == 1 ? ": " : "s: ");
        for (int i = 0; i < options.tag.size(); i++) {
            *outStream << options.tag[i] << " ";
        }
        *outStream << endl;
    }
    *outStream <<  endl;

    // Haplotype data
    bool alleles = (options.window + options.condition.size() + options.tag.size() == 1);
    // width of a haplotype when printed on screen
    //  int haplotypeWidth=options.genotype?9:(alleles?7:10);
    int haplotypeWidth = 12;

    for (int i = 0; i < nhap; i++) {
        if (!zero[i]) {
            Haplotype hap = genoCode[i];
            int width = hap.str(options.condition.size() * options.condgenotype, options.genotype, ACGT).length() + 1;
            haplotypeWidth = max(haplotypeWidth, width);
        }
    }

    // sort haplotypes into alphabetical order
    vector<Haplotype> sortedHaps;
    sortedHaps = genoCode;
    sort(sortedHaps.begin(), sortedHaps.end());
    int refIndex = reference.index(genoCode);

    // output the marginal frequencies in families
    if (haveFamilies) {
        *outStream << "ESTIMATES OF MARGINAL " << (options.genotype ? "GENOTYPE" : (alleles ? "ALLELE" : "HAPLOTYPE")) << " FREQUENCIES IN FAMILIES" << endl;
        outStream->setf(ios::left);
		if (typeOfPhenotype == "polytomous")
        	{
            *outStream << setw(haplotypeWidth) << "" << setw(12) << "Level 1"
                       << setw(12) << "Level 2";
            if (K > 2) {
            *outStream << setw(12) << "Level 3";
            if (K == 4)
            *outStream << setw(12) << "Level 4";
                       }
            *outStream << setw(12) << "Level 1"
                       << setw(12) << "Level 2";
            if (K > 2) {
            *outStream << setw(12) << "Level 3";
            if (K == 4)
            *outStream << setw(12) << "Level 4";
                       }
            *outStream << endl;
            }
        *outStream << setw(haplotypeWidth) << (options.genotype ? "Genotype" : (alleles ? "Allele" : "Haplotype"));
        if (typeOfPhenotype == "quant")
            *outStream << setw(12) << "Count"
                       << setw(12) << "Freq"
                       ;
        else { if (typeOfPhenotype == "polytomous")
        	{
            *outStream << setw(12) << "Trans"
                       << setw(12) << "Trans";
            if (K > 2) {
            *outStream << setw(12) << "Trans";
            if (K == 4)
            *outStream << setw(12) << "Trans";
                       }
            *outStream << setw(12) << "Freq"
                       << setw(12) << "Freq";
            if (K > 2) {
            *outStream << setw(12) << "Freq";
            if (K == 4)
            *outStream << setw(12) << "Freq";
                       }
            }        
        else
            *outStream << setw(12) << "Trans"
                       << setw(12) << "Untrans"
                       << setw(12) << "T-Freq"
                       << setw(12) << "U-Freq"
                       ;
		}
        if (options.rare) {
            *outStream << setw(12) << "Common";
        }
        *outStream << endl;

        // calculate the frequencies for output
        double totalFamilyCount[8];
        int nCount = 2;
        if (typeOfPhenotype == "polytomous") nCount = 2*K;
        for (int i = 0; i < nCount; i++) {
            totalFamilyCount[i] = 0;
            for (int j = 0; j < nhap; j++) if (!zero[j]) {
                    totalFamilyCount[i] += familyCount[i][j];
                }
        }
        // table of haplotype data
        outStream->precision(4);
        for (int i = 0; i < sortedHaps.size(); i++) {
            int j = sortedHaps[i].index(genoCode);
            if (!zero[j]) {
                Haplotype hap = sortedHaps[i];
                *outStream << setw(haplotypeWidth) << hap.str(options.condition.size()*options.condgenotype, options.genotype, ACGT);
                if (typeOfPhenotype == "binary")
                    *outStream << setw(11) << round(familyCount[1][j], options.epsilon) << " "
                               << setw(11) << round(familyCount[0][j], options.epsilon) << " "
                               << setw(11) << round(familyCount[1][j] / totalFamilyCount[1], options.epsilon) << " "
                               << setw(11) << round(familyCount[0][j] / totalFamilyCount[0], options.epsilon) << " "
                               ;
                else { if (typeOfPhenotype == "polytomous")
                	{
                    *outStream << setw(11) << round(familyCount[0][j], options.epsilon) << " "
                               << setw(11) << round(familyCount[1][j], options.epsilon) << " ";
                    if (K > 2) {
                    *outStream << setw(11) << round(familyCount[2][j], options.epsilon) << " ";
                    if (K == 4)
                    *outStream << setw(11) << round(familyCount[3][j], options.epsilon) << " ";
                    }
                    *outStream << setw(11) << round(familyCount[0][j] / totalFamilyCount[0], options.epsilon) << " "
                               << setw(11) << round(familyCount[1][j] / totalFamilyCount[1], options.epsilon) << " ";
                    if (K > 2) {
                    *outStream << setw(11) << round(familyCount[2][j] / totalFamilyCount[2], options.epsilon) << " ";
                    if (K == 4)
                    *outStream << setw(11) << round(familyCount[3][j] / totalFamilyCount[3], options.epsilon) << " ";
					}
					}                                               
                else
                    *outStream << setw(11) << round(familyCount[0][j], options.epsilon) << " "
                               << setw(11) << round(familyCount[0][j] / totalFamilyCount[0], options.epsilon) << " "
                               ;
                }
                if (options.rare && !rare[j]) {
                    *outStream << "+";
                }
                *outStream << endl;
            }
        }
    /* Pour le phÃ©notype polytomique, afficher les transmissions aux autres niveaux */
     if (typeOfPhenotype == "polytomous") {
    		*outStream << endl;
	        *outStream << setw(haplotypeWidth) << (options.genotype ? "Genotype" : (alleles ? "Allele" : "Haplotype"));
            *outStream << setw(12) << "Untrans"
                       << setw(12) << "Untrans";
            if (K > 2) {
            *outStream << setw(12) << "Untrans";
            if (K == 4)
            *outStream << setw(12) << "Untrans";
                       }
            *outStream << setw(12) << "Freq"
                       << setw(12) << "Freq";
            if (K > 2) {
            *outStream << setw(12) << "Freq";
            if (K == 4)
            *outStream << setw(12) << "Freq";
                       }        
        if (options.rare) {
            *outStream << setw(12) << "Common";
        }
        *outStream << endl;
        for (int i = 0; i < sortedHaps.size(); i++) {
            int j = sortedHaps[i].index(genoCode);
            if (!zero[j]) {
                Haplotype hap = sortedHaps[i];
                *outStream << setw(haplotypeWidth) << hap.str(options.condition.size()*options.condgenotype, options.genotype, ACGT);
                    *outStream << setw(11) << round(familyCount[K][j], options.epsilon) << " "
                               << setw(11) << round(familyCount[K+1][j], options.epsilon) << " ";
                    if (K > 2) {
                    *outStream << setw(11) << round(familyCount[K+2][j], options.epsilon) << " ";
                    if (K == 4)
                    *outStream << setw(11) << round(familyCount[7][j], options.epsilon) << " ";
                    }
                    *outStream << setw(11) << round(familyCount[K][j] / totalFamilyCount[K], options.epsilon) << " "
                               << setw(11) << round(familyCount[K+1][j] / totalFamilyCount[K+1], options.epsilon) << " ";
                    if (K > 2) {
                    *outStream << setw(11) << round(familyCount[K+2][j] / totalFamilyCount[K+2], options.epsilon) << " ";
                    if (K == 4)
                    *outStream << setw(11) << round(familyCount[7][j] / totalFamilyCount[7], options.epsilon) << " ";
                               }                     
     		}
            if (options.rare && !rare[j]) {
                *outStream << "+";
               } 
    		*outStream << endl;       
    	}
    }
    }

    *outStream << endl;

    // output the marginal frequencies in unrelateds
    if (haveUnrelateds) {
        *outStream << "ESTIMATES OF MARGINAL " << (options.genotype ? "GENOTYPE" : (alleles ? "ALLELE" : "HAPLOTYPE")) << " FREQUENCIES IN UNRELATEDS" << endl;
        outStream->setf(ios::left);
        *outStream << setw(haplotypeWidth) << (options.genotype ? "Genotype" : (alleles ? "Allele" : "Haplotype"));
        if (typeOfPhenotype == "quant")
            *outStream << setw(12) << "Count"
                       << setw(12) << "Freq"
                       ;
        else
            *outStream << setw(12) << "Case"
                       << setw(12) << "Control"
                       << setw(12) << "Ca-Freq"
                       << setw(12) << "Co-Freq"
                       ;

        if (options.rare) {
            *outStream << setw(12) << "Common";
        }
        *outStream << endl;

        // calculate the frequencies for output
        double totalUnrelatedCount[2];
        for (int i = 0; i < 2; i++) {
            totalUnrelatedCount[i] = 0;
            for (int j = 0; j < nhap; j++) if (!zero[j]) {
                    totalUnrelatedCount[i] += unrelatedCount[i][j];
                }
        }
        // table of haplotype data
        outStream->precision(4);
        for (int i = 0; i < sortedHaps.size(); i++) {
            int j = sortedHaps[i].index(genoCode);
            if (!zero[j]) {
                Haplotype hap = sortedHaps[i];
                *outStream << setw(haplotypeWidth) << hap.str(options.condition.size()*options.condgenotype, options.genotype, ACGT);
                if (typeOfPhenotype != "quant")
                    *outStream << setw(11) << round(unrelatedCount[1][j], options.epsilon) << " "
                               << setw(11) << round(unrelatedCount[0][j], options.epsilon) << " "
                               << setw(11) << round(unrelatedCount[1][j] / totalUnrelatedCount[1], options.epsilon) << " "
                               << setw(11) << round(unrelatedCount[0][j] / totalUnrelatedCount[0], options.epsilon) << " "
                               ;
                else
                    *outStream << setw(11) << round(unrelatedCount[0][j], options.epsilon) << " "
                               << setw(11) << round(unrelatedCount[0][j] / totalUnrelatedCount[0], options.epsilon) << " "
                               ;
                if (options.rare && !rare[j]) {
                    *outStream << "+";
                }
                *outStream << endl;
            }
        }
    }

    *outStream << endl;

    // output the main likelihood ratio test
    if (options.compare1.size() && options.compare2.size()) {
        Haplotype compare1, compare2;
        compare1 = options.compare1;
        compare2 = options.compare2;
        *outStream << "TEST OF EQUALITY OF RISKS" << endl;
        *outStream << "Comparing " << compare1.str(options.condition.size()*options.condgenotype, options.genotype, ACGT) << " with " << compare2.str(options.condition.size()*options.condgenotype, options.genotype, ACGT) << endl;
    } else {
        if (options.model == "full" || options.model == "") {
            *outStream << "TEST OF OVERALL ASSOCIATION" << endl;
        }
        if (options.model == "gxg") {
            *outStream << "TEST OF GENE-GENE INTERACTION" << endl;
        }
        if (options.model == "haplomain") {
            *outStream << "TEST OF HAPLOTYPE MAIN EFFECT" << endl;
        }
        if (options.model == "allelemain") {
            *outStream << "TEST OF ALLELE MAIN EFFECTS" << endl;
        }
    }

    outStream->precision(-1);
    *outStream << "Null log-likelihood = " << null << endl;
    *outStream << "Alternative log-likelihood = " << alternative << endl;
    *outStream << "Likelihood ratio chisq = " << 2 *(alternative - null) << " df = " << df << " p-value = "
               << pochisq(2 *(alternative - null), df) << endl << endl;
    // warning if negative chisq
    if (alternative < null) {
        *outStream << "WARNING: negative chisq test, try -neldermead option" << endl;
    }

    // output the haplotype effects at the baseline levels
    *outStream << "ESTIMATES OF " << (options.genotype ? "GENOTYPE" : (alleles ? "ALLELE" : "HAPLOTYPE")) << " EFFECTS" << endl;
    if (factor.size()) {
        *outStream << "Effects are at baseline levels: ";
        for (int i = 0; i < factor.size(); i++) {
            string thisCovariate;
            if (i < options.confounder.size()) {
                thisCovariate = options.confounder[i];
            } else {
                thisCovariate = options.modifier[i-options.confounder.size()];
            }
            *outStream << thisCovariate << " (";
            if (thisCovariate == "parsex" || thisCovariate == "sibsex") {
                *outStream << (baseline[i] == MALE ? "Male" : "Female") << ") ";
            } else {
                *outStream << baseline[i] << ") ";
            }
        }
        *outStream << endl;
    }
    *outStream << "Reference " << (options.genotype ? "genotype: " : (alleles ? "allele: " : "haplotype: ")) << reference.str(options.condition.size()*options.condgenotype, options.genotype, ACGT) << endl;
    outStream->setf(ios::left);
    if (typeOfPhenotype == "polytomous") {
    *outStream << setw(haplotypeWidth) << " ";
    for (int k = 1; k < K; k++)
    	if (options.individual) *outStream << "Level " << setw(54) << k  ;
    	else *outStream << "Level " << setw(30) << k ;
    *outStream << endl;
    	}
    *outStream << setw(haplotypeWidth) << (options.genotype ? "Genotype" : (alleles ? "Allele" : "Haplotype"));
    if (typeOfPhenotype == "quant") {
        *outStream << setw(12) << "AddVal";
    } else {
        *outStream << setw(12) << "Odds-R";
    }
    *outStream << setw(12) << "95%Lo"
               << setw(12) << "95%Hi";
    if (options.individual) *outStream << setw(12) << "Chisq"
                                           << setw(12) << "P-value";
    if (typeOfPhenotype == "polytomous") {
        for (int k = 1; k < K-1; k++) {
	        *outStream << setw(12) << "Odds-R";
		    *outStream << setw(12) << "95%Lo"
               		<< setw(12) << "95%Hi";
		    if (options.individual) *outStream << setw(12) << "Chisq"
                                           << setw(12) << "P-value";
        	}
        }

    if (options.rare) {
        *outStream << setw(12) << "Common";
    }
    *outStream << endl;

    // table of haplotype data
    outStream->precision(4);
    for (int i = 0; i < sortedHaps.size(); i++) {
        int j = sortedHaps[i].index(genoCode);
        if (!zero[j]) {
            Haplotype hap = sortedHaps[i];
            *outStream << setw(haplotypeWidth) << hap.str(options.condition.size()*options.condgenotype, options.genotype, ACGT);
            double thisbeta = beta[j] - beta[refIndex];
            if (typeOfPhenotype != "quant") {
                *outStream << setw(11) << round(exp(thisbeta), options.epsilon) << " "
                           << setw(11) << round(exp(thisbeta - 1.96 * stderror[j]), options.epsilon) << " "
                           << setw(11) << round(exp(thisbeta + 1.96 * stderror[j]), options.epsilon) << " "
                           ;
            if (options.individual) {
                // this would be for the Wald test, comparing to reference
                //pvalue[i]=pochisq(thisbeta*thisbeta/thisSE,1);
                //bestpvalue=min(bestpvalue,pvalue[i]);
                // instead using stored results from a score test vs all others
                *outStream << setw(11) << chisq[j] << " "
                           << setw(11) << pvalue[j] << " ";
            }
                if (typeOfPhenotype == "polytomous") {
                    for (int k = 1; k < K-1; k++) {
                		thisbeta = beta[j + k*nhap] - beta[refIndex + k*nhap];
		                *outStream << setw(11) << round(exp(thisbeta), options.epsilon) << " "
                           		<< setw(11) << round(exp(thisbeta - 1.96 * stderror[j + k*nhap]), options.epsilon) << " "
                           		<< setw(11) << round(exp(thisbeta + 1.96 * stderror[j + k*nhap]), options.epsilon) << " "
                           		;
		            if (options.individual) {
		                *outStream << setw(11) << chisq[j + k*nhap] << " "
                           << setw(11) << pvalue[j + k*nhap] << " ";
            }
					}	
				}
			}
            else
                *outStream << setw(11) << round(thisbeta, options.epsilon) << " "
                           << setw(11) << round(thisbeta - 1.96 * stderror[j], options.epsilon) << " "
                           << setw(11) << round(thisbeta + 1.96 * stderror[j], options.epsilon) << " "
                           ;
            if (options.rare && !rare[j]) {
                *outStream << "+";
            }
            *outStream << endl;
        }
    }

    bool normal = options.normal && typeOfPhenotype == "quant";

    // confounder tests
    int nconfounder = options.confounder.size();
    int nmodifier = options.modifier.size();
    if (nconfounder) {
        *outStream << endl << "ESTIMATES OF CONFOUNDER EFFECTS" << endl;
        *outStream << "Effects are relative to those at the baseline levels" << endl;
        *outStream << "Only main effects are estimated" << endl;
    }
    for (int i = 0; i < nconfounder; i++) {
        int ix = 0;
        for (set<double>::iterator level = factorLevels[i].begin();
                level != factorLevels[i].end(); level++) {
            if (!factor[i] || *level != baseline[i]) {
                *outStream << "Confounder: " << options.confounder[i];
                if (factor[i]) {
                    *outStream << " Level: " << *level;
                }
                if (options.confounder[i] == "parsex" || options.confounder[i] == "sibsex") {
                    *outStream << " Level: " << (baseline[i] == MALE ? "Female" : "Male");
                }
                *outStream << endl;
                // table heading
                outStream->setf(ios::left);
                outStream->precision(4);
                if (typeOfPhenotype == "polytomous") {
            		*outStream << setw(haplotypeWidth) << "" << setw(24) << "Level 1";
            		if (K > 2) {
            		*outStream << setw(24) << "Level 2";
            		if (K == 4)
            		*outStream << setw(24) << "Level 3";
                    }
                *outStream << endl;
                }
                *outStream << setw(haplotypeWidth) << (options.genotype ? "Genotype" : (alleles ? "Allele" : "Haplotype"));
                *outStream << setw(12) << "Offset"
                           << setw(12) << "StdErr";
                if (typeOfPhenotype == "polytomous") {
                    for (int h = 1; h < K-1; h++) {
                	*outStream << setw(12) << "Offset"
                           << setw(12) << "StdErr";
				}
                *outStream << endl;

                double waldchisq[3] = {0,0,0};
                int walddf = 0;
                for (int j = 0; j < sortedHaps.size(); j++) {
                    int k = sortedHaps[j].index(genoCode);
                    if (!zero[k]) {
                        Haplotype hap = sortedHaps[j];
                        /* Il y aurait un betaCovariate diffŽrent pour chaque gŽnotype du locus sur lequel on conditionne, si applicable */
                        double thisfreq = betaCovariate[i][ix][k];
                        double thisSE = (sortedHaps[j] != reference) ? stderrorCovariate[i][ix][k] : 0;
                        if (sortedHaps[j] != reference) {
                            waldchisq[0] += thisfreq * thisfreq / thisSE / thisSE;
                            walddf++;
                        }
                        *outStream << setw(haplotypeWidth) << hap.str(options.condition.size()*options.condgenotype, options.genotype, ACGT);
                        *outStream << setw(11) << round(thisfreq, options.epsilon) << " "
                                   << setw(11) << round(thisSE, options.epsilon) << " ";
                    if (typeOfPhenotype == "polytomous") {
                    	for (int h = 1; h < K-1; h++) {
							thisfreq = betaCovariate[i][ix][k + h*nhap];
                         	thisSE = (sortedHaps[j] != reference) ? stderrorCovariate[i][ix][k + h*nhap] : 0;
                        	if (sortedHaps[j] != reference) {
                            	waldchisq[h] += thisfreq * thisfreq / thisSE / thisSE;
                        	}
                        	*outStream << setw(11) << round(thisfreq, options.epsilon) << " "
                                   		<< setw(11) << round(thisSE, options.epsilon) << " ";
						}
                        *outStream << endl;
                    }
                }
            }
                if (options.testconfounders) {
					if (typeOfPhenotype == "polytomous") {
					*outStream << "Omnibus Wald test of zero offsets:" << endl;
                    	for (int h = 0; h < K-1; h++) {
							double thispvalue = pochisq(waldchisq[h], walddf);
                    		*outStream << "Level " << h+1 << "chisq = " << round(waldchisq[h], options.epsilon) << " " << " df = " << walddf << " " << "p-value = " << round(thispvalue, options.epsilon) << endl;
						}
					}
				else {
				    double thispvalue = pochisq(waldchisq[0], walddf);
                    *outStream << "Omnibus Wald test of zero offsets: chisq = " << round(waldchisq[0], options.epsilon) << " " << " df = " << walddf << " " << "p-value = " << round(thispvalue, options.epsilon) << endl;
					}
                    // do we want to count these in the multiple tests??
                    //    bestpvalue=min(bestpvalue,thispvalue);
                    //    multipleTests++;
                }
                ix++;
        	}
    	}
	}
	}
    // modifier tests 
    if (nmodifier) {
        *outStream << endl << "ESTIMATES OF MODIFIER EFFECTS" << endl;
        *outStream << "Effects are relative to those at the baseline levels" << endl;
        *outStream << "Only main effects are estimated" << endl;
    }
    for (int i = nconfounder; i < nconfounder + nmodifier; i++) {
        int ix = 0;
        for (set<double>::iterator level = factorLevels[i].begin();
                level != factorLevels[i].end(); level++) {
            if (!factor[i] || *level != baseline[i]) {
                *outStream << "Modifier: " << options.modifier[i-nconfounder];
                if (factor[i]) {
                    *outStream << " Level: " << *level;
                }
                if (options.modifier[i-nconfounder] == "parsex" || options.modifier[i-nconfounder] == "sibsex") {
                    *outStream << " Level: " << (baseline[i] == MALE ? "Female" : "Male");
                }

                *outStream << endl;
                // table heading
                outStream->setf(ios::left);
                outStream->precision(4);
                if (typeOfPhenotype == "polytomous") {
            		*outStream << setw(haplotypeWidth) << "" << setw(36) << "Level 1";
            		if (K > 2) {
            		*outStream << setw(36) << "Level 2";
            		if (K == 4)
            		*outStream << setw(36) << "Level 3";
                    }
            *outStream << endl;
            }
                *outStream << setw(haplotypeWidth) << (options.genotype ? "Genotype" : (alleles ? "Allele" : "Haplotype"))
                           //         << setw(12) << "Offset"
                           //         << setw(12) << "StdErr"
                           << setw(12) << (typeOfPhenotype == "quant" ? "AddVal" : "Odds-R")
                           << setw(12) << "95%Lo"
                           << setw(12) << "95%Hi";
                if (options.testmodifiers) {
                    *outStream << setw(12) << "Chisq"
                               << setw(12) << "P-value";
                }
                if (typeOfPhenotype == "polytomous") {
        			for (int k = 1; k < K-1; k++) {
	        			*outStream << setw(12) << "Odds-R";
		    			*outStream << setw(12) << "95%Lo"
               				<< setw(12) << "95%Hi";
		    		if (options.testmodifiers) *outStream << setw(12) << "Chisq"
                                           << setw(12) << "P-value";
        			}
        		}

                *outStream << endl;
                double thischisq[3] = {0,0,0};
                for (int j = 0; j < sortedHaps.size(); j++) {
                    int k = sortedHaps[j].index(genoCode);
                    if (!zero[k] /*&& sortedHaps[j]!=reference*/) {
                        Haplotype hap = sortedHaps[j];
                        *outStream << setw(haplotypeWidth) << hap.str(options.condition.size()*options.condgenotype, options.genotype, ACGT);

                        double thisbeta = betaCovariate[i][ix][k];
                        if (!normal) {
                            thisbeta -= betaCovariate[i][ix][refIndex];
                        }
                        double thisSE = (sortedHaps[j] != reference) ? stderrorCovariate[i][ix][k] : 0;
                        thischisq[0] = (sortedHaps[j] != reference) ? thisbeta / thisSE : 0;
                        thischisq[0] *= thischisq[0];
                        if (typeOfPhenotype != "quant")
                            *outStream << setw(11) << round(exp(thisbeta), options.epsilon) << " "
                                       << setw(11) << round(exp(thisbeta - 1.96 * thisSE), options.epsilon) << " "
                                       << setw(11) << round(exp(thisbeta + 1.96 * thisSE), options.epsilon) << " ";
                        else
                            *outStream << setw(11) << round(thisbeta, options.epsilon) << " "
                                       << setw(11) << round(thisbeta - 1.96 * thisSE, options.epsilon) << " "
                                       << setw(11) << round(thisbeta + 1.96 * thisSE, options.epsilon) << " ";
                    if (typeOfPhenotype == "polytomous") {
                    	for (int h = 1; h < K-1; h++) {
							thisbeta = betaCovariate[i][ix][k + h*nhap] - betaCovariate[i][ix][refIndex + h*nhap];
                         	thisSE = (sortedHaps[j] != reference) ? stderrorCovariate[i][ix][k + h*nhap] : 0;
                        	thischisq[h] = (sortedHaps[j] != reference) ? thisbeta / thisSE : 0;
                        	thischisq[h] *= thischisq[h];
                            *outStream << setw(11) << round(exp(thisbeta), options.epsilon) << " "
                                       << setw(11) << round(exp(thisbeta - 1.96 * thisSE), options.epsilon) << " "
                                       << setw(11) << round(exp(thisbeta + 1.96 * thisSE), options.epsilon) << " ";
							}
                        }
                        if (options.testmodifiers) {
                            double thispvalue = pochisq(thischisq[0], 1);
                            *outStream << setw(11) << round(thischisq[0], options.epsilon) << " "
                                       << setw(11) << round(thispvalue, options.epsilon) << " ";
                            // do we want to count these in the multiple tests??
                            //        bestpvalue=min(bestpvalue,thispvalue);
                            //        multipleTests++;
                    		if (typeOfPhenotype == "polytomous") {
                    			for (int h = 1; h < K-1; h++) {
                            	double thispvalue = pochisq(thischisq[h], 1);
                            	*outStream << setw(11) << round(thischisq[h], options.epsilon) << " "
                                       << setw(11) << round(thispvalue, options.epsilon) << " ";
								}
							}
							}
                        *outStream << endl;
                    }
                }
                ix++;
            }
        }
    }

    for (int i = 0; i < 40; i++) {
        *outStream << "-";
    }
    *outStream << endl;

    //  Output association analysis.
    *outStream << endl;
}



//  dump haplotypes

void UnphasedAnalysis::dumpHaplotypes(UnphasedOptions &options, ofstream &dumpfile) {

    if (!options.imputeDump) {
        // families
        for (int family = 0; family < familylist.size(); family++)
            if (familylist[family].status == "" && realisationProb[family].size()) {
                int nsib = 0;
                while (consistentHaps[family][nsib] != -1) {
                    nsib++;
                }
                vector<string> sibId(nsib);
                for (int sib = 0; sib < nsib; sib++) {
                    sibId[sib] = familylist[family].sibs[consistentHaps[family][sib]].id;
                }
                dumpfile << "Family " << familylist[family].name << " ";
                dumpfile << "Sib ";
                for (int sib = 0; sib < nsib - 1; sib++) {
                    dumpfile << sibId[sib] << ",";
                }
                dumpfile << sibId[nsib-1] << " : ";

                // find the most likely configuration
                double mostLikelyProb = 0;
                for (int i = 0; i < realisationProb[family].size(); i += 4 * nsib + 1)
                    if (realisationProb[family][i+4*nsib] > mostLikelyProb) {
                        mostLikelyProb = realisationProb[family][i+4*nsib];
                    }

                double Trisk = 0, NTrisk = 0, trisk = 0, ntrisk = 0;
                for (int i = 0; i < realisationProb[family].size(); i += 4 * nsib + 1)
                    if (!options.mostlikely ||
                            realisationProb[family][i+4*nsib] == mostLikelyProb) {
                        for (int sib = 0; sib < nsib; sib++) {
                            // ft
                            Haplotype hap = genoCode[(int)realisationProb[family][i+4*sib]];
                            dumpfile << hap.str(options.condition.size()*options.condgenotype, options.genotype, ACGT);
                            trisk = beta[(int)realisationProb[family][i+4*sib]];
                            if (!options.genotype) {
                                // mt
                                hap = genoCode[(int)realisationProb[family][i+4*sib+1]];
                                dumpfile << " / " << hap.str(options.condition.size()*options.condgenotype, options.genotype, ACGT);
                                trisk += beta[(int)realisationProb[family][i+4*sib+1]];
                            }
                            // fnt
                            hap = genoCode[(int)realisationProb[family][i+4*sib+2]];
                            dumpfile << " \\ " << hap.str(options.condition.size()*options.condgenotype, options.genotype, ACGT);
                            ntrisk = beta[(int)realisationProb[family][i+4*sib+2]];
                            if (!options.genotype) {
                                // mnt
                                hap = genoCode[(int)realisationProb[family][i+4*sib+3]];
                                dumpfile << " / " << hap.str(options.condition.size()*options.condgenotype, options.genotype, ACGT);
                                ntrisk += beta[(int)realisationProb[family][i+4*sib+3]];
                            }
                            if (sib < nsib - 1) {
                                dumpfile << " , ";
                            }
                        }
                        dumpfile << " = " << realisationProb[family][i+4*nsib] << " ; ";
                        Trisk += exp(trisk);
                        NTrisk += exp(ntrisk);
                    }
                dumpfile << Trisk << " " << NTrisk << " ";
                dumpfile << endl;
            }

        int nfamily = familylist.size();

        // unrelateds
        for (int subject = 0; subject < unrelateds.size(); subject++) {
            if (unrelateds[subject].usable && realisationProb[subject].size()) {
                dumpfile << "Subject " << unrelateds[subject].id << " : ";

                // find the most likely configuration
                double mostLikelyProb = 0;
                for (int i = 0; i < realisationProb[subject+nfamily].size(); i += 3)
                    if (realisationProb[subject+nfamily][i+2] > mostLikelyProb) {
                        mostLikelyProb = realisationProb[subject+nfamily][i+2];
                    }

                // print out probabilities
                for (int i = 0; i < realisationProb[subject+nfamily].size(); i += 3)
                    if (!options.mostlikely ||
                            realisationProb[subject+nfamily][i+2] == mostLikelyProb) {
                        Haplotype hap = genoCode[(int)realisationProb[subject+nfamily][i]];
                        dumpfile << hap.str(options.condition.size()*options.condgenotype, options.genotype, ACGT);
                        if (!options.genotype) {
                            hap = genoCode[(int)realisationProb[subject+nfamily][i+1]];
                            dumpfile << " / " << hap.str(options.condition.size()*options.condgenotype, options.genotype, ACGT);
                        }
                        dumpfile << " = " << realisationProb[subject+nfamily][i+2] << " ; ";
                    }
                dumpfile << endl;
            }
        }
    }

    else {
        int nfamily = familylist.size();
        // unrelateds only
        int minAllele = 0, maxAllele = 0;
        for (int i = 0; i < genoCode.size(); i++) {
            Haplotype hap = genoCode[i];
            if (!minAllele || hap[0] <= minAllele) {
                minAllele = hap[0];
            }
            if (!maxAllele || hap[0] >= maxAllele) {
                maxAllele = hap[0];
            }
        }
        int d = maxAllele - minAllele;
        if (d == 0) {
            d = 1;
        }
        for (int subject = 0; subject < unrelateds.size(); subject++) {
            vector<double> genofreq;
            genofreq.resize(3, 0);
            // print out probabilities
            for (int i = 0; i < realisationProb[subject+nfamily].size(); i += 3) {
                Haplotype hap1 = genoCode[(int)realisationProb[subject+nfamily][i]];
                Haplotype hap2 = genoCode[(int)realisationProb[subject+nfamily][i+1]];
                int ix = (hap1[0] + hap2[0] - 2 * minAllele) / d;
                genofreq[ix] += realisationProb[subject+nfamily][i+2];
            }
            dumpfile << unrelateds[subject].id << " " << genofreq[0] << " " << genofreq[1] << " " << genofreq[2] << " ";
        }
        dumpfile << endl;


    }

}



//  outputLD

void UnphasedAnalysis::outputLD(UnphasedOptions &options) {

    int nhap = genoCode.size();
    // D' and rsq
    double globalDprime[2][2]; // families, unrelateds, cases, controls
    valarray<double> Dprime[2][2], rsq[2][2];
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) {
            Dprime[i][j].resize(nhap, 0);
            rsq[i][j].resize(nhap, 0);
        }

    for (int ix = 0; ix < 2; ix++)
        for (int iy = 0; iy < 2; iy++) {

            globalDprime[ix][iy] = 0;
            Haplotype hap;

            valarray<double> trueFreq(genoCode.size());
            double totalCount = 0;
            for (int j = 0; j < genoCode.size(); j++) if (!zero[j]) {
                    totalCount += ix ? unrelatedCount[iy][j] : familyCount[iy][j];
                }
            for (int j = 0; j < trueFreq.size(); j++) if (!zero[j]) {
                    trueFreq[j] = (ix ? unrelatedCount[iy][j] : familyCount[iy][j]) / totalCount;
                }

            // pick up each allele or genotype at first marker
            valarray<bool> counted1(nhap);
            counted1 = false;
            for (int i = 0; i < nhap; i++)
                if (!counted1[i]) {
                    hap = genoCode[i];
                    bool genotype1 = options.condition.size() && options.condgenotype ||
                                     options.genotype;
                    double freq0 = 0;
                    for (int j = i; j < nhap; j++)
                        if (genoCode[j][0] == genoCode[i][0] &&
                                (!genotype1 || genoCode[j][1] == genoCode[i][1])) {
                            if (!zero[j]) {
                                freq0 += trueFreq[j];
                            }
                            counted1[j] = true;
                        }
                    // do the same for the the remaining markers
                    valarray<bool> counted2(nhap);
                    counted2 = false;
                    for (int k = 0; k < nhap; k++)
                        if (!counted2[k]) {
                            for (int l = 1 + genotype1; l < genoCode[k].size(); l++) {
                                hap[l] = genoCode[k][l];
                            }
                            double freq1 = 0;
                            for (int l = k; l < nhap; l++) {
                                bool match = true;
                                for (int m = 1 + genotype1; m < genoCode[k].size(); m++)
                                    if (genoCode[l][m] != genoCode[k][m]) {
                                        match = false;
                                    }
                                if (match) {
                                    if (!zero[l]) {
                                        freq1 += trueFreq[l];
                                    }
                                    counted2[l] = true;
                                }
                            }
                            double D = 0;
                            if (hap.index(genoCode) < nhap && !zero[hap.index(genoCode)]) {
                                D = trueFreq[hap.index(genoCode)];
                            }
                            D -= freq0 * freq1;
                            if (freq0 * freq1 && hap.index(genoCode) < nhap) {
                                rsq[ix][iy][hap.index(genoCode)] = D * D / (freq0 * freq1 * (1 - freq0) * (1 - freq1));
                            }
                            double Dmax;
                            if (D < 0)
                                Dmax = freq0 * freq1 < (1 - freq0) * (1 - freq1) ?
                                       freq0 * freq1 : (1 - freq0) * (1 - freq1);
                            else
                                Dmax = (1 - freq0) * freq1 < freq0 * (1 - freq1) ?
                                       (1 - freq0) * freq1 : freq0 * (1 - freq1);
                            if (hap.index(genoCode) < nhap) {
                                Dprime[ix][iy][hap.index(genoCode)] = Dmax ? D / Dmax : 0;
                            }
                            globalDprime[ix][iy] += fabs(Dprime[ix][iy][hap.index(genoCode)]) * freq0 * freq1;
                        }
                }
        }

    // sort haplotypes into alphabetical order
    vector<Haplotype> sortedHaps;
    sortedHaps = genoCode;
    sort(sortedHaps.begin(), sortedHaps.end());

    // output the LD measures
    for (int ix = 0; ix < 2; ix++) {
        if (ix == 0 && haveFamilies || ix == 1 && haveUnrelateds) {
            if (ix == 0 && haveFamilies) {
                *outStream << "ESTIMATES OF LINKAGE DISEQUILIBRIUM IN FAMILIES" << endl;
            }
            if (ix == 1 && haveUnrelateds) {
                *outStream << "ESTIMATES OF LINKAGE DISEQUILIBRIUM IN UNRELATEDS" << endl;
            }

            // width of a haplotype when printed on screen
            int haplotypeWidth = options.genotype ? 9 : 10;

            for (int i = 0; i < sortedHaps.size(); i++) {
                if (!zero[i]) {
                    Haplotype hap = genoCode[i];
                    int width = hap.str(options.condition.size() * options.condgenotype, options.genotype, ACGT).length() + 1;
                    haplotypeWidth = max(haplotypeWidth, width);
                }
            }
            outStream->setf(ios::left);
            outStream->precision(4);
            *outStream << setw(16) << (options.genotype ? "Genotype" : "Haplotype");
            if (typeOfPhenotype == "quant") {
                *outStream << setw(12) << "D'"
                           << setw(12) << "r2"
                           << endl;
            } else {
                if (ix == 0 && haveFamilies)
                    *outStream << setw(12) << "T-D'"
                               << setw(12) << "T-r2"
                               << setw(12) << "U-D'"
                               << setw(12) << "U-r2"
                               << endl;
                if (ix == 1 && haveUnrelateds)
                    *outStream << setw(12) << "ca-D'"
                               << setw(12) << "ca-r2"
                               << setw(12) << "co-D'"
                               << setw(12) << "co-r2"
                               << endl;
            }
            int ncondition = options.condgenotype * options.condition.size();
            for (int i = 0; i < nhap; i++) {
                int j = sortedHaps[i].index(genoCode);
                if (!zero[j]) {
                    Haplotype hap = genoCode[j];
                    *outStream << setw(15) << hap.str(ncondition, options.genotype, ACGT) << " ";
                    if (typeOfPhenotype == "quant") {
                        *outStream << setw(11) << Dprime[ix][0][j] << " "
                                   << setw(11) << rsq[ix][0][j] << endl;
                    } else {
                        *outStream << setw(11) << Dprime[ix][1][j] << " "
                                   << setw(11) << rsq[ix][1][j] << " "
                                   << setw(11) << Dprime[ix][0][j] << " "
                                   << setw(11) << rsq[ix][0][j] << endl;
                    }
                }
            }
            *outStream << setw(16) << "Global D'";
            if (typeOfPhenotype == "quant") {
                *outStream << setw(11) << globalDprime[ix][0] << endl;
            } else {
                *outStream << setw(11) << globalDprime[ix][1] << "             " << setw(11) << globalDprime[ix][0] << endl;
            }
        }
    }

    for (int i = 0; i < 40; i++) {
        *outStream << "-";
    }
    *outStream << endl;
}



