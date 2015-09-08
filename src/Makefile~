# Makefile for POLYUNPHASED
#
#   Copyright (c) 2006 Frank Dudbridge
#   MRC Biostatistics Unit
#  Robinson Way
#  Cambridge CB2 0SR, UK
#  frank.dudbridge@mrc-bsu.cam.ac.uk

#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
#   USA.


### You may need to give the full path to g++ here
CCC=g++

### You may need to give the full path to javac here
JAVAC=javac

### You may need to give the full path to jar here
JAR=jar

#SYNTAX=-fsyntax-only

OPTIONS =  -Wno-deprecated -O3

polyunphased: unphased.cpp UnphasedOptions.o LinkageData.o NuclearFamily.o Haplotype.o UnphasedAnalysis.o UnphasedAnalysisModels.o UnphasedAnalysisExploratory.o UnphasedAnalysisLikelihood.o UnphasedAnalysisOutput.o matinv.o pvalues.o asran.o UnphasedScore.o UnphasedScoreFamily.o UnphasedScoreUnrelated.o 
	$(CCC) -o polyunphased unphased.cpp UnphasedOptions.o LinkageData.o NuclearFamily.o Haplotype.o UnphasedAnalysis.o UnphasedAnalysisModels.o UnphasedAnalysisExploratory.o UnphasedAnalysisLikelihood.o UnphasedAnalysisOutput.o UnphasedScore.o UnphasedScoreFamily.o UnphasedScoreUnrelated.o pvalues.o asran.o matinv.o $(SYNTAX) $(OPTIONS)

UnphasedOptions.o: UnphasedOptions.cpp UnphasedOptions.h Options.h
	$(CCC) -c UnphasedOptions.cpp $(SYNTAX) $(OPTIONS)

LinkageData.o: LinkageData.cpp LinkageData.h getline.h
	$(CCC) -c LinkageData.cpp $(SYNTAX) $(OPTIONS)

NuclearFamily.o: NuclearFamily.cpp NuclearFamily.h LinkageData.h
	$(CCC) -c NuclearFamily.cpp $(SYNTAX) $(OPTIONS)

Haplotype.o: Haplotype.cpp Haplotype.h 
	$(CCC) -c Haplotype.cpp $(SYNTAX) $(OPTIONS)

UnphasedAnalysis.o: UnphasedAnalysis.cpp UnphasedAnalysis.h LinkageData.h UnphasedOptions.h NuclearFamily.h matinv.h pvalues.h asran.h
	$(CCC) -c UnphasedAnalysis.cpp $(SYNTAX) $(OPTIONS)

UnphasedAnalysisModels.o: UnphasedAnalysisModels.cpp UnphasedAnalysis.h LinkageData.h UnphasedOptions.h NuclearFamily.h matinv.h pvalues.h
	$(CCC) -c UnphasedAnalysisModels.cpp $(SYNTAX) $(OPTIONS)

UnphasedAnalysisExploratory.o: UnphasedAnalysisExploratory.cpp UnphasedAnalysis.h LinkageData.h UnphasedOptions.h NuclearFamily.h
	$(CCC) -c UnphasedAnalysisExploratory.cpp $(SYNTAX) $(OPTIONS)

UnphasedAnalysisLikelihood.o: UnphasedAnalysisLikelihood.cpp UnphasedAnalysis.h LinkageData.h UnphasedOptions.h NuclearFamily.h asran.h
	$(CCC) -c UnphasedAnalysisLikelihood.cpp $(SYNTAX) $(OPTIONS)

UnphasedAnalysisOutput.o: UnphasedAnalysisOutput.cpp UnphasedAnalysis.h LinkageData.h UnphasedOptions.h NuclearFamily.h pvalues.h
	$(CCC) -c UnphasedAnalysisOutput.cpp $(SYNTAX) $(OPTIONS)

UnphasedScore.o: UnphasedScore.cpp UnphasedAnalysis.h LinkageData.h UnphasedOptions.h NuclearFamily.h matinv.h pvalues.h asran.h
	$(CCC) -c UnphasedScore.cpp $(SYNTAX) $(OPTIONS)

UnphasedScoreFamily.o: UnphasedScoreFamily.cpp UnphasedAnalysis.h LinkageData.h UnphasedOptions.h NuclearFamily.h matinv.h pvalues.h asran.h
	$(CCC) -c UnphasedScoreFamily.cpp $(SYNTAX) $(OPTIONS)

UnphasedScoreUnrelated.o: UnphasedScoreUnrelated.cpp UnphasedAnalysis.h LinkageData.h UnphasedOptions.h NuclearFamily.h matinv.h pvalues.h asran.h
	$(CCC) -c UnphasedScoreUnrelated.cpp $(SYNTAX) $(OPTIONS)

asran.o: asran.cpp asran.h 
	$(CCC) -c asran.cpp $(SYNTAX) $(OPTIONS)

pvalues.o: pvalues.cpp pvalues.h 
	$(CCC) -c pvalues.cpp $(SYNTAX) $(OPTIONS)

matinv.o: matinv.cpp matinv.h
	$(CCC) -c matinv.cpp $(SYNTAX) $(OPTIONS)

unphased.jar: Unphased.java UnphasedManifest
	$(JAVAC) Unphased.java
	$(JAR) cmf UnphasedManifest unphased.jar *.class

clean:
	rm -f *.o *.class

cleanall:
	rm -f *.o *.class unphased unphased.jar
