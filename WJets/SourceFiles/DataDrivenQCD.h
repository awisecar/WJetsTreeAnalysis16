#ifndef _DataDrivenQCD_H_
#define _DataDrivenQCD_H_

#include <iostream>

using namespace std;

void DataDrivenQCD( string leptonFlavor = "SMu", int METcut = 0,  int doBJets = -1);
void FuncOpenAllFiles(TFile *fData[], TFile *fMC[][11], string leptonFlavor = "SMu",  int METcut = 0, bool doFlat = false, bool doVarWidth = true , int doBJets = -1 );
vector<string> getVectorOfHistoNames(TFile *fData[]);
void FuncDataDrivenQCD(string variable, TFile *fData[], TFile *fMC[][11], TFile*);

#endif

