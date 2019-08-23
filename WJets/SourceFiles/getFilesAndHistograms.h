#ifndef _getFilesAndhistograms_H_
#define _getFilesAndhistograms_H_

#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <string>
#include <vector>
using namespace std;

string getEnergy();
TFile* getFile(string, string, string, string, int JetPtMin = 30, int JetPtMax = 0, bool doFlat = false, bool doVarWidth = false, int doQCD = 0 ,  bool doSSign = false , bool doInvMassCut = false , int MET = 0 , int doBJets = 0 , string closureTest = "", string syst = "0", bool dodR = false,bool useUnfoldingFiles = false );
void getFiles(string, TFile *Files[], string, string, string, int JetPtMin = 20, int JetPtMax = 0, bool doFlat = false ,  bool doVarWidth = false , int doQCD = 0,  bool doSSign = false , bool doInvMassCut = false, int MET = 0  , int doBJets = 0 , bool useUnfoldingFiles = false );
void closeFile(TFile*);
void closeFiles(TFile *Files[]);
void closeFiles(TFile *Files[], int nFiles);
TH1D* getHisto(TFile*, string);
void getHistos(TH1D *histograms[], TFile *Files[], string, bool isDoubleLep = true );
void getStatistics(string leptonFlavor = "Muons",  int JetPtMin = 30, int JetPtMax = 0,  bool doFlat = false , bool doVarWidth = true, int doQCD = 0 , bool doSSign = false ,  bool doInvMassCut = false , int METcut = 0 , int doBjets = 0 , bool doTTScale = false );
#endif

