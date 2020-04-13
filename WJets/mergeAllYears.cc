#include <iostream>
#include <sstream>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TString.h>

#include "SourceFiles/fileNames.h"

// Grab histograms for all three years from separate HistoFiles folders
// Simply add the histograms for all three years of each BG/signal/data and write out in new file that goes into HistoFiles_Run2
// Just doing simple addition of histograms here
// Only need to do signal region (no QCD control regions or systematics)

void mergeSampleAllYears(TString lepSelection = "SMu", int systematics = 0, int jetPtCutMin = 30, int doQCD = 0, int METcut = 0, int doBJets = -1, TString sampleName = "");

void mergeAllYears(){

    // run script separately for each data/signal/BG file
    int nFiles = NFILESTTBARWJETS;

    for (unsigned short i = 0; i < nFiles; i++){
        int fileSelect = FilesTTbarWJets[i];
        TString fileNameTemp =  ProcessInfo[fileSelect].filename;

        mergeSampleAllYears("SMu", 0, 30, 0, 0, 0, fileNameTemp); // no bveto, doQCD=0
    }

}

void mergeSampleAllYears(TString lepSelection, int systematics, int jetPtCutMin, int doQCD, int METcut, int doBJets, TString sampleName){

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    cout << "\n >>>>> Running mergeSampleAllYears for " << sampleName << " <<<<< " << endl;

    //--- Read in input arguments
    ostringstream strJetPtCutMin; strJetPtCutMin << jetPtCutMin;

    ostringstream doQCDStr;       doQCDStr << doQCD;
    TString QCDStr;
    if (doQCD > 0) QCDStr = "_QCD"+doQCDStr.str();
    else QCDStr = "";

    ostringstream METcutStream; METcutStream << METcut;
    TString METcutStr;
    if (METcut > 0) METcutStr = "_MET"+METcutStream.str();
    else METcutStr = "";

    ostringstream doBVetoStream; doBVetoStream << doBJets;
    TString doBVetoStr;
    if (doBJets == 0) doBVetoStr = "";  
    else if (doBJets < 0) doBVetoStr = "_BVeto";
    else doBVetoStr = "_BJets";

    TString syst;
    if (systematics == 0) syst = "Syst_0_";
    else if (systematics == 1) syst = "Syst_1_Up_"; 
    else if (systematics == -1) syst = "Syst_1_Down_"; 
    else if (systematics == 3) syst = "Syst_3_Up_"; 
    else if (systematics == -3) syst = "Syst_3_Down_"; 
    else if (systematics == 5) syst = "Syst_5_Up_"; 
    else if (systematics == -5) syst = "Syst_5_Down_"; 
    else if (systematics == 6) syst = "Syst_6_Up_"; 
    else if (systematics == -6) syst = "Syst_6_Down_"; 

    //--- Whether muon eff. corrections are applied or not
    TString effCorrString;
    if (sampleName.Contains("Data")) effCorrString = "_EffiCorr_0";
    else effCorrString = "_EffiCorr_1";
  
    //--- Print out input arguments
    cout << "lepSelection = " << lepSelection << endl;
    cout << "syst = " << syst << endl;
    cout << "jetPtCutMin = " << jetPtCutMin << endl;
    cout << "doQCD = " << doQCD << endl;
    cout << "METcut = " << METcut << endl;
    cout << "doBJets = " << doBJets << endl;

    //--- Assemble input filenames
    TString str1, str2, str3, strf;

    str1 = "HistoFiles_2016/" + lepSelection +  "_13TeV_" + sampleName + effCorrString + "_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth" + doBVetoStr + QCDStr + METcutStr + ".root";
    str2 = "HistoFiles_2017/" + lepSelection +  "_13TeV_" + sampleName + effCorrString + "_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth" + doBVetoStr + QCDStr + METcutStr + ".root";
    str3 = "HistoFiles_2018/" + lepSelection +  "_13TeV_" + sampleName + effCorrString + "_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth" + doBVetoStr + QCDStr + METcutStr + ".root";
    strf = "HistoFiles_Run2/" + lepSelection +  "_13TeV_" + sampleName + effCorrString + "_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth" + doBVetoStr + QCDStr + METcutStr + ".root";

    cout << "\nInput files:\n" << str1 << "\n" << str2 << "\n" << str3 << std::endl;
    cout << "Output file: \n" << strf << "\n" << endl;

    //--- Create new output directory for merged files
    TString command = "mkdir -p HistoFiles_Run2";
    system(command);

    //--- Open files
    TFile *f1 = new TFile(str1, "read");
    TFile *f2 = new TFile(str2, "read");
    TFile *f3 = new TFile(str3, "read");
    TFile *ff = new TFile(strf, "RECREATE");

    //--- Loop over every histo and add all individual contributions
    int nHist = f1->GetListOfKeys()->GetEntries();
    // std::cout << "nHist = " << nHist << std::endl;
    for (int i(0); i < nHist; i++){
        string hName = f1->GetListOfKeys()->At(i)->GetName(); //the "name" should be set as the same as the pointer in the code

        // NOTE: this script is only set up to handle TH1D's at this time (hresponse histograms are TH2D)
        if (hName.find("hresponse") != string::npos) continue;
        else{
            // std::cout << "Doing histogram " << i << ": " << hName << std::endl;
            TH1D *h1 = (TH1D*) f1->Get(hName.c_str()); 
            TH1D *h2 = (TH1D*) f2->Get(hName.c_str()); 
            TH1D *h3 = (TH1D*) f3->Get(hName.c_str()); 
            TH1D *hSum = (TH1D*) h1->Clone();
            hSum->Add(h2);
            hSum->Add(h3);
            //--- Write out summed histo to outfile
            ff->cd();
            hSum->Write();
        }

    } //end loop over all histograms

    cout << "\nClosing files!" << endl;
    f1->Close();
    f2->Close();
    f3->Close();
    ff->Close();
}