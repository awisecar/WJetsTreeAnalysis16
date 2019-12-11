#include <iostream>
#include <sstream>
#include <vector>
#include <TFile.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TLegend.h>

#include "SourceFiles/fileNames.h"
#include "SourceFiles/getFilesAndHistograms.h"

using namespace std;

void runCalcBTagEffsMC(string leptonFlavor = "SMu", int year = 2017, int JetPtMin = 30,
    int doQCD = 0, bool doSSign = 0, bool doInvMassCut = 0, int METcut = 0 , int doBJets = -1, 
    int JetPtMax = 0, bool doFlat = 0, bool doVarWidth = 1);

void calcBTagEffsMC(){
    welcomeMessage();

    // 2016 --
    runCalcBTagEffsMC("SMu", 2016, 30, 0, 0, 0, 0, 0, 0, 0, 1); //no MET cut, no btag veto

 //   // 2017 --
 //   // runCalcBTagEffsMC("SMu", 2017, 30, 0, 0, 0, 0, -1, 0, 0, 1); //no MET cut
//    runCalcBTagEffsMC("SMu", 2017, 30, 0, 0, 0, 0, 0, 0, 0, 1); //no MET cut, no btag veto
}

void runCalcBTagEffsMC(string leptonFlavor, int year, int JetPtMin,
    int doQCD, bool doSSign, bool doInvMassCut, int METcut, int doBJets, 
    int JetPtMax, bool doFlat, bool doVarWidth)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    // prevents the addition of histograms in memory
    TH1::AddDirectory(kFALSE); //sets a global switch disabling the reference

    string energy = "13TeV";

    std::string yearStr;
    std::stringstream yearSStr;
    yearSStr << year;
    yearStr = yearSStr.str();

    // assuming W+Jets, leptonFlavor == "SMu" here 
    int nFiles = NFILESTTBARWJETS_NOQCD_NODATA;
    // int nFiles = 1; // just running ttBar

    TFile *file[nFiles];   
    std::vector <std::string> fileNames;
    int openedFiles(0);
    cout << "\n >>>>> Getting all files..." << endl;
    for (int i = 0; i < nFiles; i++){

        int fileSelect = FilesTTbarWJets_NoQCD_NoData[i];
        // int fileSelect = 9; // just running ttBar

        string fileNameTemp =  ProcessInfo[fileSelect].filename;
        cout << "#" << i << " -----> fileNameTemp = " << fileNameTemp << endl; 
        file[i] = getFile(FILESDIRECTORY, leptonFlavor, energy, fileNameTemp, JetPtMin, JetPtMax, doFlat, doVarWidth, doQCD, doSSign, doInvMassCut, METcut, doBJets);
        fileNames.push_back(fileNameTemp);

        // check if file is opened correctly
        if (file[i]->IsOpen()) openedFiles++;
        cout << endl;
    }
    if (openedFiles == nFiles) cout << "\n >>>>> All files found!" << endl;
    else {
        cout << "\n >>>>> Missing file(s)! Aborting!\n" << endl;
        return; //should catch error and exit the main function
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    
    TH1D *hBJetTagged = NULL;
    TH1D *hBJet = NULL;
    TH1D *hCJetTagged = NULL;
    TH1D *hCJet = NULL;
    TH1D *hLightJetTagged = NULL;
    TH1D *hLightJet = NULL;
     
    cout << "\n >>>>> Getting efficiency histos from files..." << endl;
    // Loop through all files
    for (int i = 0; i < nFiles; i++) {

        cout << "\nFile #" << i << ": " << fileNames.at(i) << endl;
        int nHist(0); 
        nHist = file[i]->GetListOfKeys()->GetEntries(); 
        
        // Loop through all histograms in each file to find the eff histos
        for (int j(0); j < nHist; j++) {

            string histoNameTemp = "";
            histoNameTemp = file[i]->GetListOfKeys()->At(j)->GetName();
            if (histoNameTemp.find("gen") != string::npos) continue; // only reco-level histos

            TH1D *histTemp = NULL;
            histTemp = (TH1D*) file[i]->Get(histoNameTemp.c_str());
            // histTemp->SetDirectory(0); // not necessary because of TH1::AddDirectory(kFALSE) statement above
            if (histTemp->GetEntries() < 1) continue; // skip empty histos
            if (!histTemp->InheritsFrom(TH1D::Class())) continue; // skip TH2's

            // get the b-tagging efficiency histos
            // structure if statements to search for the more specific string first
            // we assume that the index of the histograms is the same in each file (set in HistoSet.h)

            if (histoNameTemp.find("h_pt_b_tagged") != string::npos){
                cout << " >>> histo #" << j << ": " << histoNameTemp << std::endl;
                if (i == 0) hBJetTagged = (TH1D*) histTemp->Clone();
                else hBJetTagged->Add(histTemp);

                int nBins(histTemp->GetNbinsX());
                for (int k(1); k <= nBins; k++) {
                    cout << histTemp->GetBinContent(k) << "  ";
                }
                cout << "\n" << endl;
            }
            else if (histoNameTemp.find("h_pt_b") != string::npos){
                cout << " >>> histo #" << j << ": " << histoNameTemp << std::endl;
                if (i == 0) hBJet = (TH1D*) histTemp->Clone();
                else hBJet->Add(histTemp);

                int nBins(histTemp->GetNbinsX());
                for (int k(1); k <= nBins; k++) {
                    cout << histTemp->GetBinContent(k) << "  ";
                }
                cout << "\n" << endl;
            }
            else if (histoNameTemp.find("h_pt_c_tagged") != string::npos){
                cout << " >>> histo #" << j << ": " << histoNameTemp << std::endl;
                if (i == 0) hCJetTagged = (TH1D*) histTemp->Clone();
                else hCJetTagged->Add(histTemp);

                int nBins(histTemp->GetNbinsX());
                for (int k(1); k <= nBins; k++) {
                    cout << histTemp->GetBinContent(k) << "  ";
                }
                cout << "\n" << endl;
            }
            else if (histoNameTemp.find("h_pt_c") != string::npos){
                cout << " >>> histo #" << j << ": " << histoNameTemp << std::endl;
                if (i == 0) hCJet = (TH1D*) histTemp->Clone();
                else hCJet->Add(histTemp);

                int nBins(histTemp->GetNbinsX());
                for (int k(1); k <= nBins; k++) {
                    cout << histTemp->GetBinContent(k) << "  ";
                }
                cout << "\n" << endl;
            }
            else if (histoNameTemp.find("h_pt_udsg_tagged") != string::npos){
                cout << " >>> histo #" << j << ": " << histoNameTemp << std::endl;
                if (i == 0) hLightJetTagged = (TH1D*) histTemp->Clone();
                else hLightJetTagged->Add(histTemp);

                int nBins(histTemp->GetNbinsX());
                for (int k(1); k <= nBins; k++) {
                    cout << histTemp->GetBinContent(k) << "  ";
                }
                cout << "\n" << endl;
            }
            else if (histoNameTemp.find("h_pt_udsg") != string::npos){
                cout << " >>> histo #" << j << ": " << histoNameTemp << std::endl;
                if (i == 0) hLightJet = (TH1D*) histTemp->Clone();
                else hLightJet->Add(histTemp);

                int nBins(histTemp->GetNbinsX());
                for (int k(1); k <= nBins; k++) {
                    cout << histTemp->GetBinContent(k) << "  ";
                }
                cout << "\n" << endl;
            }
            // else cout << "Found histogram " << histoNameTemp << " at index #" << j << std::endl;
            else continue;

        }
    } 

    ///////////////////////////////////////////////////////////////////////////////////////

    cout << "\n///////////////////////////////////////////////////////////////////////////////////////" << endl;
    cout << ">>>>> Sum totals:\n" << endl;
    int nBins(hBJetTagged->GetNbinsX());

    // BJets -----
    cout << "hBJetTagged" << endl;
    for (int i(1); i <= nBins; i++){
        cout << hBJetTagged->GetBinContent(i) << "  ";
    }
    cout << "\n" << endl;
    cout << "hBJet" << endl;
    for (int i(1); i <= nBins; i++){
        cout << hBJet->GetBinContent(i) << "  ";
    }
    cout << "\n" << endl;

    // Charm Jets -----
    cout << "hCJetTagged" << endl;
    for (int i(1); i <= nBins; i++){
        cout << hCJetTagged->GetBinContent(i) << "  ";
    }
    cout << "\n" << endl;
    cout << "hCJet" << endl;
    for (int i(1); i <= nBins; i++){
        cout << hCJet->GetBinContent(i) << "  ";
    }
    cout << "\n" << endl;

    // Light Jets -----
    cout << "hLightJetTagged" << endl;
    for (int i(1); i <= nBins; i++){
        cout << hLightJetTagged->GetBinContent(i) << "  ";
    }
    cout << "\n" << endl;
    cout << "hLightJet" << endl;
    for (int i(1); i <= nBins; i++){
        cout << hLightJet->GetBinContent(i) << "  ";
    }
    cout << "\n" << endl;

    // Calculate the efficiencies
    cout << "///////////////////////////////////////////////////////////////////////////////////////" << endl;
    cout << ">>>>> Calculating efficiencies:\n" << endl;
    TH1D* hBJetEff     = (TH1D*) hBJetTagged->Clone("h_pt_bjet_eff");
    TH1D* hCJetEff     = (TH1D*) hCJetTagged->Clone("h_pt_cjet_eff");
    TH1D* hLightJetEff = (TH1D*) hLightJetTagged->Clone("h_pt_udsgjet_eff");
    // binomial errors for efficiencies
    hBJetEff->Divide(hBJetTagged, hBJet, 1, 1, "B");
    hCJetEff->Divide(hCJetTagged, hCJet, 1, 1, "B");
    hLightJetEff->Divide(hLightJetTagged, hLightJet, 1, 1, "B");

    // BJets -----
    cout << "hBJetEff" << endl;
    for (int i(1); i <= nBins; i++){
        // cout << hBJetEff->GetBinContent(i) << "   ";
        cout << hBJetEff->GetBinContent(i) << std::endl;
    }
    cout << "\n" << endl;

    // Charm Jets -----
    cout << "hCJetEff" << endl;
    for (int i(1); i <= nBins; i++){
        // cout << hCJetEff->GetBinContent(i) << "   ";
        cout << hCJetEff->GetBinContent(i) << std::endl;
    }
    cout << "\n" << endl;

    // Light Jets -----
    cout << "hLightJetEff" << endl;
    for (int i(1); i <= nBins; i++){
        // cout << hLightJetEff->GetBinContent(i) << "   ";
        cout << hLightJetEff->GetBinContent(i) << std::endl;
    }
    cout << endl;

    ///////////////////////////////////////////////////////////////////////////////////////
    TString outputDirBtag = "EfficiencyTables/bTagEffSFs/";
    // TString command = "mkdir -p "+outputDirBtag;
    // system(command);

    cout << "\n///////////////////////////////////////////////////////////////////////////////////////" << endl;
    TString outputFilename = outputDirBtag+"bTagEffsMC_DeepCSV_"+yearStr+".root";
    // TString outputFilename = "bTagEffsMC_"+yearStr+".root";
    cout << ">>>>> Writing histograms to output file at --- \n" << outputFilename << endl;

    TFile *outputFile = new TFile(outputFilename, "RECREATE");
    outputFile->cd();

    hBJet->Write();
    hBJetTagged->Write();
    hBJetEff->Write();
    hCJet->Write();
    hCJetTagged->Write();
    hCJetEff->Write();
    hLightJet->Write();
    hLightJetTagged->Write();
    hLightJetEff->Write();

    outputFile->Close();
    delete outputFile;

    ///////////////////////////////////////////////////////////////////////////////////////

    cout << "\n>>>>> Closing all files..." << endl;
    for (int i(0); i < nFiles; i++) closeFile(file[i]);

    cout << "\nExiting runCalcBTagEffsMC!" << endl;
}
