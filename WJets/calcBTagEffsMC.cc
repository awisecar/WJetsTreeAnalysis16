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
     runCalcBTagEffsMC("SMu", 2017, 30, 0, 0, 0, 0, -1, 0, 0, 1); //no MET cut
}

void runCalcBTagEffsMC(string leptonFlavor, int year, int JetPtMin,
    int doQCD, bool doSSign, bool doInvMassCut, int METcut, int doBJets, 
    int JetPtMax, bool doFlat, bool doVarWidth)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    string energy = "13TeV";

    // assuming W+Jets, leptonFlavor == "SMu" here 
    int nFiles = NFILESTTBARWJETS_NOQCD;
    TFile *file[nFiles];    
    int countFiles = 0;
    cout << "\n >>>>> Getting all files..." << endl;
    for (unsigned short i = 0; i < nFiles; i++){

        int fileSelect = FilesTTbarWJets_NoQCD[i];
        string fileNameTemp =  ProcessInfo[fileSelect].filename;
        cout << " --> fileNameTemp = " << fileNameTemp << endl; 
        file[countFiles] = getFile(FILESDIRECTORY, leptonFlavor, energy, fileNameTemp, JetPtMin, JetPtMax, doFlat, doVarWidth, doQCD, doSSign, doInvMassCut, METcut, doBJets);

        // check if file is opened correctly
        if (file[countFiles]->IsOpen()) countFiles++;
    }
    // cout << "nFiles = " << nFiles << ", countFiles = " << countFiles << endl;
    if (countFiles == nFiles) cout << "\n >>>>> All files found!" << endl;
    else {
        cout << "\n >>>>> Missing file(s)!!!!!" << endl;
        return; //should catch error and exit the main function
    }

    unsigned short nHist = file[0]->GetListOfKeys()->GetEntries();
    string histoNameTemp;
    TH1D *histTemp;
    vector<string> histoName;

    int nHistNoGen = 0;
    cout << "\n >>>>> Looking at histograms in files..." << endl;
    for (unsigned short i(0); i < nHist; i++) {

        histoNameTemp = file[0]->GetListOfKeys()->At(i)->GetName();
        if (histoNameTemp.find("gen") != string::npos) continue; // is it not a reco-level histo?

        histTemp = (TH1D*) file[0]->Get(histoNameTemp.c_str());
        if (histTemp->GetEntries() < 1) continue; // is it empty?
        if (!histTemp->InheritsFrom(TH1D::Class())) continue; // is it a TH2?

        cout << "Found histogram #" << i << ": " << histoNameTemp << endl;
        histoName.push_back(histoNameTemp);

        nHistNoGen++;
    }
    cout << "Number of relevant histograms: " << nHistNoGen << endl;

    TH1D* hist[countFiles][nHistNoGen];
    cout << "\n >>>>> Grabbing b-tag efficiency histos from each file..." << endl;

    TH1D* hBJet;
    TH1D* hBJetTagged;
    TH1D* hCJet;
    TH1D* hCJetTagged;
    TH1D* hLightJet;
    TH1D* hLightJetTagged;
    TH1D* hTemp;

    for (int i = 0; i < nFiles; i++) {
        for (int j = 0; j < nHistNoGen ; j++) {
            hist[i][j] = getHisto(file[i], histoName[j]);
            hTemp = hist[i][j];
            if (i == 0) {
                if (histoName.find("h_pt_b") != string::npos)            hBJet = (TH1D*) hist[i][j]->Clone();
                if (histoName.find("h_pt_b_tagged") != string::npos)     hBJetTagged = (TH1D*) hist[i][j]->Clone();
                if (histoName.find("h_pt_c") != string::npos)            hCJet = (TH1D*) hist[i][j]->Clone();
                if (histoName.find("h_pt_c_tagged") != string::npos)     hCJetTagged = (TH1D*) hist[i][j]->Clone();
                if (histoName.find("h_pt_udsg") != string::npos)         hLightJet = (TH1D*) hist[i][j]->Clone();
                if (histoName.find("h_pt_udsg_tagged") != string::npos)  hLightJetTagged = (TH1D*) hist[i][j]->Clone();


            }
            else {
                if (histoName.find("h_pt_b") != string::npos)            hBJet->Add(hTemp);
                if (histoName.find("h_pt_b_tagged") != string::npos)     hBJetTagged->Add(hTemp);
                if (histoName.find("h_pt_c") != string::npos)            hCJet->Add(hTemp);
                if (histoName.find("h_pt_c_tagged") != string::npos)     hCJetTagged->Add(hTemp);
                if (histoName.find("h_pt_udsg") != string::npos)         hLightJet->Add(hTemp);
                if (histoName.find("h_pt_udsg_tagged") != string::npos)  hLightJetTagged->Add(hTemp);
            }
        }

    }

    // Calculate the efficiencies
    TH1D* hBJetEff     = (TH1D*) hBJetTagged->Clone();
    TH1D* hCJetEff     = (TH1D*) hCJetTagged->Clone();
    TH1D* hLightJetEff = (TH1D*) hLightJetTagged->Clone();
    hBJetEff->Divide(hBJet);
    hCJetEff->Divide(hCJet);
    hLightJetEff->Divide(hLightJet);

    int nBins(hBJetEff->GetNbinsX());

    // BJets -----
    for (int i(1); i <= nBins; i++){
        cout << hBJetTagged->GetBinContent(i) << " ";
    }
    cout << endl;
    for (int i(1); i <= nBins; i++){
        cout << hBJet->GetBinContent(i) << " ";
    }
    cout << endl;
    for (int i(1); i <= nBins; i++){
        cout << hBJetEff->GetBinContent(i) << " ";
    }
    cout << endl;

    // Charm Jets -----
    for (int i(1); i <= nBins; i++){
        cout << hCJetTagged->GetBinContent(i) << " ";
    }
    cout << endl;
    for (int i(1); i <= nBins; i++){
        cout << hCJet->GetBinContent(i) << " ";
    }
    cout << endl;
    for (int i(1); i <= nBins; i++){
        cout << hCJetEff->GetBinContent(i) << " ";
    }
    cout << endl;

    // Light Jets -----
    for (int i(1); i <= nBins; i++){
        cout << hLightJetTagged->GetBinContent(i) << " ";
    }
    cout << endl;
    for (int i(1); i <= nBins; i++){
        cout << hLightJet->GetBinContent(i) << " ";
    }
    cout << endl;
    for (int i(1); i <= nBins; i++){
        cout << hLightJetEff->GetBinContent(i) << " ";
    }
    cout << endl;

    cout << "\n >>>>> Closing all files..." << endl;
    for (unsigned short i(0); i < nFiles; i++) closeFile(file[i]);

    cout << "\nExiting runCalcBTagEffsMC!" << endl;
}
