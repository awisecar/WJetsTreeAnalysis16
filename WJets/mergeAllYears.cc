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

void rescaleTTBarXSecs(string hName, TH1D* h2016, TH1D* h2017, TH1D* h2018);
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

    // Is the file ttBar?
    bool isTTBar(false);
    if (sampleName.Contains("TTJets")){
        isTTBar = true;
        std::cout << "\nInput file is ttBar! Will rescale below!" << std::endl;
    }

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

            // if files are ttBar, rescale histograms w/ xsec SFs
            // this rescaling is done individually for each year before adding to get full Run 2
            if (isTTBar) rescaleTTBarXSecs(hName, h1, h2, h3);

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

void rescaleTTBarXSecs(string hName, TH1D* h2016, TH1D* h2017, TH1D* h2018){

    // initialize ttBar SF values --
    double SFttbar_2016(1.), SFttbar_2017(1.), SFttbar_2018(1.);

    // list ttBar SFs as arrays --
    // 2016
    float ttBarSF_2016exc[7] = {1.058563,	1.004681,	0.995014,	0.972899,	0.948714,	0.939270,	0.947502};
    float ttBarSF_2016inc[7] = {1.002029,	0.991079,	0.981770,	0.963568,	0.946061,	0.939943,	0.941690};
    // 2017
    float ttBarSF_2017exc[7] = {1.363999,	1.113553,	1.070007,	1.071000,	1.094611,	1.132046,	1.198218};
    float ttBarSF_2017inc[7] = {1.138375,	1.091542,	1.076260,	1.085047,	1.112406,	1.154916,	1.217493};
    // 2018
    float ttBarSF_2018exc[7] = {1.123496,	0.997161,	0.955899,	0.950875,	0.947921,	0.981318,	1.032589};
    float ttBarSF_2018inc[7] = {0.997299,	0.972467,	0.955640,	0.955278,	0.963813,	1.001809,	1.058206};

    std::cout << "Rescaling: " << hName << std::endl;

    // 2016 ---
    // inclusive jet cut
    if (hName.find("Zinc2jet")!= string::npos)      SFttbar_2016 = ttBarSF_2016inc[0];
    else if (hName.find("Zinc3jet")!= string::npos) SFttbar_2016 = ttBarSF_2016inc[1];
    else if (hName.find("Zinc4jet")!= string::npos) SFttbar_2016 = ttBarSF_2016inc[2];
    else if (hName.find("Zinc5jet")!= string::npos) SFttbar_2016 = ttBarSF_2016inc[3];
    else if (hName.find("Zinc6jet")!= string::npos) SFttbar_2016 = ttBarSF_2016inc[4];
    else if (hName.find("Zinc7jet")!= string::npos) SFttbar_2016 = ttBarSF_2016inc[5];
    else if (hName.find("Zinc8jet")!= string::npos) SFttbar_2016 = ttBarSF_2016inc[6];
    // exclusive jet cut
    else if (hName.find("Zexc2jet")!= string::npos) SFttbar_2016 = ttBarSF_2016exc[0];
    else if (hName.find("Zexc3jet")!= string::npos) SFttbar_2016 = ttBarSF_2016exc[1];
    else if (hName.find("Zexc4jet")!= string::npos) SFttbar_2016 = ttBarSF_2016exc[2];
    else if (hName.find("Zexc5jet")!= string::npos) SFttbar_2016 = ttBarSF_2016exc[3];
    else if (hName.find("Zexc6jet")!= string::npos) SFttbar_2016 = ttBarSF_2016exc[4];
    else if (hName.find("Zexc7jet")!= string::npos) SFttbar_2016 = ttBarSF_2016exc[5];
    else if (hName.find("Zexc8jet")!= string::npos) SFttbar_2016 = ttBarSF_2016exc[6];
    // if not >= 2 jets cut, set xsec SF to 1
    else SFttbar_2016 = 1.;

    // 2017 --
    // inclusive jet cut
    if (hName.find("Zinc2jet")!= string::npos)      SFttbar_2017 = ttBarSF_2017inc[0];
    else if (hName.find("Zinc3jet")!= string::npos) SFttbar_2017 = ttBarSF_2017inc[1];
    else if (hName.find("Zinc4jet")!= string::npos) SFttbar_2017 = ttBarSF_2017inc[2];
    else if (hName.find("Zinc5jet")!= string::npos) SFttbar_2017 = ttBarSF_2017inc[3];
    else if (hName.find("Zinc6jet")!= string::npos) SFttbar_2017 = ttBarSF_2017inc[4];
    else if (hName.find("Zinc7jet")!= string::npos) SFttbar_2017 = ttBarSF_2017inc[5];
    else if (hName.find("Zinc8jet")!= string::npos) SFttbar_2017 = ttBarSF_2017inc[6];
    // exclusive jet cut
    else if (hName.find("Zexc2jet")!= string::npos) SFttbar_2017 = ttBarSF_2017exc[0];
    else if (hName.find("Zexc3jet")!= string::npos) SFttbar_2017 = ttBarSF_2017exc[1];
    else if (hName.find("Zexc4jet")!= string::npos) SFttbar_2017 = ttBarSF_2017exc[2];
    else if (hName.find("Zexc5jet")!= string::npos) SFttbar_2017 = ttBarSF_2017exc[3];
    else if (hName.find("Zexc6jet")!= string::npos) SFttbar_2017 = ttBarSF_2017exc[4];
    else if (hName.find("Zexc7jet")!= string::npos) SFttbar_2017 = ttBarSF_2017exc[5];
    else if (hName.find("Zexc8jet")!= string::npos) SFttbar_2017 = ttBarSF_2017exc[6];
    // if not >= 2 jets cut, set xsec SF to 1
    else SFttbar_2017 = 1.;

    // 2018 --
    // inclusive jet cut
    if (hName.find("Zinc2jet")!= string::npos)      SFttbar_2018 = ttBarSF_2018inc[0];
    else if (hName.find("Zinc3jet")!= string::npos) SFttbar_2018 = ttBarSF_2018inc[1];
    else if (hName.find("Zinc4jet")!= string::npos) SFttbar_2018 = ttBarSF_2018inc[2];
    else if (hName.find("Zinc5jet")!= string::npos) SFttbar_2018 = ttBarSF_2018inc[3];
    else if (hName.find("Zinc6jet")!= string::npos) SFttbar_2018 = ttBarSF_2018inc[4];
    else if (hName.find("Zinc7jet")!= string::npos) SFttbar_2018 = ttBarSF_2018inc[5];
    else if (hName.find("Zinc8jet")!= string::npos) SFttbar_2018 = ttBarSF_2018inc[6];
    // exclusive jet cut
    else if (hName.find("Zexc2jet")!= string::npos) SFttbar_2018 = ttBarSF_2018exc[0];
    else if (hName.find("Zexc3jet")!= string::npos) SFttbar_2018 = ttBarSF_2018exc[1];
    else if (hName.find("Zexc4jet")!= string::npos) SFttbar_2018 = ttBarSF_2018exc[2];
    else if (hName.find("Zexc5jet")!= string::npos) SFttbar_2018 = ttBarSF_2018exc[3];
    else if (hName.find("Zexc6jet")!= string::npos) SFttbar_2018 = ttBarSF_2018exc[4];
    else if (hName.find("Zexc7jet")!= string::npos) SFttbar_2018 = ttBarSF_2018exc[5];
    else if (hName.find("Zexc8jet")!= string::npos) SFttbar_2018 = ttBarSF_2018exc[6];
    // if not >= 2 jets cut, set xsec SF to 1
    else SFttbar_2018 = 1.;

    //the two sections below scale the jet multiplicities
    //bin number 3 is 2 jet mult. (histogram bin number 1 is 0 jet bin)
    if (hName.find("ZNGoodJets_Zinc") != string::npos){
        // 2016 ---
        // bin contents
        h2016->SetBinContent(3, h2016->GetBinContent(3)*ttBarSF_2016inc[0]);
        h2016->SetBinContent(4, h2016->GetBinContent(4)*ttBarSF_2016inc[1]);
        h2016->SetBinContent(5, h2016->GetBinContent(5)*ttBarSF_2016inc[2]);
        h2016->SetBinContent(6, h2016->GetBinContent(6)*ttBarSF_2016inc[3]);
        h2016->SetBinContent(7, h2016->GetBinContent(7)*ttBarSF_2016inc[4]);
        h2016->SetBinContent(8, h2016->GetBinContent(8)*ttBarSF_2016inc[5]);
        h2016->SetBinContent(9, h2016->GetBinContent(9)*ttBarSF_2016inc[6]);
        // bin errors 
        h2016->SetBinError(3, h2016->GetBinError(3)*ttBarSF_2016inc[0]);
        h2016->SetBinError(4, h2016->GetBinError(4)*ttBarSF_2016inc[1]);
        h2016->SetBinError(5, h2016->GetBinError(5)*ttBarSF_2016inc[2]);
        h2016->SetBinError(6, h2016->GetBinError(6)*ttBarSF_2016inc[3]);
        h2016->SetBinError(7, h2016->GetBinError(7)*ttBarSF_2016inc[4]);
        h2016->SetBinError(8, h2016->GetBinError(8)*ttBarSF_2016inc[5]);
        h2016->SetBinError(9, h2016->GetBinError(9)*ttBarSF_2016inc[6]);
        
        // 2017 --
        // bin contents
        h2017->SetBinContent(3, h2017->GetBinContent(3)*ttBarSF_2017inc[0]);
        h2017->SetBinContent(4, h2017->GetBinContent(4)*ttBarSF_2017inc[1]);
        h2017->SetBinContent(5, h2017->GetBinContent(5)*ttBarSF_2017inc[2]);
        h2017->SetBinContent(6, h2017->GetBinContent(6)*ttBarSF_2017inc[3]);
        h2017->SetBinContent(7, h2017->GetBinContent(7)*ttBarSF_2017inc[4]);
        h2017->SetBinContent(8, h2017->GetBinContent(8)*ttBarSF_2017inc[5]);
        h2017->SetBinContent(9, h2017->GetBinContent(9)*ttBarSF_2017inc[6]);
        // bin errors 
        h2017->SetBinError(3, h2017->GetBinError(3)*ttBarSF_2017inc[0]);
        h2017->SetBinError(4, h2017->GetBinError(4)*ttBarSF_2017inc[1]);
        h2017->SetBinError(5, h2017->GetBinError(5)*ttBarSF_2017inc[2]);
        h2017->SetBinError(6, h2017->GetBinError(6)*ttBarSF_2017inc[3]);
        h2017->SetBinError(7, h2017->GetBinError(7)*ttBarSF_2017inc[4]);
        h2017->SetBinError(8, h2017->GetBinError(8)*ttBarSF_2017inc[5]);
        h2017->SetBinError(9, h2017->GetBinError(9)*ttBarSF_2017inc[6]);
        
        // 2018 --
        // bin contents
        h2018->SetBinContent(3, h2018->GetBinContent(3)*ttBarSF_2018inc[0]);
        h2018->SetBinContent(4, h2018->GetBinContent(4)*ttBarSF_2018inc[1]);
        h2018->SetBinContent(5, h2018->GetBinContent(5)*ttBarSF_2018inc[2]);
        h2018->SetBinContent(6, h2018->GetBinContent(6)*ttBarSF_2018inc[3]);
        h2018->SetBinContent(7, h2018->GetBinContent(7)*ttBarSF_2018inc[4]);
        h2018->SetBinContent(8, h2018->GetBinContent(8)*ttBarSF_2018inc[5]);
        h2018->SetBinContent(9, h2018->GetBinContent(9)*ttBarSF_2018inc[6]);
        // bin errors 
        h2018->SetBinError(3, h2018->GetBinError(3)*ttBarSF_2018inc[0]);
        h2018->SetBinError(4, h2018->GetBinError(4)*ttBarSF_2018inc[1]);
        h2018->SetBinError(5, h2018->GetBinError(5)*ttBarSF_2018inc[2]);
        h2018->SetBinError(6, h2018->GetBinError(6)*ttBarSF_2018inc[3]);
        h2018->SetBinError(7, h2018->GetBinError(7)*ttBarSF_2018inc[4]);
        h2018->SetBinError(8, h2018->GetBinError(8)*ttBarSF_2018inc[5]);
        h2018->SetBinError(9, h2018->GetBinError(9)*ttBarSF_2018inc[6]);
        
    }
    else if (hName.find("ZNGoodJets_Zexc") != string::npos){
        // 2016 --
        // bin contents
        h2016->SetBinContent(3, h2016->GetBinContent(3)*ttBarSF_2016exc[0]);
        h2016->SetBinContent(4, h2016->GetBinContent(4)*ttBarSF_2016exc[1]);
        h2016->SetBinContent(5, h2016->GetBinContent(5)*ttBarSF_2016exc[2]);
        h2016->SetBinContent(6, h2016->GetBinContent(6)*ttBarSF_2016exc[3]);
        h2016->SetBinContent(7, h2016->GetBinContent(7)*ttBarSF_2016exc[4]);
        h2016->SetBinContent(8, h2016->GetBinContent(8)*ttBarSF_2016exc[5]);
        h2016->SetBinContent(9, h2016->GetBinContent(9)*ttBarSF_2016exc[6]);
        // bin errors 
        h2016->SetBinError(3, h2016->GetBinError(3)*ttBarSF_2016exc[0]);
        h2016->SetBinError(4, h2016->GetBinError(4)*ttBarSF_2016exc[1]);
        h2016->SetBinError(5, h2016->GetBinError(5)*ttBarSF_2016exc[2]);
        h2016->SetBinError(6, h2016->GetBinError(6)*ttBarSF_2016exc[3]);
        h2016->SetBinError(7, h2016->GetBinError(7)*ttBarSF_2016exc[4]);
        h2016->SetBinError(8, h2016->GetBinError(8)*ttBarSF_2016exc[5]);
        h2016->SetBinError(9, h2016->GetBinError(9)*ttBarSF_2016exc[6]);
        
        // 2017 --
        // bin contents
        h2017->SetBinContent(3, h2017->GetBinContent(3)*ttBarSF_2017exc[0]);
        h2017->SetBinContent(4, h2017->GetBinContent(4)*ttBarSF_2017exc[1]);
        h2017->SetBinContent(5, h2017->GetBinContent(5)*ttBarSF_2017exc[2]);
        h2017->SetBinContent(6, h2017->GetBinContent(6)*ttBarSF_2017exc[3]);
        h2017->SetBinContent(7, h2017->GetBinContent(7)*ttBarSF_2017exc[4]);
        h2017->SetBinContent(8, h2017->GetBinContent(8)*ttBarSF_2017exc[5]);
        h2017->SetBinContent(9, h2017->GetBinContent(9)*ttBarSF_2017exc[6]);
        // bin errors 
        h2017->SetBinError(3, h2017->GetBinError(3)*ttBarSF_2017exc[0]);
        h2017->SetBinError(4, h2017->GetBinError(4)*ttBarSF_2017exc[1]);
        h2017->SetBinError(5, h2017->GetBinError(5)*ttBarSF_2017exc[2]);
        h2017->SetBinError(6, h2017->GetBinError(6)*ttBarSF_2017exc[3]);
        h2017->SetBinError(7, h2017->GetBinError(7)*ttBarSF_2017exc[4]);
        h2017->SetBinError(8, h2017->GetBinError(8)*ttBarSF_2017exc[5]);
        h2017->SetBinError(9, h2017->GetBinError(9)*ttBarSF_2017exc[6]);
        
        // 2018 --
        // bin contents
        h2018->SetBinContent(3, h2018->GetBinContent(3)*ttBarSF_2018exc[0]);
        h2018->SetBinContent(4, h2018->GetBinContent(4)*ttBarSF_2018exc[1]);
        h2018->SetBinContent(5, h2018->GetBinContent(5)*ttBarSF_2018exc[2]);
        h2018->SetBinContent(6, h2018->GetBinContent(6)*ttBarSF_2018exc[3]);
        h2018->SetBinContent(7, h2018->GetBinContent(7)*ttBarSF_2018exc[4]);
        h2018->SetBinContent(8, h2018->GetBinContent(8)*ttBarSF_2018exc[5]);
        h2018->SetBinContent(9, h2018->GetBinContent(9)*ttBarSF_2018exc[6]);
        // bin errors 
        h2018->SetBinError(3, h2018->GetBinError(3)*ttBarSF_2018exc[0]);
        h2018->SetBinError(4, h2018->GetBinError(4)*ttBarSF_2018exc[1]);
        h2018->SetBinError(5, h2018->GetBinError(5)*ttBarSF_2018exc[2]);
        h2018->SetBinError(6, h2018->GetBinError(6)*ttBarSF_2018exc[3]);
        h2018->SetBinError(7, h2018->GetBinError(7)*ttBarSF_2018exc[4]);
        h2018->SetBinError(8, h2018->GetBinError(8)*ttBarSF_2018exc[5]);
        h2018->SetBinError(9, h2018->GetBinError(9)*ttBarSF_2018exc[6]);
        
    }
    else {
        //it is this line where we scale the ttbar MC for each of the histograms
        //if the name of the histogram does not contain one of the phrases above, it is simply scaled by 1
        h2016->Scale(SFttbar_2016);
        h2017->Scale(SFttbar_2017);
        h2018->Scale(SFttbar_2018);
    }

}