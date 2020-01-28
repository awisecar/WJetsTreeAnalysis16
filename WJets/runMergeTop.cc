#include <iostream>
#include <sstream>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TString.h>

void mergeTop(TString lepSelection = "SMu", int systematics = 0, int jetPtCutMin = 30, int doQCD = 0, int METcut = 0, int doBJets = -1);

void runMergeTop(){

    //central ---
    // mergeTop_BVeto("SMu", 0, 30, 0, 0, -1); // bveto of 1 jet 
    //mergeTop("SMu", 0, 30, 0, 0, 0); // no bveto
    mergeTop("SMu", 0, 30, 0, 0, 2); // >= 2 btags required

     // mergeTop("SMu", 0, 30, 0, 0, 0); // no bveto, doQCD=0
     // mergeTop("SMu", 0, 30, 1, 0, 0); // no bveto, doQCD=1
     // mergeTop("SMu", 0, 30, 2, 0, 0); // no bveto, doQCD=2
     // mergeTop("SMu", 0, 30, 3, 0, 0); // no bveto, doQCD=3

    //systematics ---
//     runMergeTop_BVeto("SMu", 1, 30, 0, 0);
//     runMergeTop_BVeto("SMu", -1, 30, 0, 0);
//     runMergeTop_BVeto("SMu", 3, 30, 0, 0);
//     runMergeTop_BVeto("SMu", -3, 30, 0, 0);
//   //  runMergeTop_BVeto("SMu", 5, 30, 0, 0);
//   //  runMergeTop_BVeto("SMu", -5, 30, 0, 0);
//     runMergeTop_BVeto("SMu", 6, 30, 0, 0);
//     runMergeTop_BVeto("SMu", -6, 30, 0, 0);

}

void mergeTop(TString lepSelection, int systematics, int jetPtCutMin, int doQCD, int METcut, int doBJets)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    cout << "\n-----> Running mergeTop!" << endl;
    //--- Read in input arguments
    ostringstream strJetPtCutMin; strJetPtCutMin << jetPtCutMin;
    ostringstream doQCDStr;       doQCDStr << doQCD;

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
  
    cout << "lepSelection = " << lepSelection << endl;
    cout << "syst = " << syst << endl;
    cout << "jetPtCutMin = " << jetPtCutMin << endl;
    cout << "doQCD = " << doQCD << endl;
    cout << "METcut = " << METcut << endl;
    cout << "doBJets = " << doBJets << endl;

    TString str1, str2, str3, str4, str5, strf;

    //--- Form file strings from input arguments    
    if (doQCD == 0) {
        str1 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_s_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth"+doBVetoStr+METcutStr+".root";
        str2 = "HistoFiles/"+ lepSelection + "_13TeV_ST_t_antitop_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth"+doBVetoStr+METcutStr+".root";
        str3 = "HistoFiles/"+ lepSelection + "_13TeV_ST_t_top_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth"+doBVetoStr+METcutStr+".root";
        str4 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_tW_top_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth"+doBVetoStr+METcutStr+".root";
        str5 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_tW_antitop_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth"+doBVetoStr+METcutStr+".root";
        strf = "HistoFiles/"+ lepSelection +  "_13TeV_Top_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth"+doBVetoStr+METcutStr+".root";
    }

    if (doQCD > 0){
        str1 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_s_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth"+doBVetoStr+"_QCD" + doQCDStr.str() + METcutStr+ ".root";
        str2 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_t_antitop_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth"+doBVetoStr+"_QCD" + doQCDStr.str() + METcutStr+ ".root";
        str3 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_t_top_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth"+doBVetoStr+"_QCD" + doQCDStr.str() + METcutStr+ ".root";
        str4 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_tW_top_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth"+doBVetoStr+"_QCD" + doQCDStr.str() + METcutStr+ ".root";
        str5 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_tW_antitop_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth"+doBVetoStr+"_QCD" + doQCDStr.str() + METcutStr+ ".root";
        strf = "HistoFiles/"+ lepSelection +  "_13TeV_Top_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth"+doBVetoStr+"_QCD" + doQCDStr.str() + METcutStr+ ".root";
    }

    cout << "\nInput files:\n" << str1 << "\n" << str2 << "\n" << str3 << "\n" << str4 << "\n" << str5 << std::endl;
    cout << "Output file: " << strf << "\n" << endl;

    //--- Open files
    TFile *f1 = new TFile(str1, "read");
    TFile *f2 = new TFile(str2, "read");
    TFile *f3 = new TFile(str3, "read");
    TFile *f4 = new TFile(str4, "read");
    TFile *f5 = new TFile(str5, "read");
    TFile *ff = new TFile(strf, "RECREATE");

    //--- Loop over every histo and add all indiv. single top contributions
    int nHist = f1->GetListOfKeys()->GetEntries();
    for (int i(0); i < nHist; i++){
        string hName = f1->GetListOfKeys()->At(i)->GetName();
        if (hName.find("hresponse") != string::npos){
            continue;
        }
        else {
         //   std::cout << "Doing histogram: " << hName << std::endl;
            TH1D *h1 = (TH1D*) f1->Get(hName.c_str()); 
            TH1D *h2 = (TH1D*) f2->Get(hName.c_str()); 
            TH1D *h3 = (TH1D*) f3->Get(hName.c_str()); 
            TH1D *h4 = (TH1D*) f4->Get(hName.c_str()); 
            TH1D *h5 = (TH1D*) f5->Get(hName.c_str()); 

            TH1D *hSum = (TH1D*) h1->Clone();
            hSum->Add(h2);
            hSum->Add(h3);
            hSum->Add(h4);
            hSum->Add(h5);

            //--- Write out summed histo to outfile
            ff->cd();
            hSum->Write();
        }
    }

    cout << "\nClosing single top files!" << endl;
    f1->Close();
    f2->Close();
    f3->Close();
    f4->Close();
    f5->Close();
    ff->Close();
}
