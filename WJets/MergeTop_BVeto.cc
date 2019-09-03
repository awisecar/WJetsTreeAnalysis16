#include <iostream>
#include <sstream>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TString.h>

void runMergeTop_BVeto(TString lepSelection = "SMu", int systematics = 0, int jetPtCutMin = 30, int doQCD = 0, int METcut = 0);

void MergeTop_BVeto(){

    //central
   runMergeTop_BVeto("SMu", 0, 30, 0, 0);
//    runMergeTop_BVeto("SMu", 0, 30, 0, 30); // METcut of 30 GeV

    //systematics
//     runMergeTop_BVeto("SMu", 1, 30, 0, 0);
//     runMergeTop_BVeto("SMu", -1, 30, 0, 0);
//     runMergeTop_BVeto("SMu", 3, 30, 0, 0);
//     runMergeTop_BVeto("SMu", -3, 30, 0, 0);
//   //  runMergeTop_BVeto("SMu", 5, 30, 0, 0);
//   //  runMergeTop_BVeto("SMu", -5, 30, 0, 0);
//     runMergeTop_BVeto("SMu", 6, 30, 0, 0);
//     runMergeTop_BVeto("SMu", -6, 30, 0, 0);

}

void runMergeTop_BVeto(TString lepSelection, int systematics, int jetPtCutMin, int doQCD, int METcut)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    cout << "\n-----> Running runMergeTop_BVeto!" << endl;
    //--- Read in input arguments
    ostringstream strJetPtCutMin; strJetPtCutMin << jetPtCutMin;
    ostringstream doQCDStr;       doQCDStr << doQCD;
    ostringstream METcutStream;   METcutStream << METcut;
    TString METcutStr;
    if (METcut > 0) {
        METcutStr = "_MET"+METcutStream.str();
    }
    else {
        METcutStr = "";
    }
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

    TString str1, str2, str3, str4, str5, strf;

    // int nDYfiles = 2 ;
    string sstrDY[10];

    //--- Form file strings from input arguments    
    if (doQCD == 0) {
        str1 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_s_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto"+METcutStr+".root";
        str2 = "HistoFiles/"+ lepSelection + "_13TeV_ST_t_antitop_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto"+METcutStr+".root";
        str3 = "HistoFiles/"+ lepSelection + "_13TeV_ST_t_top_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto"+METcutStr+".root";
        str4 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_tW_top_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto"+METcutStr+".root";
        str5 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_tW_antitop_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto"+METcutStr+".root";
        strf = "HistoFiles/"+ lepSelection +  "_13TeV_Top_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto"+METcutStr+".root";
        
        sstrDY[0] = "HistoFiles/"+ lepSelection +  "_13TeV_DYJets50toInf_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto"+METcutStr+".root";
        sstrDY[1] = "HistoFiles/"+ lepSelection +  "_13TeV_DYJets10toInf3_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto"+METcutStr+".root";
    }

    if (doQCD > 0){
        str1 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_s_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto_QCD" + doQCDStr.str() + METcutStr+ ".root";
        str2 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_t_antitop_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto_QCD" + doQCDStr.str() + METcutStr+ ".root";
        str3 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_t_top_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto_QCD" + doQCDStr.str() + METcutStr+ ".root";
        str4 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_tW_top_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto_QCD" + doQCDStr.str() + METcutStr+ ".root";
        str5 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_tW_antitop_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto_QCD" + doQCDStr.str() + METcutStr+ ".root";
        strf = "HistoFiles/"+ lepSelection +  "_13TeV_Top_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto_QCD" + doQCDStr.str() + METcutStr+ ".root";
        
        sstrDY[0] = "HistoFiles/"+ lepSelection +  "_13TeV_DYJets50toInf_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str()+ "_VarWidth_BVeto_QCD" + doQCDStr.str() + METcutStr+ ".root";
        sstrDY[1] = "HistoFiles/"+ lepSelection +  "_13TeV_DYJets10toInf3_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto_QCD" + doQCDStr.str() + METcutStr+ ".root";
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
