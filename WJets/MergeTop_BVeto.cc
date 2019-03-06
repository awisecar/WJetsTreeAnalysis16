#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <sstream>

void runMergeTop_BVeto(string lepSelection = "DE", int systematics =0  , int jetPtCutMin = 30 , int doQCD = 0 );

void MergeTop_BVeto(){
    //central
  //  runMergeTop_BVeto("SMu", 0, 30, 0);
    //systematics
    runMergeTop_BVeto("SMu", 1, 30, 0);
    runMergeTop_BVeto("SMu", -1, 30, 0);
    runMergeTop_BVeto("SMu", 3, 30, 0);
    runMergeTop_BVeto("SMu", -3, 30, 0);
  //  runMergeTop_BVeto("SMu", 5, 30, 0);
  //  runMergeTop_BVeto("SMu", -5, 30, 0);
  //  runMergeTop_BVeto("SMu", 6, 30, 0);
  //  runMergeTop_BVeto("SMu", -6, 30, 0);
}

void runMergeTop_BVeto(string lepSelection, int systematics, int jetPtCutMin, int doQCD)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    cout << "\n-----> Running runMergeTop_BVeto!"
    //--- Read in input arguments
    ostringstream strJetPtCutMin; strJetPtCutMin << jetPtCutMin;
    ostringstream doQCDStr;     doQCDStr << doQCD ;
    string syst;
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

    string str1, str2, str3, str4, str5, strf;

    int nDYfiles = 2 ;
    string sstrDY[10];

    //--- Form file strings from input arguments    
    if (doQCD == 0) {
        str1 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_s_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto.root";
        str2 = "HistoFiles/"+ lepSelection + "_13TeV_ST_t_antitop_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto.root";
        str3 = "HistoFiles/"+ lepSelection + "_13TeV_ST_t_top_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto.root";
        str4 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_tW_top_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto.root";
        str5 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_tW_antitop_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto.root";
        strf = "HistoFiles/"+ lepSelection +  "_13TeV_Top_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto.root";
        
        sstrDY[0] = "HistoFiles/"+ lepSelection +  "_13TeV_DYJets50toInf_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto.root";
        sstrDY[1] = "HistoFiles/"+ lepSelection +  "_13TeV_DYJets10toInf3_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto.root";
    }

    if (doQCD > 0){
        str1 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_s_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto_QCD" + doQCDStr.str() + ".root";
        str2 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_t_antitop_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto_QCD" + doQCDStr.str() + ".root";
        str3 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_t_top_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto_QCD" + doQCDStr.str() + ".root";
        str4 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_tW_top_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto_QCD" + doQCDStr.str() + ".root";
        str5 = "HistoFiles/"+ lepSelection +  "_13TeV_ST_tW_antitop_channel_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto_QCD" + doQCDStr.str() + ".root";
        strf = "HistoFiles/"+ lepSelection +  "_13TeV_Top_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto_QCD" + doQCDStr.str() + ".root";
        
        sstrDY[0] = "HistoFiles/"+ lepSelection +  "_13TeV_DYJets50toInf_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str()+ "_VarWidth_BVeto_QCD" + doQCDStr.str() + ".root";
        sstrDY[1] = "HistoFiles/"+ lepSelection +  "_13TeV_DYJets10toInf3_dR_5311_List_EffiCorr_1_TrigCorr_1_" + syst + "JetPtMin_" + strJetPtCutMin.str() + "_VarWidth_BVeto_QCD" + doQCDStr.str() + ".root";
    }

    cout << "\nInput files: " << str1 << "\n" << str2 << "\n" << str3 << "\n" << str4 << "\n" << str5 << std::endl;
    cout << "Output file: " << strf << endl;

    //--- Open files
    TFile *f1 = new TFile(str1.c_str(), "read");
    TFile *f2 = new TFile(str2.c_str(), "read");
    TFile *f3 = new TFile(str3.c_str(), "read");
    TFile *f4 = new TFile(str4.c_str(), "read");
    TFile *f5 = new TFile(str5.c_str(), "read");
    TFile *ff = new TFile(strf.c_str(), "RECREATE");

    //--- Loop over every histo and add all indiv. single top contributions
    int nHist = f1->GetListOfKeys()->GetEntries();
    for (int i(0); i < nHist; i++){
        string hName = f1->GetListOfKeys()->At(i)->GetName();
        if (hName.find("hresponse") != string::npos){
            continue;
        }
        else {
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

//    //--- merge DY files
//    int countHist(0);
//    if (lepSelection == "SMu"){
//        TFile *fDY[10] = {NULL};
//        for ( int i = 0 ; i < nDYfiles ; i++){
//            if ( i == nDYfiles - 1 )  fDY[i] =  new TFile(sstrDY[i].c_str(), "recreate");
//            else fDY[i] =  new TFile(sstrDY[i].c_str(), "read");
//        }
//        cout << "Output file: " << sstrDY[nDYfiles - 1] << endl;
//
//        nHist = fDY[0]->GetListOfKeys()->GetEntries();
//        for (int i(0); i < nHist; i++){
//            
//            string hName = fDY[0]->GetListOfKeys()->At(i)->GetName();
//            if (hName.find("hresponse") != string::npos){
//                continue;
//                /*
//                cout << i << " TH2D " << hName << "  " << nHist << endl;
//                TH2D* hrSum = NULL; TH2D* hrDY[10] = {NULL};
//                for ( int j = 0 ; j < nDYfiles -1 ; j++){
//                    hrDY[j] = (TH2D*) fDY[j]->Get(hName.c_str());
//                    if ( j == 0 ) hrSum = (TH2D*) hrDY[j]->Clone();
//                    else hrSum->Add(hrDY[j]);
//                }
//                fDY[nDYfiles -1]->cd();
//                hrSum->Write();
//                */
//            }
//            else {
//                countHist++;
//                cout << countHist << " " << i << " TH1D " << hName << "  " << nHist << endl;
//                TH1D* hSum = NULL; TH1D* hDY[10] = {NULL};
//                for ( int j = 0 ; j < nDYfiles -1 ; j++){
//                    hDY[j] = (TH1D*) fDY[j]->Get(hName.c_str());
//                    if ( j == 0 ) hSum = (TH1D*) hDY[j]->Clone();
//                    else hSum->Add(hDY[j]);
//                }
//                fDY[nDYfiles -1]->cd();
//                hSum->Write();
//            }
//        }
//        
//        cout << "closing DY files" << endl;
//        for ( int i = 0 ; i < nDYfiles ; i++){
//            fDY[i] ->Close();
//        }
//    }
    
}
