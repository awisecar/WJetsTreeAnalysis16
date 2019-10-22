#include <vector>
#include <TDirectory.h>
#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TPaveStats.h>
#include <TText.h>
#include <TLegend.h>
#include <iostream>
#include <sstream>
#include <TLorentzVector.h>
#include <TSVDUnfold.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TDatime.h>
#include <TRandom3.h>
#include <TMatrixD.h>

#include "getFilesAndHistograms.h"
#include "DataDrivenQCD.h"
#include "fileNames.h"

using namespace std;

const int NQCD = 4;
const int NMC = 13;

//string energy = getEnergy();
string energy = "13TeV";

int JetPtMin(30);
int JetPtMax(0);

//-------------------------------------------------------------------------------------------
void DataDrivenQCD(string leptonFlavor, int METcut , int doBJets){
    
    TH1::SetDefaultSumw2();
    
    // Opening all files
    TFile *fData[NQCD] = {NULL};
    TFile *fMC[NQCD][NMC] = {{NULL}};
    FuncOpenAllFiles(fData, fMC, leptonFlavor, METcut, false, true, doBJets);
    vector<string> histoNameRun = getVectorOfHistoNames(fData);
    
    // Creating output file
    string nameQCDout = fData[0]->GetName();
    nameQCDout.insert(nameQCDout.find("Data") + 4,"QCD");
    TFile *fOut = new TFile(nameQCDout.c_str(), "RECREATE");
    
    // Derive QCD BG for each histogram
    std::cout << "\n-----> Start QCD BG derivation!" << std::endl;
    for (int i(0); i < int(histoNameRun.size()) ; i++){
        
        if ( (histoNameRun[i].find("ZNGoodJets_") == string::npos) ) continue; // 22 oct 2019 -- just looking at jet mult for right now

        cout << "\n------------------------------------------------------------------------------" << endl;
        cout << " >>>>>>> Processing histogram #" << i << " : " << histoNameRun[i] << endl;
        FuncDataDrivenQCD(histoNameRun[i], fData, fMC, fOut);
    }
    
    //-- Close all the files ------------------------------
    cout << "\n\nClosing files..." << endl;
    for (int i(0); i < NQCD; i++) {
        closeFile(fData[i]);
        for (int j(0); j < NMC; j++){
            closeFile(fMC[i][j]);
        }
    }
    fOut ->Close();
    cout << "\nQCD calculation is done! All files are closed." << endl;
    //-----------------------------------------------------
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void FuncOpenAllFiles(TFile *fData[], TFile *fMC[][NMC], string leptonFlavor, int METcut, bool doFlat , bool doVarWidth, int doBJets){
    // Get data files
    std::cout << "\n --- Getting data files --- " << std::endl;
    for ( int i = 0 ; i < NQCD ; i++){
        std::cout << "Getting " << ProcessInfo[DATAFILENAME].filename << " for NQCD = " << i << std::endl;
        fData[i] = getFile(FILESDIRECTORY, leptonFlavor, energy, ProcessInfo[DATAFILENAME].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, i, 0, 0, METcut, doBJets, "", "0", false, false);
    }
    /// get MC files
    std::cout << "\n --- Getting MC files --- " << std::endl;
    for ( int i=0 ; i < NQCD ; i++){
        std::cout << "\n-----> Doing NQCD = " << i << std::endl;
        for ( int j = 0 ; j < NMC ; j++){
            string FilenameTemp;
            //if (j == 0) FilenameTemp = "WJetsALL_MIX_UNFOLDING_dR_5311_List";
            //note: is this the wjets file that needs to be used? (other is MLM)
            //comment: believe so, MLM is matching/merching for LO and FxFx is for NLO

            //if (j == 0) FilenameTemp = "WJets_FxFx_dR_5311_List";
            if (j == 0) FilenameTemp = "WJets_FxFx_012J_dR_5311_List";
            // if (j == 0) FilenameTemp = "WJets_FxFx_Wpt_dR_5311_List";
            if (j == 1) FilenameTemp = "DYJets50toInf_dR_5311_List";
            if (j == 2) FilenameTemp = "TT_FullHad_dR_5311_List";
            if (j == 3) FilenameTemp = "TT_SemiLep_dR_5311_List";
            if (j == 4) FilenameTemp = "TT_2L2Nu_dR_5311_List";
            if (j == 5) FilenameTemp = "ST_s_channel_dR_5311_List";
            if (j == 6) FilenameTemp = "ST_t_antitop_channel_dR_5311_List";
            if (j == 7) FilenameTemp = "ST_t_top_channel_dR_5311_List";
            if (j == 8) FilenameTemp = "ST_tW_top_channel_dR_5311_List";
            if (j == 9) FilenameTemp = "ST_tW_antitop_channel_dR_5311_List";
            if (j == 10) FilenameTemp = "WW_dR_5311_List";
            if (j == 11) FilenameTemp = "WZ_dR_5311_List";
            if (j == 12) FilenameTemp = "ZZ_dR_5311_List";

            std::cout << "Getting " << FilenameTemp << std::endl;
            
            fMC[i][j] = getFile(FILESDIRECTORY, leptonFlavor, energy, FilenameTemp, JetPtMin, JetPtMax, doFlat, doVarWidth, i , 0, 0, METcut, doBJets, "", "0", false, false);
            TH1D *hTemp2 = getHisto(fMC[i][j], "ZNGoodJets_Zexc");
            cout << "Checking the integral of ZNGoodJets_Zexc in this file: " << hTemp2 ->Integral() << endl;
        }
    }
    cout << endl;
    cout << "\n --- Opened all data and MC files! ---"<< endl;
}

vector<string> getVectorOfHistoNames(TFile *fData[]){
    cout << "\n --- Checking the histogram names to be studied ---" << endl;
    unsigned short nHist = fData[0]->GetListOfKeys()->GetEntries();
    vector<string> histoName;
    int countHist = 0 ;
    for (unsigned short i(0); i < nHist; i++) {
        string histoNameTemp = fData[0]->GetListOfKeys()->At(i)->GetName();
        TH1D* histTemp = (TH1D*) fData[0]->Get(histoNameTemp.c_str());
        if (histTemp->GetEntries() < 1) continue;
        if (histTemp->InheritsFrom(TH1D::Class())) {
            // Add histogram name to the list if it is not empty and if it inherits from TH1D
            histoName.push_back(fData[0]->GetListOfKeys()->At(i)->GetName());
            countHist++;
        }
    }
    cout << "We will produce " << histoName.size() << " QCD histograms! " << endl;
    // cout << "--------------------------------------------------" << endl;
    return histoName;
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FuncDataDrivenQCD(string variable, TFile *fData[], TFile *fMC[][NMC], TFile *fOut){
    TH1D  *hData[NQCD], *hSignal[NQCD], *hBack[NQCD];

    cout << "\n-----> Opening histograms from Data files" << endl;
    for (int i=0 ; i < NQCD ; i++){
        TH1D *hTemp = getHisto(fData[i], variable);
        hData[i] = (TH1D *) hTemp->Clone();
        cout << "Got data of histo name: " << variable << " for QCD region: " << i << " " << ", Integral: " << hData[i]->Integral() << endl;
        delete hTemp;
    }
    cout << " --- Retrieved Data histos ---" << endl;
    
    cout << "\n-----> Opening histograms from all the MC files" << endl;
    for ( int i=0 ; i < NQCD ; i++){
        cout << "   --- QCD region #" << i << " --- " << endl;
        for ( int j = 0 ; j < NMC ; j++){
            string FilenameTemp;
            // if (j == 0) FilenameTemp =  "WJets_FxFx_Wpt_dR_5311_List";
            if (j == 0) FilenameTemp =  "WJets_FxFx_012J_dR_5311_List";
            if (j == 1) FilenameTemp =  "DYJets50toInf_dR_5311_List";
            if (j == 2) FilenameTemp =  "TT_FullHad_dR_5311_List";
            if (j == 3) FilenameTemp =  "TT_SemiLep_dR_5311_List";
            if (j == 4) FilenameTemp =  "TT_2L2Nu_dR_5311_List";
            if (j == 5) FilenameTemp =  "ST_s_channel_dR_5311_List";
            if (j == 6) FilenameTemp =  "ST_t_antitop_channel_dR_5311_List";
            if (j == 7) FilenameTemp =  "ST_t_top_channel_dR_5311_List";
            if (j == 8) FilenameTemp =  "ST_tW_top_channel_dR_5311_List";
            if (j == 9) FilenameTemp =  "ST_tW_antitop_channel_dR_5311_List";
            if (j == 10) FilenameTemp =  "WW_dR_5311_List";
            if (j == 11) FilenameTemp =  "WZ_dR_5311_List";
            if (j == 12) FilenameTemp =  "ZZ_dR_5311_List";
            
            TH1D *hTemp1 = getHisto(fMC[i][j], variable);
            cout << "File: " << FilenameTemp << endl;
            if (FilenameTemp.find("WJets") != string::npos) {
                cout << "   This is signal: " << FilenameTemp << endl;
                hSignal[i] = (TH1D *) hTemp1->Clone();
            }
            else{
                if (j == 1) {
                    cout << "   This is background: " << FilenameTemp << endl;
                    hBack[i] = (TH1D *) hTemp1->Clone();
                }
                else {
                    cout << "   Add another background: " << FilenameTemp << endl;
                    hBack[i]->Add(hTemp1);
                }
            }
            cout << "   Got MC histograms for QCD region #" << i << "; MC type: " << j << "; Integral: " << hTemp1->Integral() << endl;
            delete hTemp1;
        }
    }
    cout << " --- Retrieved MC histos ---" << endl;
    
    //--- QCD estiamtion procedure -------------
    // Note: Currently QCD regions are indexed as such:
    // QCD signal region A = 0
    // QCD control region B (phase space where MT < 50 GeV and Iso < 0.15) = 1
    // QCD control region C (phase space where MT > 50 GeV and Iso > 0.20) = 2
    // QCD control region D (phase space where MT < 50 GeV and Iso > 0.20) = 3
    // Region C is used to derive shape, and regions D and B are used to determine normalization

    cout << "\n --- Start QCD background estimation! --- " << endl;
    
    // Followed as described in AN2012_331_v20 pages 7-9
    TH1D* scaledMC[NQCD], *hQCD[NQCD];
    double NormFactor_o[15], NormFactorISO_o[15], NormFactIsoError_o[15];
    
    // For all histograms that are not jet multiplicity-binned
    if (variable.find("ZNGoodJets") == string::npos){
        
        // Step 0 : initial normalization of Wjets and data
        // this factor should stabilize towards 1 as QCD BG is iteratively derived
        double NormFactor = (hData[0]->Integral() - hBack[0]->Integral()) / hSignal[0]->Integral();
        cout << "QCD region #0: (non-QCD) BG-subtracted data divided by signal (# of events) = " << NormFactor << endl;
        
        // // for a quick look at QCD contribution in each region
        // for ( int i = 0  ; i < NQCD ; i++){
        //     cout << " >> Initial integral (Data - all MC), QCD region #" << i << ": " << hData[i]->Integral() - (hBack[i]->Integral() + hSignal[i]->Integral()) << ", Integral of signal: " << hSignal[i]->Integral() << endl;
        // }
        
        // 12 iterations to make sure the stability of NormFactor is achieved
        for ( int j = 0 ; j < 12; j++){
            cout << " >>> Iteration #" << j << endl;

            // recalculate NormFactor by including QCD[0]
            if ( j > 0 )  NormFactor = (hData[0]->Integral() - hBack[0]->Integral() - hQCD[0]->Integral()) / hSignal[0]->Integral();
            cout << " NormFactor f_W: " << NormFactor << endl;
            
            // Step 1: calculate QCD BG in all 4 regions (1 signal and 3 control) ---
            for ( int i = 0  ; i < NQCD ; i++){
                // currently using the NormFactor derived in QCD region 0 for the calculation in all regions
                scaledMC[i] = (TH1D*) hSignal[i]->Clone();
                scaledMC[i]->Scale(NormFactor) ;
                // QCD BG is the data subtracted by all non-QCD MC
                // as NormFactor decreases towards 1, scaledMC should just become hSignal
                hQCD[i] = (TH1D*) hData[i]->Clone();
                hQCD[i]->Add(scaledMC[i],-1);
                hQCD[i]->Add(hBack[i],-1);
                // set the negative bins to 0
                for ( int k = 0 ; k < hData[i]->GetNbinsX()+1 ; k++ ){
                    if ( hQCD[i]->GetBinContent(k) <= 0 ) {
                        hQCD[i]->SetBinContent(k, 0.);
                        hQCD[i]->SetBinError(k, 0.);
                    }
                }
            }
            
            // step 2: Isolation Fake Rate ---
            double NormFactorISO(0);
            double NormFactIsoError(0);
            double errhqcd1(0), errhqcd3(0), inthqcd1(0), inthqcd3(0);
            // i=1,3 are the regions below the MT=50 GeV cut
            if ( hQCD[3]->Integral() > 0 && hQCD[1]->Integral() > 0){
                inthqcd1 = hQCD[1]->IntegralAndError(1, hData[0]->GetNbinsX(), errhqcd1);
                inthqcd3 = hQCD[3]->IntegralAndError(1, hData[0]->GetNbinsX(), errhqcd3);
                
                // NormFactorISO is the ratio of integral(region B) to integral(region D)
                NormFactorISO = inthqcd1/inthqcd3;
                NormFactIsoError = NormFactorISO * sqrt( pow((errhqcd1/inthqcd1), 2) + pow((errhqcd3/inthqcd3), 2) );
            }

            cout << " Ratio of QCD contributions (region B to D), f_B/D: " << NormFactorISO << " +- " << NormFactIsoError << ", Percent error: " << (NormFactIsoError/NormFactorISO)*100. << endl;
            
            // step 3 : isolation fake-rate from step 2 is aplied to QCD[2] to get QCD[0]
            for (int m = 1; m <= hQCD[0]->GetNbinsX(); m++){
                hQCD[0]->SetBinContent(m, NormFactorISO * hQCD[2]->GetBinContent(m));
                if ( hQCD[0]->GetBinContent(m) > 0){
                    // this condition ensures that both NormFactorISO and hQCD[2] != 0
                    // error assigned to bins is propagated from the isolation transfer factor f_B/D and the statistics of the distribution that gives the shape
                    hQCD[0]->SetBinError(m, hQCD[0]->GetBinContent(m) * sqrt( pow((NormFactIsoError/NormFactorISO), 2) + pow((hQCD[2]->GetBinError(m)/hQCD[2]->GetBinContent(m)), 2) ));
                }
                else {
                    hQCD[0]->SetBinError(m, 0.);
                }
            }

            // Note: the only quantity that keeps updating is NormFactor f_W and we treat it as constant without additional uncertainty. 
            // So, the uncertainty in QCD[i] is determined by the operation in the last loop, although we do calculate it in every loop. 
            // For QCD[0], the uncertainty also include the error in NormFactorISO.
            
            NormFactorISO_o[0] = NormFactorISO;
            NormFactIsoError_o[0] = NormFactIsoError;
        }
        
        NormFactor_o[0] = NormFactor;
        cout << "\n >>>>> Final value of W+Jets normalization factor f_W: " << NormFactor_o[0] << endl;
        cout << " >>>>> Print out final error on normalization factor f_B/D for " << variable << endl;
        cout << " f_B/D = " << NormFactorISO_o[0] << " +- " << NormFactIsoError_o[0] << endl;
        
    }
    else{
        std::cout << "\n          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n          !!!!!!  Processing jet-binned histogram!  !!!!!! \n          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" << std::endl;
        
        ////////////////////////////////////////////////////////////////////////
        std::cout << "\n >>>>> Bin contents for hData histos" << std::endl;
        for (int i = 0; i < NQCD; i++){
            std::cout << "QCD region " << i << std::endl;
            for (int m = 1; m <= hData[0]->GetNbinsX(); m++){
                std::cout << hData[i]->GetBinContent(m) << "   " << std::endl;
            }
        }

        std::cout << "\n >>>>> Bin contents for hSignal histos" << std::endl;
        for (int i = 0; i < NQCD; i++){
            std::cout << "QCD region " << i << std::endl;
            for (int m = 1; m <= hSignal[0]->GetNbinsX(); m++){
                std::cout << hSignal[i]->GetBinContent(m) << "   " << std::endl;
            }
        }

        std::cout << "\n >>>>> Bin contents for hBack histos" << std::endl;
        for (int i = 0; i < NQCD; i++){
            std::cout << "QCD region " << i << std::endl;
            for (int m = 1; m <= hBack[0]->GetNbinsX(); m++){
                std::cout << hBack[i]->GetBinContent(m) << "   " << std::endl;
            }
        }
        ////////////////////////////////////////////////////////////////////////

        // for ZNGoodJets dist.
        for (int i = 0; i < NQCD; i++){
            scaledMC[i] = (TH1D*) hSignal[i]->Clone();
            hQCD[i] = (TH1D*) hData[i]->Clone();
        }
        
        // Computing QCD separately for each jet multiplicity bin
        for (int m = 1; m <= hData[0]->GetNbinsX(); m++){
            cout << " -------- Processing bin number: " << m << endl;

            // step 0 : initial normalization of Wjets and data
            double NormFactor(1);
            if(hSignal[0]->GetBinContent(m) > 0){
                NormFactor = (hData[0]->GetBinContent(m) - hBack[0]->GetBinContent(m)) / hSignal[0]->GetBinContent(m);
            }
            else NormFactor = 1.;
            
            for (int j = 0; j < 12; j++){
                if ( j > 0 )  {
                    if(hSignal[0]->GetBinContent(m) > 0 ){
                        NormFactor = (hData[0]->GetBinContent(m) - hBack[0]->GetBinContent(m) - hQCD[0]->GetBinContent(m) ) / hSignal[0]->GetBinContent(m) ;
                    }
                    else NormFactor = 1. ;
                }
                
                // step 1:
                for ( int i = 0  ; i < NQCD ; i++){
                    scaledMC[i]->SetBinContent(m, NormFactor * hSignal[i]->GetBinContent(m));
                    scaledMC[i]->SetBinError(m, NormFactor * hSignal[i]->GetBinError(m));
                    
                    hQCD[i]->SetBinContent(m, hData[i]->GetBinContent(m) - (scaledMC[i]->GetBinContent(m) + hBack[i]->GetBinContent(m)) );
                    hQCD[i]->SetBinError(m,  sqrt(pow(hData[i]->GetBinError(m), 2) + pow(scaledMC[i]->GetBinError(m), 2) + pow(hBack[i]->GetBinError(m), 2)) );
                    
                    // set the negative bins to 0
                    if (hQCD[i]->GetBinContent(m) <= 0) {
                        hQCD[i]->SetBinContent(m, 0. ) ;
                        hQCD[i]->SetBinError(m, 0. ) ;
                    }
                }

                // step 2: isoaltion fake rate
                double NormFactorISO(0);
                double NormFactIsoError(0);
                if ( hQCD[3]->GetBinContent(m) > 0 && hQCD[1]->GetBinContent(m) > 0 ) {
                    NormFactorISO = hQCD[1]->GetBinContent(m) / hQCD[3]->GetBinContent(m);
                    NormFactIsoError = NormFactorISO * sqrt( pow((hQCD[1]->GetBinError(m)/hQCD[1]->GetBinContent(m) ), 2) + pow((hQCD[3]->GetBinError(m)/hQCD[3]->GetBinContent(m)), 2) );
                }
                cout << " >>>>> NormFactor f_W: " << NormFactor << endl;
                cout << " Ratio of QCD regions B to D: " << NormFactorISO <<  " +- " << NormFactIsoError << ", % Error = " << NormFactIsoError*100/NormFactorISO << endl;
                
                // step 3 : isolation fake-rate from step 2 is aplied to QCD[2] to get QCD[0]
                hQCD[0]->SetBinContent(m, NormFactorISO * hQCD[2]->GetBinContent(m));
                if ( hQCD[0]->GetBinContent(m) > 0){
                    hQCD[0]->SetBinError(m, hQCD[0]->GetBinContent(m) * sqrt( pow((NormFactIsoError/NormFactorISO), 2) + pow((hQCD[2]->GetBinError(m)/hQCD[2]->GetBinContent(m)), 2)));
                }
                else hQCD[0]->SetBinError(m, 0.);
                
                // Note: the only quantity that keeps updating is NormFactor and we treat it as constant without uncertainty. 
                // So, the uncertainty in QCD[i] is determined by the operation in the last loop, although we do calculate it in every loop. 
                // For QCD[0], the uncertainty also include the error in NormFactorISO.
                
                NormFactorISO_o[m] = NormFactorISO;
                NormFactIsoError_o[m] = NormFactIsoError;
            }
            NormFactor_o[m] = NormFactor;
        }
        
        cout << "\n >>>>> Print out QCD BG statistics for " << variable << endl;
        cout << " \t" << "# QCD Events" << " \t\t" << "Error " << " \t\t" << "Percent Error" << " \t\t" << "Final f_W" << " \t\t" << "Final transfer factor f_B/D" << " \t\t" << "Percent Error in f_B/D" << endl;
        for (int i = 1; i <= hQCD[0]->GetNbinsX(); i++){
            cout << "Bin #" << i << endl;
            cout << "\t" <<  hQCD[0]->GetBinContent(i) << "\t" << hQCD[0]->GetBinError(i) << "\t" <<  (hQCD[0]->GetBinError(i)*100)/hQCD[0]->GetBinContent(i) << "\t" << NormFactor_o[i] << "\t" << NormFactorISO_o[i] << "\t" << (NormFactIsoError_o[i]*100)/NormFactorISO_o[i] << endl;
        }
    }
    
    // ---- Save derived QCD BG (signal region A) to file
    cout << " --- Now save results to file for: " << variable << endl;
    fOut->cd();
    hQCD[0]->Write();
    
}