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

void calcTheoryUncert_NLOFxFx_Ratios(){
    
    // ----------------------------------------------------------------
    // Setup

    // TString year       = "2016";
    TString year       = "2017";

    TString hNameNUM   = "LepPtPlusLeadingJetAK8Pt_Zinc2jet_TUnfold";
    TString hNameDENOM = "LepPtPlusLeadingJetAK8Pt_Zinc1jet_TUnfold";

    // ----------------------------------------------------------------

    // base filename for the NLO FxFx pt-binned samples
    TString filename_central = "HistoFiles_theoryUncert/"+year+"/"+"SMu_13TeV_WJets_FxFx_Wpt_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth";

    // Open nominal file and all theoretical uncertainty variations to get input histograms
    // -- nominal
    TFile *file_central  = TFile::Open(filename_central+".root", "READ");
    // -- alpha-s variations
    TFile *file_alphaS_1 = TFile::Open(filename_central+"_TheoryUncert_31.root", "READ");
    TFile *file_alphaS_2 = TFile::Open(filename_central+"_TheoryUncert_32.root", "READ");
    // -- PDF set variations
    TFile *file_PDF_1    = TFile::Open(filename_central+"_TheoryUncert_21.root", "READ");
    TFile *file_PDF_2    = TFile::Open(filename_central+"_TheoryUncert_22.root", "READ");
    // -- QCD scale variations
    TFile *file_scales_1 = TFile::Open(filename_central+"_TheoryUncert_11.root", "READ");
    TFile *file_scales_2 = TFile::Open(filename_central+"_TheoryUncert_12.root", "READ");
    TFile *file_scales_3 = TFile::Open(filename_central+"_TheoryUncert_13.root", "READ");
    TFile *file_scales_4 = TFile::Open(filename_central+"_TheoryUncert_14.root", "READ");
    TFile *file_scales_5 = TFile::Open(filename_central+"_TheoryUncert_16.root", "READ");
    TFile *file_scales_6 = TFile::Open(filename_central+"_TheoryUncert_18.root", "READ");

    // Open output file
    TFile *file_output = TFile::Open("HistoFiles_theoryUncert/theoryUncert_NLOFxFx_pTbinned_"+year+"_Ratios.root", "RECREATE");

    // ----------------------------------------------------------------

    std::cout << "\n\n >>> Analyzing ratio: " << hNameNUM << "_TO_" << hNameDENOM << std::endl;

    // grab central histogram first
    // NOTE: for error histograms, clone the main distribution to get the binning
    TH1D *h_central_NUM   = (TH1D*) file_central->Get("gen"+hNameNUM);
    TH1D *h_central_DENOM = (TH1D*) file_central->Get("gen"+hNameDENOM);

    // ----------------------------------------------------------------
    // alpha-s Uncertainty

    // -- set up the main error histograms
    // up variation
    TH1D *h_Err_AlphaS_Up = (TH1D*) h_central_NUM->Clone(hNameNUM+"_TO_"+hNameDENOM+"_ErrAlphaS_Up");
    h_Err_AlphaS_Up->Reset();
    h_Err_AlphaS_Up->SetTitle(hNameNUM+"_TO_"+hNameDENOM+"_ErrAlphaS_Up");
    h_Err_AlphaS_Up->GetYaxis()->SetRangeUser(0., 0.1);
    h_Err_AlphaS_Up->GetYaxis()->SetTitle("Relative Uncertainty (#alpha_{s}), Up Variation");
    // down variation
    TH1D *h_Err_AlphaS_Down = (TH1D*) h_central_NUM->Clone(hNameNUM+"_TO_"+hNameDENOM+"_ErrAlphaS_Down");
    h_Err_AlphaS_Down->Reset();
    h_Err_AlphaS_Down->SetTitle(hNameNUM+"_TO_"+hNameDENOM+"_ErrAlphaS_Down");
    h_Err_AlphaS_Down->GetYaxis()->SetRangeUser(0., 0.1);
    h_Err_AlphaS_Down->GetYaxis()->SetTitle("Relative Uncertainty (#alpha_{s}), Down Variation");

    if (file_alphaS_1 && file_alphaS_2){
        std::cout << "\n>>> alpha-s Uncertainty:" << std::endl;

        // get alpha-s uncertainty variations
        TH1D *h_alphaS_1_NUM   = (TH1D*) file_alphaS_1->Get("gen"+hNameNUM);
        TH1D *h_alphaS_1_DENOM = (TH1D*) file_alphaS_1->Get("gen"+hNameDENOM);
        TH1D *h_alphaS_2_NUM   = (TH1D*) file_alphaS_2->Get("gen"+hNameNUM);
        TH1D *h_alphaS_2_DENOM = (TH1D*) file_alphaS_2->Get("gen"+hNameDENOM);

        for (int iBin(1); iBin < h_central_NUM->GetNbinsX()+1; iBin++){

            double contentCentral = h_central_NUM->GetBinContent(iBin) / h_central_DENOM->GetBinContent(iBin);
            double contentVar1    = h_alphaS_1_NUM->GetBinContent(iBin) / h_alphaS_1_DENOM->GetBinContent(iBin);
            double contentVar2    = h_alphaS_2_NUM->GetBinContent(iBin) / h_alphaS_2_DENOM->GetBinContent(iBin);

            double fracErr1(0.), fracErr2(0.);
            if (contentCentral > 0){
                fracErr1   = fabs( (contentVar1 - contentCentral)/contentCentral ); // down variation
                fracErr2   = fabs( (contentVar2 - contentCentral)/contentCentral ); // up variation
            }
                    
            std::cout << ">>>" << std::endl;
            std::cout << "bin #" << iBin << ": fracErr1 = "   << fracErr1 << std::endl;
            std::cout << "bin #" << iBin << ": fracErr2 = "   << fracErr2 << std::endl;
            
            h_Err_AlphaS_Up->SetBinContent(iBin, fracErr2);
            h_Err_AlphaS_Down->SetBinContent(iBin, fracErr1);
        }
    }

    // ----------------------------------------------------------------
    // PDF Set Uncertainty

    // setup main error histograms
    // up variation
    TH1D *h_Err_PDF_Up = (TH1D*) h_central_NUM->Clone(hNameNUM+"_TO_"+hNameDENOM+"_ErrPDF_Up");
    h_Err_PDF_Up->Reset();
    h_Err_PDF_Up->SetTitle(hNameNUM+"_TO_"+hNameDENOM+"_ErrPDF_Up");
    h_Err_PDF_Up->GetYaxis()->SetRangeUser(0., 0.1);
    h_Err_PDF_Up->GetYaxis()->SetTitle("Relative Uncertainty (PDF), Up Variation");
    // down variation
    TH1D *h_Err_PDF_Down = (TH1D*) h_central_NUM->Clone(hNameNUM+"_TO_"+hNameDENOM+"_ErrPDF_Down");
    h_Err_PDF_Down->Reset();
    h_Err_PDF_Down->SetTitle(hNameNUM+"_TO_"+hNameDENOM+"_ErrPDF_Down");
    h_Err_PDF_Down->GetYaxis()->SetRangeUser(0., 0.1);
    h_Err_PDF_Down->GetYaxis()->SetTitle("Relative Uncertainty (PDF), Down Variation");

    if (file_PDF_1 && file_PDF_2){
        std::cout << "\n>>> PDF Uncertainty:" << std::endl;

        // get PDF uncertainty variations
        TH1D *h_PDF_1_NUM   = (TH1D*) file_PDF_1->Get("gen"+hNameNUM);
        TH1D *h_PDF_1_DENOM = (TH1D*) file_PDF_1->Get("gen"+hNameDENOM);
        TH1D *h_PDF_2_NUM   = (TH1D*) file_PDF_2->Get("gen"+hNameNUM);
        TH1D *h_PDF_2_DENOM = (TH1D*) file_PDF_2->Get("gen"+hNameDENOM);

        for (int iBin(1); iBin < h_central_NUM->GetNbinsX()+1; iBin++){

            double contentCentral = h_central_NUM->GetBinContent(iBin) / h_central_DENOM->GetBinContent(iBin);
            double contentVar1    = h_PDF_1_NUM->GetBinContent(iBin) / h_PDF_1_DENOM->GetBinContent(iBin);
            double contentVar2    = h_PDF_2_NUM->GetBinContent(iBin) / h_PDF_2_DENOM->GetBinContent(iBin);

            double fracErr1(0.), fracErr2(0.);
            if (contentCentral > 0){
                fracErr1   = fabs( (contentVar1 - contentCentral)/contentCentral ); // up variation
                fracErr2   = fabs( (contentVar2 - contentCentral)/contentCentral ); // down variation
            }
    
            std::cout << ">>>" << std::endl;
            std::cout << "bin #" << iBin << ": fracErr1 = " << fracErr1 << std::endl;
            std::cout << "bin #" << iBin << ": fracErr2 = " << fracErr2 << std::endl;

            h_Err_PDF_Up->SetBinContent(iBin, fracErr1);
            h_Err_PDF_Down->SetBinContent(iBin, fracErr2);

        }
    }

    // ----------------------------------------------------------------
    // QCD Scale Uncertainty

    // setup main error histograms
    // up variation
    TH1D *h_Err_Scale_Up = (TH1D*) h_central_NUM->Clone(hNameNUM+"_TO_"+hNameDENOM+"_ErrScale_Up");
    h_Err_Scale_Up->Reset();
    h_Err_Scale_Up->SetTitle(hNameNUM+"_TO_"+hNameDENOM+"_ErrScale_Up");
    h_Err_Scale_Up->GetYaxis()->SetRangeUser(0., 0.1);
    h_Err_Scale_Up->GetYaxis()->SetTitle("Relative Uncertainty (QCD Scales), Up Variation");
    // down variation
    TH1D *h_Err_Scale_Down = (TH1D*) h_central_NUM->Clone(hNameNUM+"_TO_"+hNameDENOM+"_ErrScale_Down");
    h_Err_Scale_Down->Reset();
    h_Err_Scale_Down->SetTitle(hNameNUM+"_TO_"+hNameDENOM+"_ErrScale_Down");
    h_Err_Scale_Down->GetYaxis()->SetRangeUser(0., 0.1);
    h_Err_Scale_Down->GetYaxis()->SetTitle("Relative Uncertainty (QCD Scales), Down Variation");

    if (file_scales_1 && file_scales_2 && file_scales_3 && file_scales_4 && file_scales_5 && file_scales_6){
        std::cout <<"\n>>> Scale Uncertainty:" << std::endl;

        // get scale uncertainty variations
        // --- NUM
        TH1D *h_scales_1_NUM   = (TH1D*) file_scales_1->Get("gen"+hNameNUM);
        TH1D *h_scales_2_NUM   = (TH1D*) file_scales_2->Get("gen"+hNameNUM);
        TH1D *h_scales_3_NUM   = (TH1D*) file_scales_3->Get("gen"+hNameNUM);
        TH1D *h_scales_4_NUM   = (TH1D*) file_scales_4->Get("gen"+hNameNUM);
        TH1D *h_scales_5_NUM   = (TH1D*) file_scales_5->Get("gen"+hNameNUM);
        TH1D *h_scales_6_NUM   = (TH1D*) file_scales_6->Get("gen"+hNameNUM);
        // --- DENOM
        TH1D *h_scales_1_DENOM = (TH1D*) file_scales_1->Get("gen"+hNameDENOM);
        TH1D *h_scales_2_DENOM = (TH1D*) file_scales_2->Get("gen"+hNameDENOM);
        TH1D *h_scales_3_DENOM = (TH1D*) file_scales_3->Get("gen"+hNameDENOM);
        TH1D *h_scales_4_DENOM = (TH1D*) file_scales_4->Get("gen"+hNameDENOM);
        TH1D *h_scales_5_DENOM = (TH1D*) file_scales_5->Get("gen"+hNameDENOM);
        TH1D *h_scales_6_DENOM = (TH1D*) file_scales_6->Get("gen"+hNameDENOM);

        // out of all of the 6 scale variations, need to take the 
        // outer borders of the total envelope as the uncertainty
        for (int iBin(1); iBin < h_central_NUM->GetNbinsX()+1; iBin++){

            // --- do first for NUM
            double contentCentral = h_central_NUM->GetBinContent(iBin);
            double contentVar1    = h_scales_1_NUM->GetBinContent(iBin);
            double contentVar2    = h_scales_2_NUM->GetBinContent(iBin);
            double contentVar3    = h_scales_3_NUM->GetBinContent(iBin);
            double contentVar4    = h_scales_4_NUM->GetBinContent(iBin);
            double contentVar5    = h_scales_5_NUM->GetBinContent(iBin);
            double contentVar6    = h_scales_6_NUM->GetBinContent(iBin);

            double envelopeMax_NUM(0.), envelopeMin_NUM(0.);
            if (contentCentral > 0){

                std::vector<double> scaleVar;
                scaleVar.push_back(contentVar1/contentCentral);
                scaleVar.push_back(contentVar2/contentCentral);
                scaleVar.push_back(contentVar3/contentCentral);
                scaleVar.push_back(contentVar4/contentCentral);
                scaleVar.push_back(contentVar5/contentCentral);
                scaleVar.push_back(contentVar6/contentCentral);

                envelopeMax_NUM = *max_element(scaleVar.begin(), scaleVar.end());
                envelopeMin_NUM = *min_element(scaleVar.begin(), scaleVar.end());
            }

            // --- then do for DENOM
            contentCentral = h_central_DENOM->GetBinContent(iBin);
            contentVar1    = h_scales_1_DENOM->GetBinContent(iBin);
            contentVar2    = h_scales_2_DENOM->GetBinContent(iBin);
            contentVar3    = h_scales_3_DENOM->GetBinContent(iBin);
            contentVar4    = h_scales_4_DENOM->GetBinContent(iBin);
            contentVar5    = h_scales_5_DENOM->GetBinContent(iBin);
            contentVar6    = h_scales_6_DENOM->GetBinContent(iBin);

            double envelopeMax_DENOM(0.), envelopeMin_DENOM(0.);
            if (contentCentral > 0){

                std::vector<double> scaleVar;
                scaleVar.push_back(contentVar1/contentCentral);
                scaleVar.push_back(contentVar2/contentCentral);
                scaleVar.push_back(contentVar3/contentCentral);
                scaleVar.push_back(contentVar4/contentCentral);
                scaleVar.push_back(contentVar5/contentCentral);
                scaleVar.push_back(contentVar6/contentCentral);

                envelopeMax_DENOM = *max_element(scaleVar.begin(), scaleVar.end());
                envelopeMin_DENOM = *min_element(scaleVar.begin(), scaleVar.end());
            }

            double envelopeMax_RATIO(0.), envelopeMin_RATIO(0.);
            if ( (envelopeMax_DENOM > 0.) && (envelopeMin_DENOM > 0.) ){
                envelopeMax_RATIO = envelopeMax_NUM/envelopeMax_DENOM;
                envelopeMin_RATIO = envelopeMin_NUM/envelopeMin_DENOM;
            }

            // std::cout << envelopeMax_NUM << " " << envelopeMax_DENOM << " " << envelopeMax_RATIO << std::endl;
            // std::cout << envelopeMin_NUM << " " << envelopeMin_DENOM << " " << envelopeMin_RATIO << std::endl;

            std::cout << ">>>"  << std::endl;
            std::cout << "bin #" << iBin << ": envelopeMax_RATIO - 1 = " << (envelopeMax_RATIO - 1.) << std::endl;
            std::cout << "bin #" << iBin << ": envelopeMin_RATIO - 1 = " << (envelopeMin_RATIO - 1.) << std::endl;
            
            if ( (envelopeMax_RATIO > 0.) && (envelopeMin_RATIO > 0.) ){
                h_Err_Scale_Up->SetBinContent(iBin, fabs(envelopeMax_RATIO - 1.));
                h_Err_Scale_Down->SetBinContent(iBin, fabs(envelopeMin_RATIO - 1.));
            }
            else{
                h_Err_Scale_Up->SetBinContent(iBin, 0.);
                h_Err_Scale_Down->SetBinContent(iBin, 0.);
            }

        }
    }

    // --------------------------------------------------------------------
    // Total Theoretical Uncertainty 

    // setup main error histograms
    // up variation
    TH1D *h_Err_Total_Up = (TH1D*) h_central_NUM->Clone(hNameNUM+"_TO_"+hNameDENOM+"_ErrTotal_Up");
    h_Err_Total_Up->Reset();
    h_Err_Total_Up->SetTitle(hNameNUM+"_TO_"+hNameDENOM+"_ErrTotal_Up");
    h_Err_Total_Up->GetYaxis()->SetRangeUser(0., 0.1);
    h_Err_Total_Up->GetYaxis()->SetTitle("Relative Uncertainty (Total), Up Variation");
    // down variation
    TH1D *h_Err_Total_Down = (TH1D*) h_central_NUM->Clone(hNameNUM+"_TO_"+hNameDENOM+"_ErrTotal_Down");
    h_Err_Total_Down->Reset();
    h_Err_Total_Down->SetTitle(hNameNUM+"_TO_"+hNameDENOM+"_ErrTotal_Down");
    h_Err_Total_Down->GetYaxis()->SetRangeUser(0., 0.1);
    h_Err_Total_Down->GetYaxis()->SetTitle("Relative Uncertainty (Total), Down Variation");

    if (h_Err_AlphaS_Up && h_Err_PDF_Up && h_Err_Scale_Up){
        std::cout << "\n>>> Total Theoretical Uncertainty:" << std::endl;

        // add all three sources of (relative) uncertainty in quadrature
        for (int iBin(1); iBin < h_central_NUM->GetNbinsX()+1; iBin++){

            double errAlphaS_Up = h_Err_AlphaS_Up->GetBinContent(iBin);
            double errScale_Up  = h_Err_Scale_Up->GetBinContent(iBin);
            double errPDF_Up    = h_Err_PDF_Up->GetBinContent(iBin);
            double errTotal_Up = sqrt( (errAlphaS_Up*errAlphaS_Up) + (errScale_Up*errScale_Up) + (errPDF_Up*errPDF_Up) );

            double errAlphaS_Down = h_Err_AlphaS_Down->GetBinContent(iBin);
            double errScale_Down  = h_Err_Scale_Down->GetBinContent(iBin);
            double errPDF_Down    = h_Err_PDF_Down->GetBinContent(iBin);
            double errTotal_Down = sqrt( (errAlphaS_Down*errAlphaS_Down) + (errScale_Down*errScale_Down) + (errPDF_Down*errPDF_Down) );
            
            std::cout << "bin #" << iBin << ": errTotal_Up   = " << errTotal_Up << std::endl;
            std::cout << "bin #" << iBin << ": errTotal_Down = " << errTotal_Down << std::endl;

            h_Err_Total_Up->SetBinContent(iBin, errTotal_Up);
            h_Err_Total_Down->SetBinContent(iBin, errTotal_Down);
        }
    }

    // ----------------------------------------------------------------
    // Write results out to output file

    file_output->cd();

    h_Err_AlphaS_Up->Write();
    h_Err_AlphaS_Down->Write();

    h_Err_PDF_Up->Write();
    h_Err_PDF_Down->Write();

    h_Err_Scale_Up->Write();
    h_Err_Scale_Down->Write();

    h_Err_Total_Up->Write();
    h_Err_Total_Down->Write();

    
    // ----------------------------------------------------------------
    std::cout << "\n\n----------------------------------------------------------------------------" << std::endl;
    std::cout << "\nClosing all files!" << std::endl;
    // -- central
    if (file_central->IsOpen()) file_central->Close();
    // -- alpha-s variations
    if (file_alphaS_1->IsOpen()) file_alphaS_1->Close();
    if (file_alphaS_2->IsOpen()) file_alphaS_2->Close();
    // -- PDF set variations
    if (file_PDF_1->IsOpen()) file_PDF_1->Close();
    if (file_PDF_2->IsOpen()) file_PDF_2->Close();
    // -- QCD scale variations
    if (file_scales_1->IsOpen()) file_scales_1->Close();
    if (file_scales_2->IsOpen()) file_scales_2->Close();
    if (file_scales_3->IsOpen()) file_scales_3->Close();
    if (file_scales_4->IsOpen()) file_scales_4->Close();
    if (file_scales_5->IsOpen()) file_scales_5->Close();
    if (file_scales_6->IsOpen()) file_scales_6->Close();
    // -- output
    if (file_output->IsOpen()) file_output->Close();
    // ----------------------------------------------------------------
    std::cout << "\nFinished!\n" << std::endl;
}
