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

void calcTheoryUncert_NLOFxFx(){
    welcomeMessage();
    // ----------------------------------------------------------------
    // Setup

    TString year = "2016";

    // base filename for the NLO FxFx inclusive sample
    TString filename_central = "HistoFiles_theoryUncert/SMu_13TeV_WJets_FxFx_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth";

    // Open nominal file and all theoretical uncertainty variations to get input histograms
    // -- nominal
    TFile *file_central = TFile::Open(filename_central+".root", "READ");
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
    TFile *file_output = TFile::Open("HistoFiles_theoryUncert/theoryUncertainties_nlofxfxinclusive_"+year+".root", "RECREATE");

    // ----------------------------------------------------------------
    // Loop over every histogram and do analysis
    // NOTE: Currently we are symmetrizing the up/down uncertainty variations into one uncertainty value

    // get total number of histograms in the file
    int nHist = file_central->GetListOfKeys()->GetEntries();
    std::cout << "\nnHist = " << nHist << "\n" << std::endl;

    // loop over all of the histograms
    for (int i(0); i < nHist; i++){

        // the "name" should be set as the same as the pointer in the code
        TString hName = file_central->GetListOfKeys()->At(i)->GetName();
        // do analysis only for the "gen" && "TUnfold" histograms
        if ( !hName.Contains("gen") || !hName.Contains("_TUnfold") ) continue;
        std::cout << "\n\n\n\n>>>>>>>>>>>>>>>>>>>>>>>> Analyzing: " << hName << std::endl;

        // grab central histogram first
        // NOTE: for error histograms, clone the main distribution to get the binning
        TH1D *h_central = (TH1D*) file_central->Get(hName);
        h_central->SetDirectory(0);

        // ----------------------------------------------------------------
        // alpha-s Uncertainty

        // set up the main error histogram
        TH1D *h_Err_AlphaS = (TH1D*) h_central->Clone(hName+"_ErrAlphaS");
        h_Err_AlphaS->Reset();
        h_Err_AlphaS->SetTitle(hName+"_ErrAlphaS");
        h_Err_AlphaS->GetYaxis()->SetRangeUser(0., 0.37);
        h_Err_AlphaS->GetYaxis()->SetTitle("Relative Uncertainty (#alpha_{s})");

        if (file_alphaS_1 && file_alphaS_2){
            std::cout << "\n>>> alpha-s Uncertainty:" << std::endl;

            // get alpha-s uncertainty variations
            TH1D *h_alphaS_1 = (TH1D*) file_alphaS_1->Get(hName);
            TH1D *h_alphaS_2 = (TH1D*) file_alphaS_2->Get(hName);

            for (int iBin(1); iBin < h_central->GetNbinsX()+1; iBin++){

                double contentCentral = h_central->GetBinContent(iBin);
                double contentVar1    = h_alphaS_1->GetBinContent(iBin);
                double contentVar2    = h_alphaS_2->GetBinContent(iBin);

                double fracErr1(0.), fracErr2(0.), fracErrAvg(0.);
                if (contentCentral > 0){
                    fracErr1   = (contentVar1 - contentCentral)/contentCentral;
                    fracErr2   = (contentVar2 - contentCentral)/contentCentral;
                    fracErrAvg = ( fabs(fracErr1) + fabs(fracErr2) )/2.0;
                }
                      
                std::cout << ">>>" << std::endl;
                std::cout << "bin #" << iBin << ": fracErr1 = "   << fracErr1 << std::endl;
                std::cout << "bin #" << iBin << ": fracErr2 = "   << fracErr2 << std::endl;
                std::cout << "bin #" << iBin << ": fracErrAvg = " << fracErrAvg << std::endl;
                h_Err_AlphaS->SetBinContent(iBin, fracErrAvg);
            }
        }

        // ----------------------------------------------------------------
        // PDF Set Uncertainty

        // setup main error histogram
        TH1D *h_Err_PDF = (TH1D*) h_central->Clone(hName+"_ErrPDF");
        h_Err_PDF->Reset();
        h_Err_PDF->SetTitle(hName+"_ErrPDF");
        h_Err_PDF->GetYaxis()->SetRangeUser(0., 0.37);
        h_Err_PDF->GetYaxis()->SetTitle("Relative Uncertainty (PDF)");

        if (file_PDF_1 && file_PDF_2){
            std::cout << "\n>>> PDF Uncertainty:" << std::endl;

            // get PDF uncertainty variations
            TH1D *h_PDF_1 = (TH1D*) file_PDF_1->Get(hName);
            TH1D *h_PDF_2 = (TH1D*) file_PDF_2->Get(hName);

            for (int iBin(1); iBin < h_central->GetNbinsX()+1; iBin++){

                double contentCentral = h_central->GetBinContent(iBin);
                double contentVar1    = h_PDF_1->GetBinContent(iBin);
                double contentVar2    = h_PDF_2->GetBinContent(iBin);

                double fracErr1(0.), fracErr2(0.), fracErrAvg(0.);
                if (contentCentral > 0){
                    fracErr1   = (contentVar1 - contentCentral)/contentCentral;
                    fracErr2   = (contentVar2 - contentCentral)/contentCentral;
                    fracErrAvg = ( fabs(fracErr1) + fabs(fracErr2) )/2.0;
                }
        
                std::cout << ">>>" << std::endl;
                std::cout << "bin #" << iBin << ": fracErr1 = "   << fracErr1 << std::endl;
                std::cout << "bin #" << iBin << ": fracErr2 = "   << fracErr2 << std::endl;
                std::cout << "bin #" << iBin << ": fracErrAvg = " << fracErrAvg << std::endl;
                h_Err_PDF->SetBinContent(iBin, fracErrAvg);
            }
        }

        // ----------------------------------------------------------------
        // QCD Scale Uncertainty

        // setup main error histogram
        TH1D *h_Err_Scale = (TH1D*) h_central->Clone(hName+"_ErrScale");
        h_Err_Scale->Reset();
        h_Err_Scale->SetTitle(hName+"_ErrScale");
        h_Err_Scale->GetYaxis()->SetRangeUser(0., 0.37);
        h_Err_Scale->GetYaxis()->SetTitle("Relative Uncertainty (QCD Scales)");

        if (file_scales_1 && file_scales_2 && file_scales_3 && file_scales_4 && file_scales_5 && file_scales_6){
            std::cout <<"\n>>> Scale Uncertainty:" << std::endl;

            // get scale uncertainty variations
            TH1D *h_scales_1 = (TH1D*) file_scales_1->Get(hName);
            TH1D *h_scales_2 = (TH1D*) file_scales_2->Get(hName);
            TH1D *h_scales_3 = (TH1D*) file_scales_3->Get(hName);
            TH1D *h_scales_4 = (TH1D*) file_scales_4->Get(hName);
            TH1D *h_scales_5 = (TH1D*) file_scales_5->Get(hName);
            TH1D *h_scales_6 = (TH1D*) file_scales_6->Get(hName);

            // out of all of the 6 scale variations, have to take the 
            // outer borders of the total envelope as the uncertainty
            for (int iBin(1); iBin < h_central->GetNbinsX()+1; iBin++){

                double contentCentral = h_central->GetBinContent(iBin);
                double contentVar1    = h_scales_1->GetBinContent(iBin);
                double contentVar2    = h_scales_2->GetBinContent(iBin);
                double contentVar3    = h_scales_3->GetBinContent(iBin);
                double contentVar4    = h_scales_4->GetBinContent(iBin);
                double contentVar5    = h_scales_5->GetBinContent(iBin);
                double contentVar6    = h_scales_6->GetBinContent(iBin);

                double envelopeMax(0.), envelopeMin(0.), fracErrAvg(0.);
                if (contentCentral > 0){

                    std::vector<double> scaleVar;
                    scaleVar.push_back(contentVar1/contentCentral);
                    scaleVar.push_back(contentVar2/contentCentral);
                    scaleVar.push_back(contentVar3/contentCentral);
                    scaleVar.push_back(contentVar4/contentCentral);
                    scaleVar.push_back(contentVar5/contentCentral);
                    scaleVar.push_back(contentVar6/contentCentral);

                    envelopeMax = *max_element(scaleVar.begin(), scaleVar.end());
                    envelopeMin = *min_element(scaleVar.begin(), scaleVar.end());
                    fracErrAvg  = ( fabs(envelopeMax - 1.) + fabs(envelopeMin - 1.) )/2.0;
                }
               
                std::cout << ">>>"  << std::endl;
                std::cout << "bin #" << iBin << ": envelopeMax - 1 = " << (envelopeMax - 1.) << std::endl;
                std::cout << "bin #" << iBin << ": envelopeMin - 1 = " << (envelopeMin - 1.) << std::endl;
                std::cout << "bin #" << iBin << ": fracErrAvg = " << fracErrAvg << std::endl;
                h_Err_Scale->SetBinContent(iBin, fracErrAvg);
            }
        }

        // --------------------------------------------------------------------
        // Total Theoretical Uncertainty 

        // setup main error histogram
        TH1D *h_Err_Total = (TH1D*) h_central->Clone(hName+"_ErrTotal");
        h_Err_Total->Reset();
        h_Err_Total->SetTitle(hName+"_ErrTotal");
        h_Err_Total->GetYaxis()->SetRangeUser(0., 0.37);
        h_Err_Total->GetYaxis()->SetTitle("Relative Uncertainty (Total)");

        if (h_Err_AlphaS && h_Err_PDF && h_Err_Scale){
            std::cout << "\n>>> Total Theoretical Uncertainty:" << std::endl;

            // add all three sources of (relative) uncertainty in quadrature
            for (int iBin(1); iBin < h_central->GetNbinsX()+1; iBin++){

                double errAlphaS = h_Err_AlphaS->GetBinContent(iBin);
                double errScale  = h_Err_Scale->GetBinContent(iBin);
                double errPDF    = h_Err_PDF->GetBinContent(iBin);
                double errTotal = sqrt( (errAlphaS*errAlphaS) + (errScale*errScale) + (errPDF*errPDF) );

                std::cout << "bin #" << iBin << ": errTotal = " << errTotal << std::endl;
                h_Err_Total->SetBinContent(iBin, errTotal);
            }
        }

        // ----------------------------------------------------------------
        // Write results out to output file

        file_output->cd();

        h_Err_AlphaS->Write();
        h_Err_PDF->Write();
        h_Err_Scale->Write();
        h_Err_Total->Write();

    }
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
