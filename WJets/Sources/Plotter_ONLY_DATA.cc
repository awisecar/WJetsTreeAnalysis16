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
#include "getFilesAndHistograms.h"

//--  Setting global variables --------------------------------------------------------------
#include "fileNames.h"
//-------------------------------------------------------------------------------------------

using namespace std;

void Plotter_ONLY_DATA(string leptonFlavor = "SMu", int year = 2017, int JetPtMin = 30,
    int doQCD = 0, bool doSSign = 0, bool doInvMassCut = 0, int METcut = 0 , int doBJets = -1 , 
    int JetPtMax = 0, int ZEtaMin = -999999, int ZEtaMax = 999999, 
    bool doRoch = 0, bool doFlat = 0, bool doVarWidth = 1)
{

    std::string yearStr;
    std::stringstream yearSStr;
    yearSStr << year;
    yearStr = yearSStr.str();

    std::cout << "\nPlotting data/MC for year = " << yearStr << "!" << std::endl;
   
    // string energy = getEnergy();
    string energy = "13TeV";

    cout << endl << "\n-----> Running the Plotter with the following options as input: " << endl;
    cout << "doQCD: " << doQCD << ", METcut: " << METcut << " doBJets: " << doBJets << endl;

    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    ostringstream JetPtMinStr;  JetPtMinStr << JetPtMin;
    ostringstream JetPtMaxStr;  JetPtMaxStr << JetPtMax;
    ostringstream ZEtaMinStr;   ZEtaMinStr << abs(ZEtaMin); 
    ostringstream ZEtaMaxStr;   ZEtaMaxStr << abs(ZEtaMax); 
    ostringstream doQCDStr;     doQCDStr << doQCD ;
    ostringstream METStr;   METStr << METcut ; 
    
    int Colors[] = {kBlack, kSpring+5, kOrange, kOrange-3, kRed+1, kPink-6, kViolet+5, kPink, kAzure+4, kBlue, kCyan+1, kCyan+1, kCyan+1}; 
    string legendNames[] = {
        " #mu#mu ", " ZZJets2L2Nu", " ZZJets4L", " ZZJets2L2Q", 
        " WZJets3LNu", " WZJets2L2Q", " WWJets2L2Nu", " Single Top", 
        " DYtautau", " TTJets", " DYJets MD", " DYJets Po", " DYJets Sh"
    };
    if (leptonFlavor == "Electrons") legendNames[0] = " ee ";
    else if  (leptonFlavor == "Muons") legendNames[0] = " #mu#mu ";
    else if  (leptonFlavor == "SMuE") legendNames[0] = " #mue ";
    else if  (leptonFlavor == "SMu" || leptonFlavor == "Muon" ) legendNames[0] = " #mu#nu ";
    else if  (leptonFlavor == "SE" || leptonFlavor == "Electron") legendNames[0] = " e#nu ";

    //-----------------------------------------------------
    // Getting files ---
    int nFiles = NFILESDYJETS;
    bool isDoubleLep(1);
    if ( leptonFlavor == "SMuE" || leptonFlavor == "SMu" || leptonFlavor == "Muon" || leptonFlavor == "Electron") {
        isDoubleLep = 0;
        // nFiles = NFILESTTBARWJETS; 
        nFiles = 1; //setting just to one file for only data
    }
    TFile *file[nFiles];
    int countFiles = 0 ;
    for (unsigned short i = 0; i < nFiles; i++){
        int fileSelect = FilesDYJets[i] ;
        // if (!isDoubleLep) fileSelect = FilesTTbarWJets[i];
        if (!isDoubleLep) fileSelect = 0; //setting just to one file for only data
        if (leptonFlavor == "SMuE") fileSelect = FilesTTbar[i] ;

        // FilesTTbarWJets will select the correct files for W+jets (defined in fileNames.h)
        string fileNameTemp =  ProcessInfo[fileSelect].filename ; 

        cout << "isDoubleLep = " << isDoubleLep << ", leptonFlavor = " << leptonFlavor << ", fileNameTemp = " << fileNameTemp << endl;
        if ((doQCD > 0 || doInvMassCut || doSSign ) && fileNameTemp.find("QCD") != string::npos) continue;
        //if (fileNameTemp.find("QCD") != string::npos) continue; // use this line if you do not want to plot QCD
        file[countFiles] = getFile(FILESDIRECTORY, leptonFlavor, energy, fileNameTemp, JetPtMin, JetPtMax, doFlat, doVarWidth, doQCD, doSSign, doInvMassCut, METcut, doBJets);

        if ( i == 0 ){
            if (leptonFlavor == "Electrons") legendNames[0] = " ee ";
            else if  (leptonFlavor == "Muons") legendNames[0] = " #mu#mu ";
            else if  (leptonFlavor == "SMuE") legendNames[0] = " #mue ";
            //else if  (leptonFlavor == "SMu" || leptonFlavor == "Muon" ) legendNames[0] = " #mu#nu ";
            else if  (leptonFlavor == "SMu" || leptonFlavor == "Muon" ) legendNames[0] = "";
            else if  (leptonFlavor == "SE" || leptonFlavor == "Electron") legendNames[0] = " e#nu ";
            legendNames[0] += "Data";
        }
        // for W+jets, if not Data, then grab the legend and color parts of the ProcessInfo struct for each file (in fileNames.h)
        else legendNames[countFiles] = ProcessInfo[fileSelect].legend; 
        Colors[countFiles] = ProcessInfo[fileSelect].color;    
        countFiles++;
    }
    nFiles = countFiles ;

    //-----------------------------------------------------
    // Creating output file ---
    string outputFileName = "PNGFiles/Comparison_" + leptonFlavor + "_" + energy + "_Data_All_MC_";

    outputFileName += "JetPtMin_" + JetPtMinStr.str();
    if (JetPtMax > JetPtMin) outputFileName += "_JetPtMax_" + JetPtMaxStr.str();
    if (ZEtaMin > -999999 && ZEtaMin <  0 ) outputFileName += "_ZEtaMin_m" + ZEtaMinStr.str();
    if (ZEtaMin > -999999 && ZEtaMin >= 0 ) outputFileName += "_ZEtaMin_" + ZEtaMinStr.str();
    if (ZEtaMax <  999999 && ZEtaMax >= 0 ) outputFileName += "_ZEtaMax_" + ZEtaMaxStr.str();
    if (ZEtaMax <  999999 && ZEtaMax <  0 ) outputFileName += "_ZEtaMax_m" + ZEtaMaxStr.str();
    if (doRoch) outputFileName += "_rochester";
    if (doFlat) outputFileName += "_Flat";
    if (doVarWidth) outputFileName += "_VarWidth";
    if (doInvMassCut) outputFileName +=  "_InvMass";
    if (doSSign )   outputFileName += "_SS";
    if (doBJets > 0 ) outputFileName += "_BJets";
    if (doBJets < 0 ) outputFileName += "_BVeto";
    if (doQCD > 0) outputFileName += "_QCD" + doQCDStr.str();
    if ( METcut > 0 ) outputFileName += "_MET"+METStr.str();
    outputFileName += "_SFInvers";
    if (doInvMassCut) outputFileName += "_InvMass";

    // create the directory if it doesn't exist
    string command = "mkdir -p " + outputFileName;
    system(command.c_str());
    string outputFileRoot = outputFileName + ".root";

    cout << "\nOutput directory is: " << outputFileName << endl;
    cout << "Output root file is: " << outputFileRoot << endl;

    TFile *outputFile = new TFile(outputFileRoot.c_str(), "RECREATE");
    outputFile->cd();

    //-----------------------------------------------------
    // Preparing to loop over all hists ---
    unsigned short nHist = file[0]->GetListOfKeys()->GetEntries();
    vector<string> histoName;
    vector<string> histoTitle;
    string histoNameTemp;
    TCanvas *canvas[nHist];
    TPad *pad1[nHist];
    TPad *pad2[nHist];
    TH1D *hist[25][nHist];
    TH1D *histTemp;
    THStack *histSumMC[nHist];
    TLegend *legend[nHist];
    TLatex *cmsColl[nHist];
    TLatex *cmsPre[nHist];
    TLatex *jetAlgo[nHist];
    TLatex *jetCuts[nHist];
    TLatex *intLumi[nHist];

    int nHistNoGen=0;
    for (unsigned short i(0); i < nHist; i++) {
        // Doublechecking that this histo is useable
        histoNameTemp = file[0]->GetListOfKeys()->At(i)->GetName();
        //looks at all of the histograms without double-counting the gen ones
        if (histoNameTemp.find("gen") != string::npos) continue;
        histTemp = (TH1D*) file[0]->Get(histoNameTemp.c_str());
        if (histTemp->GetEntries() < 1) continue;
        //only TH1's considered
        if (!histTemp->InheritsFrom(TH1D::Class())) continue;

        histoName.push_back(file[0]->GetListOfKeys()->At(i)->GetName());
        histoTitle.push_back(file[0]->GetListOfKeys()->At(i)->GetTitle());
        histSumMC[nHistNoGen] = new THStack(histoName[nHistNoGen].c_str(), histoTitle[nHistNoGen].c_str());

        double xLowLeg(0.66), xHighLeg(0.78);
        legend[nHistNoGen] = new TLegend(xLowLeg, 0.54, xHighLeg, 0.91);
        legend[nHistNoGen]->SetFillStyle(0);
        legend[nHistNoGen]->SetBorderSize(0);
        legend[nHistNoGen]->SetTextSize(0.040);
        legend[nHistNoGen]->SetTextFont(42);

        cmsColl[nHistNoGen] = new TLatex();
        cmsColl[nHistNoGen]->SetTextSize(0.054);
        cmsColl[nHistNoGen]->SetTextFont(61);
        cmsColl[nHistNoGen]->SetLineWidth(2);
        cmsColl[nHistNoGen]->SetTextColor(kBlack);
        cmsColl[nHistNoGen]->SetNDC();
        cmsColl[nHistNoGen]->SetTextAlign(11);

        cmsPre[nHistNoGen] = new TLatex();
        cmsPre[nHistNoGen]->SetTextSize(0.038);
        cmsPre[nHistNoGen]->SetTextFont(52);
        cmsPre[nHistNoGen]->SetLineWidth(1);
        cmsPre[nHistNoGen]->SetTextColor(kBlack);
        cmsPre[nHistNoGen]->SetNDC();
        cmsPre[nHistNoGen]->SetTextAlign(11);

        intLumi[nHistNoGen] = new TLatex();
        intLumi[nHistNoGen]->SetTextSize(0.040);
        intLumi[nHistNoGen]->SetTextFont(42);
        intLumi[nHistNoGen]->SetLineWidth(2);
        intLumi[nHistNoGen]->SetTextColor(kBlack);
        intLumi[nHistNoGen]->SetNDC();
        intLumi[nHistNoGen]->SetTextAlign(11);

        jetAlgo[nHistNoGen] = new TLatex();
        jetAlgo[nHistNoGen]->SetTextSize(0.040);
        jetAlgo[nHistNoGen]->SetTextFont(42);
        jetAlgo[nHistNoGen]->SetLineWidth(2);
        jetAlgo[nHistNoGen]->SetTextColor(kBlack);
        jetAlgo[nHistNoGen]->SetNDC();
        jetAlgo[nHistNoGen]->SetTextAlign(11);

        jetCuts[nHistNoGen] = new TLatex();
        jetCuts[nHistNoGen]->SetTextSize(0.040);
        jetCuts[nHistNoGen]->SetTextFont(42);
        jetCuts[nHistNoGen]->SetLineWidth(2);
        jetCuts[nHistNoGen]->SetTextColor(kBlack);
        jetCuts[nHistNoGen]->SetNDC();
        jetCuts[nHistNoGen]->SetTextAlign(11);

        nHistNoGen++;
    } 
    nHist=nHistNoGen; 
    cout <<"\nNumber of histograms to plot: " << nHistNoGen << endl;

    ////////////////////this is where we start using ttbar rescaling option


    //andrew -- ttbar SFs for 2016 data -- 8 April 2019
    double SFttbar(1);
	
    //looping over files
    // for (int i = 0; i < nFiles; i++) { // all input files
    for (int i = 0; i < 1; i++) { // just loop over data
        cout << "File #" << i << ": " << legendNames[i] << endl;

        //Looping over histograms (nHistNoGen is all of the found histograms, without double counting the gen ones)
        for (int j = 0; j < nHistNoGen ; j++) {
            hist[i][j] = getHisto(file[i], histoName[j]);
            hist[i][j]->SetTitle(histoTitle[j].c_str());
			
            //factors applied for jet multiplicities of 2 to 6
            //here I need to make sure we apply SFs to only ttbar distribution, not to other MC distributions i.e. file[5]=ttbar (see fileNames.h)
            //note: I believe the scaling of the histograms is actually not saved back into the original file; 
            //the getFile function used to grab the source file only "reads" the rootfile, doesn't "recreate" it

            // inclusive
            if (histoName[j].find("Zinc2jet")!= string::npos) SFttbar = 1.04904374;
            else if (histoName[j].find("Zinc3jet")!= string::npos) SFttbar = 0.99624372;
            else if (histoName[j].find("Zinc4jet")!= string::npos) SFttbar = 0.97243212;
            else if (histoName[j].find("Zinc5jet")!= string::npos) SFttbar = 0.94750793;
            else if (histoName[j].find("Zinc6jet")!= string::npos) SFttbar = 0.92506162;
            // exclusive
            else if (histoName[j].find("Zexc2jet")!= string::npos) SFttbar = 1.28532061;
            else if (histoName[j].find("Zexc3jet")!= string::npos) SFttbar = 1.02880218;
            else if (histoName[j].find("Zexc4jet")!= string::npos) SFttbar = 0.98994466;
            else if (histoName[j].find("Zexc5jet")!= string::npos) SFttbar = 0.95915763;
            else if (histoName[j].find("Zexc6jet")!= string::npos) SFttbar = 0.93695869;
			
            if ( i == 0) {
                //Don't need to do anything to the data, just plotting the points on top of the stacked MC
                hist[i][j]->SetMarkerStyle(20);
                hist[i][j]->SetMarkerColor(Colors[i]);
                hist[i][j]->SetLineColor(Colors[i]);
            }

            //andrew -- comment out this next block if you want to turn SFs off
	        //else if (i == 4){ // ttbar is file #4 when QCD is turned off (ttbar study) -- see fileNames.h
	        // else if (i == 5){ //ttbar is file #5 when running usual distributions with bveto == -1
            //     hist[i][j]->SetFillColor(Colors[i]);
            //     hist[i][j]->SetLineColor(Colors[i]);
            //     legend[j]->AddEntry(hist[i][j], legendNames[i].c_str(), "f");
	
            //     //the two sections below scale the jet multiplicities
            //     if (histoName[j].find("ZNGoodJets_Zinc") != string::npos){
            //                     //bin number 3 is 2 jet mult. (histogram bin number 1 is 0 jet bin)
            //         hist[i][j]->SetBinContent(3, hist[i][j]->GetBinContent(3)*1.04904374);
            //         hist[i][j]->SetBinContent(4, hist[i][j]->GetBinContent(4)*0.99624372);
            //         hist[i][j]->SetBinContent(5, hist[i][j]->GetBinContent(5)*0.97243212);
            //         hist[i][j]->SetBinContent(6, hist[i][j]->GetBinContent(6)*0.94750793);
            //         hist[i][j]->SetBinContent(7, hist[i][j]->GetBinContent(7)*0.92506162);
                
            //         hist[i][j]->SetBinError(3, hist[i][j]->GetBinError(3)*1.04904374);
            //         hist[i][j]->SetBinError(4, hist[i][j]->GetBinError(4)*0.99624372);
            //         hist[i][j]->SetBinError(5, hist[i][j]->GetBinError(5)*0.97243212);
            //         hist[i][j]->SetBinError(6, hist[i][j]->GetBinError(6)*0.94750793);
            //         hist[i][j]->SetBinError(7, hist[i][j]->GetBinError(7)*0.92506162);
            //     }
            //     else if (histoName[j].find("ZNGoodJets_Zexc") != string::npos){
            //         hist[i][j]->SetBinContent(3, hist[i][j]->GetBinContent(3)*1.28532061);
            //         hist[i][j]->SetBinContent(4, hist[i][j]->GetBinContent(4)*1.02880218);
            //         hist[i][j]->SetBinContent(5, hist[i][j]->GetBinContent(5)*0.98994466);
            //         hist[i][j]->SetBinContent(6, hist[i][j]->GetBinContent(6)*0.95915763);
            //         hist[i][j]->SetBinContent(7, hist[i][j]->GetBinContent(7)*0.93695869);
                    
            //         hist[i][j]->SetBinError(3, hist[i][j]->GetBinError(3)*1.28532061);
            //         hist[i][j]->SetBinError(4, hist[i][j]->GetBinError(4)*1.02880218);
            //         hist[i][j]->SetBinError(5, hist[i][j]->GetBinError(5)*0.98994466);
            //         hist[i][j]->SetBinError(6, hist[i][j]->GetBinError(6)*0.95915763);
            //         hist[i][j]->SetBinError(7, hist[i][j]->GetBinError(7)*0.93695869);
            //     }
            //     else {
            //         //it is this line where we scale the ttbar MC for each of the histograms
            //         //if the name of the histogram does not contain one of the phrases above, it is simply scaled by 1
            //         hist[i][j]->Scale(SFttbar);
            //     }

	        // } // END SCALING OF TTBAR

            else {
                //if it's any of the other MC, set fill/line colors and legends like normal
                hist[i][j]->SetFillColor(Colors[i]);
                hist[i][j]->SetLineColor(Colors[i]);
                legend[j]->AddEntry(hist[i][j], legendNames[i].c_str(), "f");
            }

        }
    }
    
    ////////////////////this is where we end using ttbar rescaling option  

    for (unsigned short i(0); i < nHistNoGen; i++) {

        if (!file[0]->Get(histoName[i].c_str())->InheritsFrom(TH1D::Class())) continue;
        unsigned short nBins(hist[0][i]->GetNbinsX());
        legend[i]->AddEntry(hist[0][i], legendNames[0].c_str(), "ep");
        canvas[i] = new TCanvas(histoName[i].c_str(), histoName[i].c_str(), 700, 900);

        // PAD 1 ------------------------------------------------------
        pad1[i] = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
        pad1[i]->SetTopMargin(0.055);
        pad1[i]->SetBottomMargin(0.);
        pad1[i]->SetRightMargin(0.04);
        pad1[i]->SetLeftMargin(0.12);
        pad1[i]->SetTicks();
        pad1[i]->SetLogy();
        pad1[i]->Draw();
        pad1[i]->cd();

        // Need to draw MC Stack first other wise
        // cannot access Xaxis !!!
        // histSumMC[i]->Draw("HIST"); 
        hist[0][i]->Draw("e"); // Drawing here instead!!!!

        if (histoName[i].find("ZNGoodJets") != string::npos){
            hist[0][i]->GetXaxis()->SetRangeUser(0,6);
            // histSumMC[i]->GetXaxis()->SetRangeUser(0,6);
        }
            
        //andrew -- use the below option if plotting the ttbar control region
        // if (histoName[i].find("ZNGoodJets_") != string::npos){
        //     hist[0][i]->GetXaxis()->SetRangeUser(2,6);
        //     histSumMC[i]->GetXaxis()->SetRangeUser(2,6);
        // }

        hist[0][i]->SetTitle("");
        hist[0][i]->SetTitle(""); 
        hist[0][i]->GetYaxis()->SetLabelSize(0.048); //0.04
        hist[0][i]->GetYaxis()->SetLabelOffset(0.002); 
        hist[0][i]->GetYaxis()->SetLabelFont(42); 
        hist[0][i]->GetYaxis()->SetTitle("Event count"); 
        hist[0][i]->GetYaxis()->SetTitleFont(42); 
        hist[0][i]->GetYaxis()->SetTitleSize(0.051); //0.04
        hist[0][i]->GetYaxis()->SetTitleOffset(1.07); //1.2
        // hist[0][i]->SetMinimum(8);
        hist[0][i]->SetMinimum(20);
        hist[0][i]->SetMaximum(110*hist[0][i]->GetMaximum()); 
        if (histoName[i].find("AbsRapidity") != string::npos){
        hist[0][i]->SetMaximum(2100*hist[0][i]->GetMaximum()); 
        }
        if (histoName[i].find("dPhiLepJet") != string::npos){
        hist[0][i]->SetMaximum(6300*hist[0][i]->GetMaximum());  
        }
        if (histoName[i].find("LepPtPlusHT") != string::npos){
         hist[0][i]->SetMaximum(5000*hist[0][i]->GetMaximum());
        }

        /// first pad plots
        // hist[0][i]->DrawCopy("e same");
        legend[i]->Draw();
        cmsColl[i]->DrawLatex(0.17,0.87, "CMS");
        cmsPre[i]->DrawLatex(0.27,0.87, " Work in Progress"); //uncomment later on
       
        if (yearStr == "2016") intLumi[i]->DrawLatex(0.71,0.955, "35.9 fb^{-1} (13 TeV)");
        else if (yearStr == "2017") intLumi[i]->DrawLatex(0.71,0.955, "41.5 fb^{-1} (13 TeV)");
        else intLumi[i]->DrawLatex(0.71,0.955, "? fb^{-1} (13 TeV)");


        if ( histoName[i].find("inc0") == string::npos){
            ostringstream ptLegend;
            ptLegend << "p_{T}^{jet} > " << JetPtMin << " GeV,  |y^{jet}| < 2.4";
            //ptLegend << "p_{T}^{jet} > 100 GeV, |y^{jet}| < 2.4";  //uncomment for DR plot
            jetCuts[i]->DrawLatex(0.17,0.75, ptLegend.str().c_str());
            jetAlgo[i]->DrawLatex(0.17,0.80, "anti-k_{T} jets,  R = 0.4");
            //jetAlgo[i]->DrawLatex(0.17,0.70, "Leading jet p_{T} > 300 GeV");  //uncomment for DR plot
            pad1[i]->Draw();
        }  
   
        // PAD 2 ------------------------------------------------------
        //Drawing the data-to-MC ratio plot on the bottom of the canvas
        canvas[i]->cd();
        pad2[i] = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
        pad2[i]->SetTopMargin(0.);
        pad2[i]->SetBottomMargin(0.32);  //0.26
        pad2[i]->SetRightMargin(0.04);
        pad2[i]->SetLeftMargin(0.12);
        pad2[i]->SetGridy();
        pad2[i]->SetTicks();
        pad2[i]->Draw();
        pad2[i]->cd();

        TH1D* hDataTemp = (TH1D *) hist[0][i]->Clone();

        hDataTemp->SetStats(0);
        hDataTemp->SetTitle("");

        hDataTemp->GetXaxis()->SetTickLength(0.03);
        hDataTemp->GetXaxis()->SetTitleSize(0.126);  //0.11
        hDataTemp->GetXaxis()->SetTitleOffset(1.11); //1.05
        hDataTemp->GetXaxis()->SetTitleFont(42); 
        hDataTemp->GetXaxis()->SetLabelSize(0.113);  //0.12
        hDataTemp->GetXaxis()->SetLabelOffset(0.018);
        hDataTemp->GetXaxis()->SetLabelFont(42); 

        if (histoName[i].find("ZNGoodJets") != string::npos){
            hDataTemp->GetXaxis()->SetLabelSize(0.16);  //0.12
        }

        hDataTemp->GetYaxis()->SetRangeUser(0.51,1.49);
        hDataTemp->GetYaxis()->SetNdivisions(5,5,0);
        hDataTemp->GetYaxis()->SetTitle("Simulation/Data");
        hDataTemp->GetYaxis()->SetTitleFont(42);
        hDataTemp->GetYaxis()->SetTitleSize(0.104); //0.1
        hDataTemp->GetYaxis()->SetTitleOffset(0.48); //0.5
        hDataTemp->GetYaxis()->CenterTitle();
        hDataTemp->GetYaxis()->SetLabelSize(0.09);  //0.08
        hDataTemp->GetYaxis()->SetLabelFont(42); 

        for (unsigned short j(1); j <= nBins; j++){
            hDataTemp->SetBinContent(j, 1.);
            hDataTemp->SetBinError(j, 0.);
        }

        hDataTemp->SetMarkerSize(0); 

        hDataTemp->Draw("EP");

        canvas[i]->cd();
        string outputFilePDF = outputFileName + "/" + histoName[i] + ".pdf";
        canvas[i]->Print(outputFilePDF.c_str());
        outputFile->cd();
        canvas[i]->Write();

    }
    outputFile->cd();
    outputFile->Close();

    cout << "\nClosing all files..." << endl;
    for (unsigned short i(0); i < nFiles; i++) closeFile(file[i]);
    cout << "\nPlotting finished!" << endl;
}