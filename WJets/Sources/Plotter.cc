// History
//---- 2015_05_17
// Add, as a commented line, the one to use if you do not want to plot QCD contribution in the Signal region.

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

void Plotter(string leptonFlavor = "Muons", int JetPtMin = 30,
        int doQCD = 0, bool doSSign = 0, bool doInvMassCut = 0, int MET = 0 , int doBJets = -1 , 
        int JetPtMax = 0, int ZEtaMin = -999999, int ZEtaMax = 999999, 
        bool doRoch = 0, bool doFlat = 0, bool doVarWidth = 1)
{
    string energy = getEnergy();
    energy = "13TeV";

    cout << endl << "Running the Plotter with the following files as input: " <<doQCD << "   " << MET << "   " << doBJets << endl;
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    ostringstream JetPtMinStr;  JetPtMinStr << JetPtMin;
    ostringstream JetPtMaxStr;  JetPtMaxStr << JetPtMax;
    ostringstream ZEtaMinStr;   ZEtaMinStr << abs(ZEtaMin); 
    ostringstream ZEtaMaxStr;   ZEtaMaxStr << abs(ZEtaMax); 
    ostringstream doQCDStr;     doQCDStr << doQCD ;
    ostringstream METStr;   METStr << MET ; 
    string year = "2016"; 
    //if (energy == "13TeV") year = "2015";
    if (energy == "13TeV") year = "2016";
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

    //    legendNames[0] += "Data";

    //int ColorsTTbar[] = {kBlack, kOrange,  kRed+1,  kViolet+5, kPink, kBlue, kCyan+1, kCyan+1, kCyan+1,kGreen};

    //-- get the files ------------------------------------
    int nFiles = NFILESDYJETS;
    bool isDoubleLep(1);
    if ( leptonFlavor == "SMuE" || leptonFlavor == "SMu" || leptonFlavor == "Muon" || leptonFlavor == "Electron") {
        isDoubleLep = 0;
        nFiles = NFILESTTBARWJETS; 
    }
    TFile *file[nFiles];
    int countFiles = 0 ;
    for (unsigned short i = 0; i < nFiles; i++){
        int fileSelect = FilesDYJets[i] ;
        if (!isDoubleLep) fileSelect = FilesTTbarWJets[i] ;
        if ( leptonFlavor == "SMuE"  ) fileSelect = FilesTTbar[i] ;

        //if (!isDoubleLep) fileSelect = FilesTTbarWJets[i] ;
        string fileNameTemp =  ProcessInfo[fileSelect].filename ; 
        cout << "Is double lepton:" << isDoubleLep << "   " << leptonFlavor <<"   " << fileNameTemp << endl;
        if ((doQCD > 0 || doInvMassCut || doSSign ) && fileNameTemp.find("QCD") != string::npos) continue;
        //if (fileNameTemp.find("QCD") != string::npos) continue; // use this line if you do not want to plot QCD
        file[countFiles] = getFile(FILESDIRECTORY, leptonFlavor, energy, fileNameTemp, JetPtMin, JetPtMax, doFlat, doVarWidth, doQCD, doSSign,    doInvMassCut, MET, doBJets);

        if ( i == 0 ){
            if (leptonFlavor == "Electrons") legendNames[0] = " ee ";
            else if  (leptonFlavor == "Muons") legendNames[0] = " #mu#mu ";
            else if  (leptonFlavor == "SMuE") legendNames[0] = " #mue ";
            //else if  (leptonFlavor == "SMu" || leptonFlavor == "Muon" ) legendNames[0] = " #mu#nu ";
            else if  (leptonFlavor == "SMu" || leptonFlavor == "Muon" ) legendNames[0] = "";
            else if  (leptonFlavor == "SE" || leptonFlavor == "Electron") legendNames[0] = " e#nu ";
            legendNames[0] += "Data";
        }
        else legendNames[countFiles] = ProcessInfo[fileSelect].legend; 
        Colors[countFiles] = ProcessInfo[fileSelect].color;    
        countFiles++;
    }
    nFiles = countFiles ;
    //-----------------------------------------------------

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

    if (doQCD>0) outputFileName += "_QCD" + doQCDStr.str();
    if ( MET > 0 ) outputFileName += "_MET"+METStr.str();

    //if (doBJets) nameStr << "_BJets";

    /*    if (doQCD > 0 ){
          if ( leptonFlavor == "SMuE") outputFileName +="_SS";
          else  outputFileName += "_QCD" + doQCDStr.str();
          }
          */
    outputFileName += "_SFInvers";
    if (doInvMassCut) outputFileName += "_InvMass";
    // create the directory if it doesn't exist
    string command = "mkdir -p " + outputFileName;
    system(command.c_str());
    string outputFileRoot = outputFileName + ".root";

    cout << "Output directory is: " << outputFileName << endl;
    cout << "Output root file is: " << outputFileRoot << endl;

    TFile *outputFile = new TFile(outputFileRoot.c_str(), "RECREATE");
    outputFile->cd();

    unsigned short nHist = file[0]->GetListOfKeys()->GetEntries();
    //vector<string> histoName(nHist);
    //vector<string> histoTitle(nHist);
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
        histoNameTemp = file[0]->GetListOfKeys()->At(i)->GetName();
        //looks at all of the histograms without double-counting the gen ones
        if (histoNameTemp.find("gen") != string::npos) continue;
        histTemp = (TH1D*) file[0]->Get(histoNameTemp.c_str());
        if (histTemp->GetEntries() < 1) continue;
        //only TH1's considered
        if (!histTemp->InheritsFrom(TH1D::Class())) continue;

        //        histoName[nHistNoGen] = file[0]->GetListOfKeys()->At(i)->GetName();
        //        histoTitle[nHistNoGen] = file[0]->GetListOfKeys()->At(i)->GetTitle();
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
    cout <<"Number of histograms:" << nHist << " and we plot :" << nHistNoGen << endl;
    //nHist=4;
    /*for (int i = 0; i < nFiles; i++) {
        cout << i <<"  "<<legendNames[i]  << "   "<<  endl;
        for (int j = 0; j < nHistNoGen ; j++) {
            // cout << i <<"  "<<legendNames[i]  << "   "<< j << "   " <<  histoName[j] << endl;
            hist[i][j] = getHisto(file[i], histoName[j]);
            hist[i][j]->SetTitle(histoTitle[j].c_str());
            if ( i == 0) {
                hist[i][j]->SetMarkerStyle(20);
                hist[i][j]->SetMarkerColor(Colors[i]);
                hist[i][j]->SetLineColor(Colors[i]);
            }
            else {
                hist[i][j]->SetFillColor(Colors[i]);
                hist[i][j]->SetLineColor(Colors[i]);
                legend[j]->AddEntry(hist[i][j], legendNames[i].c_str(), "f");
            }
        }
    }
    
    for (int i = 1; i < nFiles; i++) {
        for (int j = 0; j < nHistNoGen ; j++) {
            if (doBJets <= 0 ){
                histSumMC[j]->Add(hist[i][j]);
            }
            else {
                if (i == nFiles - 2) histSumMC[j]->Add(hist[nFiles - 1][j]);
                else if (i == nFiles - 1) histSumMC[j]->Add(hist[nFiles - 2][j]);
                else histSumMC[j]->Add(hist[i][j]);
            }
        }
    } */  //use this part if you do not want to rescale ttbar 
	

    ////////////////////this is where we start using ttbar rescaling option


    //andrew -- ttbar SFs before MET filters and B-tag SFs applied -- 6.12.2017
    double SFttbar(1);
	
    //looping over files
    for (int i = 0; i < nFiles; i++) {
        cout << i <<"  "<<legendNames[i]  << "   "<<  endl;
        //looping over histograms (nHistNoGen is all of the found histograms, without double counting the gen ones)
        for (int j = 0; j < nHistNoGen ; j++) {
            // cout << i <<"  "<<legendNames[i]  << "   "<< j << "   " <<  histoName[j] << endl;
            hist[i][j] = getHisto(file[i], histoName[j]);
            hist[i][j]->SetTitle(histoTitle[j].c_str());
			
            //need factors for exclusive and inclusive distributions
            //factors applied for jet multiplicities of 2 to 6
            //here I need to make sure we apply SFs to only ttbar distribution, not to other MC distributions i.e. file[5]=ttbar (see fileNames.h)
            //note: I believe the scaling of the histograms is actually not saved back into the original file; the getFile function used to grab the source file only "reads" the rootfile, doesn't "recreate" it
	    if (histoName[j].find("Zinc2jet")!= string::npos) SFttbar = 1.100549991;
	    else if (histoName[j].find("Zinc3jet")!= string::npos) SFttbar = 1.016080563;
	    else if (histoName[j].find("Zinc4jet")!= string::npos) SFttbar = 0.979844961;
            else if (histoName[j].find("Zinc5jet")!= string::npos) SFttbar = 0.941202003;
            else if (histoName[j].find("Zinc6jet")!= string::npos) SFttbar = 0.927091125;
            else if (histoName[j].find("Zexc2jet")!= string::npos) SFttbar = 1.461640581;
            else if (histoName[j].find("Zexc3jet")!= string::npos) SFttbar = 1.06483934;
            else if (histoName[j].find("Zexc4jet")!= string::npos) SFttbar = 1.006837223;
            else if (histoName[j].find("Zexc5jet")!= string::npos) SFttbar = 0.948498974;
            else if (histoName[j].find("Zexc6jet")!= string::npos) SFttbar = 0.92190308;
			
            if ( i == 0) {
                //don't need to do anything to the data, just plotting the points on top of the stacked MC
                hist[i][j]->SetMarkerStyle(20);
                hist[i][j]->SetMarkerColor(Colors[i]);
                hist[i][j]->SetLineColor(Colors[i]);
            }
            //andrew -- comment out this next block if you want to turn SFs off
	    //else if (i == 4){ // ttbar is file #4 when QCD is turned off (ttbar study) -- see fileNames.h
//	    else if (i == 5){ //ttbar is file #5 when running usual distributions with bveto == -1
//            //    //it is this line where we scale the ttbar MC for each of the histograms
//            //    //if the name of the histogram does not contain one of the phrases above, it is simply scaled by 1
//	    //    hist[i][j]->Scale(SFttbar);
//	
//                hist[i][j]->SetFillColor(Colors[i]);
//                hist[i][j]->SetLineColor(Colors[i]);
//                legend[j]->AddEntry(hist[i][j], legendNames[i].c_str(), "f");
//	
//                //the two sections below scale the jet multiplicities (each of the bins is a different number of jets)
//		if (histoName[j].find("ZNGoodJets_Zinc") != string::npos){
//                        //bin number 3 is 2 jet mult. (histogram bin number 0 is underflow)
//			hist[i][j]->SetBinContent(3, hist[i][j]->GetBinContent(3)*1.100549991);
//			hist[i][j]->SetBinContent(4, hist[i][j]->GetBinContent(4)*1.016080563);
//			hist[i][j]->SetBinContent(5, hist[i][j]->GetBinContent(5)*0.979844961);
//			hist[i][j]->SetBinContent(6, hist[i][j]->GetBinContent(6)*0.941202003);
//			hist[i][j]->SetBinContent(7, hist[i][j]->GetBinContent(7)*0.927091125);
//			hist[i][j]->SetBinContent(8, hist[i][j]->GetBinContent(8)*0.939210833);
//		
//			hist[i][j]->SetBinError(3, hist[i][j]->GetBinError(3)*1.100549991);
//			hist[i][j]->SetBinError(4, hist[i][j]->GetBinError(4)*1.016080563);
//			hist[i][j]->SetBinError(5, hist[i][j]->GetBinError(5)*0.979844961);
//			hist[i][j]->SetBinError(6, hist[i][j]->GetBinError(6)*0.941202003);
//			hist[i][j]->SetBinError(7, hist[i][j]->GetBinError(7)*0.927091125);
//			hist[i][j]->SetBinError(8, hist[i][j]->GetBinError(8)*0.939210833);
//		
//		}
//		else if (histoName[j].find("ZNGoodJets_Zexc") != string::npos){
//			hist[i][j]->SetBinContent(3, hist[i][j]->GetBinContent(3)*1.461640581);
//			hist[i][j]->SetBinContent(4, hist[i][j]->GetBinContent(4)*1.06483934);
//			hist[i][j]->SetBinContent(5, hist[i][j]->GetBinContent(5)*1.006837223);
//			hist[i][j]->SetBinContent(6, hist[i][j]->GetBinContent(6)*0.948498974);
//			hist[i][j]->SetBinContent(7, hist[i][j]->GetBinContent(7)*0.92190308);
//			hist[i][j]->SetBinContent(8, hist[i][j]->GetBinContent(8)*0.959601048);
//			
//			hist[i][j]->SetBinError(3, hist[i][j]->GetBinError(3)*1.461640581);
//			hist[i][j]->SetBinError(4, hist[i][j]->GetBinError(4)*1.06483934);
//			hist[i][j]->SetBinError(5, hist[i][j]->GetBinError(5)*1.006837223);
//			hist[i][j]->SetBinError(6, hist[i][j]->GetBinError(6)*0.948498974);
//			hist[i][j]->SetBinError(7, hist[i][j]->GetBinError(7)*0.92190308);
//			hist[i][j]->SetBinError(8, hist[i][j]->GetBinError(8)*0.959601048);
//		}
//                else {
//                        //it is this line where we scale the ttbar MC for each of the histograms
//                        //if the name of the histogram does not contain one of the phrases above, it is simply scaled by 1
//                        hist[i][j]->Scale(SFttbar);
//                }
//
//        	//hist[i][j]->SetFillColor(Colors[i]);
//        	//hist[i][j]->SetLineColor(Colors[i]);
//        	//legend[j]->AddEntry(hist[i][j], legendNames[i].c_str(), "f");
//	    }
            else {
                //if it's any of the other MC, set fill/line colors and legends like normal
                hist[i][j]->SetFillColor(Colors[i]);
                hist[i][j]->SetLineColor(Colors[i]);
                legend[j]->AddEntry(hist[i][j], legendNames[i].c_str(), "f");
            }
        }
    }
    
    //here is where the THStack histSumMC is filled/stacked
    //loop starts with i=1 here because hist[0] refers to the data file (see fileNames.h)
    for (int i = 1; i < nFiles; i++) {
        for (int j = 0; j < nHistNoGen ; j++) {
            if (doBJets <= 0 ){
                histSumMC[j]->Add(hist[i][j]);
            }
            else {
                //if no bveto, then the ttbar background changes in size
                //the following lines switch the order in which ttbar and dy+jets is stacked
                if (i == nFiles - 2) histSumMC[j]->Add(hist[nFiles - 1][j]);
                else if (i == nFiles - 1) histSumMC[j]->Add(hist[nFiles - 2][j]);
                else histSumMC[j]->Add(hist[i][j]);
            }
        }
    }
    
    ////////////////////this is where we end using ttbar rescaling option  


    cout << " added all histograms " << endl;

    for (unsigned short i(0); i < nHistNoGen; i++) {
        if (!file[0]->Get(histoName[i].c_str())->InheritsFrom(TH1D::Class())) continue;
        //cout << histoName[i] << endl;
        unsigned short nBins(hist[0][i]->GetNbinsX());
        legend[i]->AddEntry(hist[0][i], legendNames[0].c_str(), "ep");
        canvas[i] = new TCanvas(histoName[i].c_str(), histoName[i].c_str(), 700, 900);
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
        histSumMC[i]->Draw("HIST"); 
        if ( (leptonFlavor == "Muons" || leptonFlavor == "DMu" || leptonFlavor == "Electrons" ) && !doInvMassCut ) {  
            if (histoName[i].find("ZMass_Z") != string::npos){
                hist[0][i]->GetXaxis()->SetRangeUser(71,111);
                histSumMC[i]->GetXaxis()->SetRangeUser(71,111);
            }
            if (histoName[i].find("JetEta") != string::npos){
                hist[0][i]->GetXaxis()->SetRangeUser(-2.4,2.4);
                histSumMC[i]->GetXaxis()->SetRangeUser(-2.4,2.4);

            }
        }

         //   if (histoName[i].find("ZNGoodJets") != string::npos){
         //       hist[0][i]->GetXaxis()->SetRangeUser(0,6);
         //       histSumMC[i]->GetXaxis()->SetRangeUser(0,6);
         //   }
            
            //andrew -- use the below option if plotting the ttbar control region
            if (histoName[i].find("ZNGoodJets_") != string::npos){
                 hist[0][i]->GetXaxis()->SetRangeUser(2,6);
                 histSumMC[i]->GetXaxis()->SetRangeUser(2,6);
             }



        //andrew -- x-axis titles are set in HistoSet.cc (?)
        hist[0][i]->SetTitle("");
        histSumMC[i]->SetTitle(""); 
        histSumMC[i]->GetYaxis()->SetLabelSize(0.048); //0.04
        histSumMC[i]->GetYaxis()->SetLabelOffset(0.002); 
        histSumMC[i]->GetYaxis()->SetLabelFont(42); 
        histSumMC[i]->GetYaxis()->SetTitle("Number of events"); 
        histSumMC[i]->GetYaxis()->SetTitleFont(42); 
        histSumMC[i]->GetYaxis()->SetTitleSize(0.051); //0.04
        histSumMC[i]->GetYaxis()->SetTitleOffset(1.07); //1.2
        histSumMC[i]->SetMinimum(8);
        histSumMC[i]->SetMaximum(110*histSumMC[i]->GetMaximum()); 
        //andrew -- setting special y-axis bounds for select plots(?)
        if (histoName[i].find("AbsRapidity") != string::npos){
        histSumMC[i]->SetMaximum(2100*histSumMC[i]->GetMaximum()); 
            }
        if (histoName[i].find("dPhiLepJet") != string::npos){
        histSumMC[i]->SetMaximum(6300*histSumMC[i]->GetMaximum());  
            }
        if (histoName[i].find("LepPtPlusHT") != string::npos){
         histSumMC[i]->SetMaximum(5000*histSumMC[i]->GetMaximum());
            }

        /// first pad plots
        hist[0][i]->DrawCopy("e same");
        legend[i]->Draw();
        cmsColl[i]->DrawLatex(0.17,0.87, "CMS");
        cmsPre[i]->DrawLatex(0.27,0.87, "Preliminary"); //uncomment later on
       
        //if (energy == "8TeV") intLumi[i]->DrawLatex(0.20,0.83, "#int L dt = 19.6 fb^{-1},  #sqrt{s} = 8 TeV");
        //if (energy == "13TeV") intLumi[i]->DrawLatex(0.73,0.955, "2.2 fb^{-1} (13 TeV)");


        if ( histoName[i].find("inc0") == string::npos){
            ostringstream ptLegend;
            ptLegend << "p_{T}^{jet} > " << JetPtMin << " GeV,  |y^{jet}| < 2.4";
            //ptLegend << "p_{T}^{jet} > 100 GeV, |y^{jet}| < 2.4";  //uncomment for DR plot
            jetCuts[i]->DrawLatex(0.17,0.75, ptLegend.str().c_str());
            jetAlgo[i]->DrawLatex(0.17,0.80, "anti-k_{T} jets,  R = 0.4");
            //jetAlgo[i]->DrawLatex(0.17,0.70, "Leading jet p_{T} > 300 GeV");  //uncomment for DR plot
            pad1[i]->Draw();
        }  
   
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

        //andrew -- remember, the loop over i here is a loop over nHistNoGen
        //hist[0][i] means the ith histogram from the data file
        hist[0][i]->SetStats(0);
        hist[0][i]->SetTitle("");

        hist[0][i]->GetXaxis()->SetTickLength(0.03);
        hist[0][i]->GetXaxis()->SetTitleSize(0.126);  //0.11
        hist[0][i]->GetXaxis()->SetTitleOffset(1.11); //1.05
        hist[0][i]->GetXaxis()->SetTitleFont(42); 
        hist[0][i]->GetXaxis()->SetLabelSize(0.113);  //0.12
        hist[0][i]->GetXaxis()->SetLabelOffset(0.018);
        hist[0][i]->GetXaxis()->SetLabelFont(42); 

            if (histoName[i].find("ZNGoodJets") != string::npos){
                hist[0][i]->GetXaxis()->SetLabelSize(0.16);  //0.12
            }

        hist[0][i]->GetYaxis()->SetRangeUser(0.51,1.49);
        hist[0][i]->GetYaxis()->SetNdivisions(5,5,0);
        hist[0][i]->GetYaxis()->SetTitle("Simulation/Data");
        hist[0][i]->GetYaxis()->SetTitleFont(42);
        hist[0][i]->GetYaxis()->SetTitleSize(0.104); //0.1
        hist[0][i]->GetYaxis()->SetTitleOffset(0.48); //0.5
        hist[0][i]->GetYaxis()->CenterTitle();
        hist[0][i]->GetYaxis()->SetLabelSize(0.09);  //0.08
        hist[0][i]->GetYaxis()->SetLabelFont(42); 

        //andrew -- dividing data (hist[0]) by stacked MC (histSumMC)
        hist[0][i]->Divide((TH1D*) histSumMC[i]->GetStack()->Last());
        for (unsigned short j(1); j <= nBins; j++){

/*
            double content(hist[0][i]->GetBinContent(j));
            double error(hist[0][i]->GetBinError(j));

           double binW(hist[0][i]->GetBinWidth(j));

            if (content > 0){
                hist[0][i]->SetBinContent(j, content*1./(binW));
                //hist[0][i]->SetBinError(j, error*1./(binW));
            }
     */
            //andrew -- these lines flip the bin content to Sim/Data
            //no treatment of errors here??? - need to SetBinError?
            double content(hist[0][i]->GetBinContent(j));
            if (content > 0){
                hist[0][i]->SetBinContent(j, 1./content);
            }
        }

 

        hist[0][i]->DrawCopy("EP");
        canvas[i]->cd();

        //string outputFilePNG = outputFileName + "/" + histoName[i] + ".png";
        string outputFilePDF = outputFileName + "/" + histoName[i] + ".pdf";
        //canvas[i]->Print(outputFilePNG.c_str());
        canvas[i]->Print(outputFilePDF.c_str());
        outputFile->cd();
        canvas[i]->Write();

        histSumMC[i]->SetMaximum(1.5*histSumMC[i]->GetMaximum());
        TCanvas *tmpCanvas = (TCanvas*) canvas[i]->Clone();
        tmpCanvas->cd();
        tmpCanvas->Draw();
        TPad *tmpPad = (TPad*) tmpCanvas->GetPrimitive("pad1");
        tmpPad->SetLogy(0);
        histoName[i] += "_Lin";
        tmpCanvas->SetTitle(histoName[i].c_str());
        tmpCanvas->SetName(histoName[i].c_str());
        //string outputFileLinPNG = outputFileName + "/" + histoName[i] + ".png";
        string outputFileLinPDF = outputFileName + "/" + histoName[i] + ".pdf";
        //tmpCanvas->Print(outputFileLinPNG.c_str());
        tmpCanvas->Print(outputFileLinPDF.c_str());
        outputFile->cd();
        tmpCanvas->Write();
    }
    outputFile->cd();
    outputFile->Close();

    cout << "I m closing the files" << endl;
    //-- Close all the files ------------------------------
    for (unsigned short i(0); i < nFiles; i++) closeFile(file[i]);
    //-----------------------------------------------------
    cout << "Everything went fine -- Plotting done" << endl;

}
