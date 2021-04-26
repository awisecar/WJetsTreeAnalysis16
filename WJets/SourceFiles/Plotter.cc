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

#include "fileNames.h"
#include "getFilesAndHistograms.h"

using namespace std;

void Plotter(string leptonFlavor = "SMu", int year = 2016, int JetPtMin = 30,
    int doQCD = 0, bool doSSign = 0, bool doInvMassCut = 0, int METcut = 0 , int doBJets = -1, 
    int JetPtMax = 0, int ZEtaMin = -999999, int ZEtaMax = 999999, 
    bool doRoch = 0, bool doFlat = 0, bool doVarWidth = 1)
{
    std::string yearStr;
    std::stringstream yearSStr;
    yearSStr << year;
    yearStr = yearSStr.str();

    std::cout << "\nPlotting data/MC for year = " << yearStr << "!" << std::endl;
   
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
    if (leptonFlavor == "SMu" || leptonFlavor == "Muon" ) legendNames[0] = " #mu#nu ";

    //-----------------------------------------------------
    // Getting files ---
    int nFiles = 0;
    bool isDoubleLep(1);
    if (leptonFlavor == "SMu") {
        isDoubleLep = 0;

        nFiles = NFILESTTBARWJETS; // the nominal switch incl. QCD
         //nFiles = NFILESTTBARWJETS_NOQCD; // // turn off QCD, used for ttbar studies
        // nFiles = NFILESTTBARWJETS_NOQCD_NOTTBAR; // 
        //nFiles = NFILESTTBARWJETS_DATAQCDVV; // just run QCD and VV

    }
    TFile *file[nFiles];
    int countFiles = 0 ;
    for (unsigned short i = 0; i < nFiles; i++){
        int fileSelect = 0;

        if (!isDoubleLep) fileSelect = FilesTTbarWJets[i]; // the nominal switch incl. QCD
         //if (!isDoubleLep) fileSelect = FilesTTbarWJets_NoQCD[i]; // turn off QCD, used for ttbar studies
        // if (!isDoubleLep) fileSelect = FilesTTbarWJets_NoQCD_NoTTBar[i]; // andrew -- 2 sept 2019 -- turn off QCD for now
        // if (!isDoubleLep) fileSelect = FilesTTbarWJets_DataQCDVV[i]; // just run QCD and VV

        // FilesTTbarWJets will select the correct files for W+jets (defined in fileNames.h)
        string fileNameTemp =  ProcessInfo[fileSelect].filename ; 

        cout << "isDoubleLep = " << isDoubleLep << ", leptonFlavor = " << leptonFlavor << ", fileNameTemp = " << fileNameTemp << endl;
        if ((doQCD > 0 || doInvMassCut || doSSign ) && fileNameTemp.find("QCD") != string::npos) continue;
        //if (fileNameTemp.find("QCD") != string::npos) continue; // use this line if you do not want to plot QCD
        file[countFiles] = getFile(FILESDIRECTORY, leptonFlavor, energy, fileNameTemp, JetPtMin, JetPtMax, doFlat, doVarWidth, doQCD, doSSign, doInvMassCut, METcut, doBJets);

        if ( i == 0 ){
            //if  (leptonFlavor == "SMu" || leptonFlavor == "Muon" ) legendNames[0] = " #mu#nu ";
            if  (leptonFlavor == "SMu" || leptonFlavor == "Muon" ) legendNames[0] = "";
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
    else if (doBJets < 0 ) outputFileName += "_BVeto";
    else outputFileName += "";
    
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

    // NOTE: commenting out outputFile for now
    // TFile *outputFile = new TFile(outputFileRoot.c_str(), "RECREATE");
    // outputFile->cd();

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

    int nHistNoGen = 0;
    std::cout << "\nLooping over " << nHist << " histograms" << std::endl;
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
    nHist = nHistNoGen; 
    cout << "Number of histograms to plot: " << nHistNoGen << endl;

    ////////////////////////////////////////////////////////////////////////////////////////

    // initialize ttBar SF value --
    double SFttbar(1.);

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
	
    //looping over files
    for (int i = 0; i < nFiles; i++) {
        //cout << "File #" << i << ": " << legendNames[i] << endl;

        //Looping over histograms (nHistNoGen is all of the found histograms, without double counting the gen ones)
        for (int j = 0; j < nHistNoGen ; j++) {

            //std::cout << histoName[j] << std::endl;

            hist[i][j] = getHisto(file[i], histoName[j]);
            hist[i][j]->SetTitle(histoTitle[j].c_str());

            if (i == 0){
                //Don't need to do anything to the data, just plotting the points on top of the stacked MC
                hist[i][j]->SetMarkerStyle(20);
                hist[i][j]->SetMarkerColor(Colors[i]);
                hist[i][j]->SetLineColor(Colors[i]);
            }
            
			// -----------------------------------------------------------------------------------------
			
            //andrew -- comment out this next block if you want to turn ttBar SFs off
	        // // else if (i == 5){ // ttbar is file 5 when QCD is turned off (ttbar study) -- see fileNames.h
            else if (i == 6){ //ttbar is file 6 when running usual distributions with bveto == -1

               hist[i][j]->SetFillColor(Colors[i]);
               hist[i][j]->SetLineColor(Colors[i]);
               legend[j]->AddEntry(hist[i][j], legendNames[i].c_str(), "f");

               //factors applied for jet multiplicities of 2 to 8
               //here I need to make sure we apply SFs to only ttbar distribution, not to other MC distributions i.e. file[5]=ttbar (see fileNames.h)
               //note: I believe the scaling of the histograms is actually not saved back into the original file; 
               //the getFile function used to grab the source file only "reads" the rootfile, doesn't "recreate" it
               if (year == 2016){
                   // inclusive jet cut
                   if (histoName[j].find("Zinc2jet")!= string::npos)      SFttbar = ttBarSF_2016inc[0];
                   else if (histoName[j].find("Zinc3jet")!= string::npos) SFttbar = ttBarSF_2016inc[1];
                   else if (histoName[j].find("Zinc4jet")!= string::npos) SFttbar = ttBarSF_2016inc[2];
                   else if (histoName[j].find("Zinc5jet")!= string::npos) SFttbar = ttBarSF_2016inc[3];
                   else if (histoName[j].find("Zinc6jet")!= string::npos) SFttbar = ttBarSF_2016inc[4];
                   else if (histoName[j].find("Zinc7jet")!= string::npos) SFttbar = ttBarSF_2016inc[5];
                   else if (histoName[j].find("Zinc8jet")!= string::npos) SFttbar = ttBarSF_2016inc[6];
                   // exclusive jet cut
                   else if (histoName[j].find("Zexc2jet")!= string::npos) SFttbar = ttBarSF_2016exc[0];
                   else if (histoName[j].find("Zexc3jet")!= string::npos) SFttbar = ttBarSF_2016exc[1];
                   else if (histoName[j].find("Zexc4jet")!= string::npos) SFttbar = ttBarSF_2016exc[2];
                   else if (histoName[j].find("Zexc5jet")!= string::npos) SFttbar = ttBarSF_2016exc[3];
                   else if (histoName[j].find("Zexc6jet")!= string::npos) SFttbar = ttBarSF_2016exc[4];
                   else if (histoName[j].find("Zexc7jet")!= string::npos) SFttbar = ttBarSF_2016exc[5];
                   else if (histoName[j].find("Zexc8jet")!= string::npos) SFttbar = ttBarSF_2016exc[6];
                   // if not >= 2 jets cut, set xsec SF to 1
                   else SFttbar = 1.;
               }
               else if (year == 2017){
                   // inclusive jet cut
                   if (histoName[j].find("Zinc2jet")!= string::npos)      SFttbar = ttBarSF_2017inc[0];
                   else if (histoName[j].find("Zinc3jet")!= string::npos) SFttbar = ttBarSF_2017inc[1];
                   else if (histoName[j].find("Zinc4jet")!= string::npos) SFttbar = ttBarSF_2017inc[2];
                   else if (histoName[j].find("Zinc5jet")!= string::npos) SFttbar = ttBarSF_2017inc[3];
                   else if (histoName[j].find("Zinc6jet")!= string::npos) SFttbar = ttBarSF_2017inc[4];
                   else if (histoName[j].find("Zinc7jet")!= string::npos) SFttbar = ttBarSF_2017inc[5];
                   else if (histoName[j].find("Zinc8jet")!= string::npos) SFttbar = ttBarSF_2017inc[6];
                   // exclusive jet cut
                   else if (histoName[j].find("Zexc2jet")!= string::npos) SFttbar = ttBarSF_2017exc[0];
                   else if (histoName[j].find("Zexc3jet")!= string::npos) SFttbar = ttBarSF_2017exc[1];
                   else if (histoName[j].find("Zexc4jet")!= string::npos) SFttbar = ttBarSF_2017exc[2];
                   else if (histoName[j].find("Zexc5jet")!= string::npos) SFttbar = ttBarSF_2017exc[3];
                   else if (histoName[j].find("Zexc6jet")!= string::npos) SFttbar = ttBarSF_2017exc[4];
                   else if (histoName[j].find("Zexc7jet")!= string::npos) SFttbar = ttBarSF_2017exc[5];
                   else if (histoName[j].find("Zexc8jet")!= string::npos) SFttbar = ttBarSF_2017exc[6];
                   // if not >= 2 jets cut, set xsec SF to 1
                   else SFttbar = 1.;
               }
               else{
                   // inclusive jet cut
                   if (histoName[j].find("Zinc2jet")!= string::npos)      SFttbar = ttBarSF_2018inc[0];
                   else if (histoName[j].find("Zinc3jet")!= string::npos) SFttbar = ttBarSF_2018inc[1];
                   else if (histoName[j].find("Zinc4jet")!= string::npos) SFttbar = ttBarSF_2018inc[2];
                   else if (histoName[j].find("Zinc5jet")!= string::npos) SFttbar = ttBarSF_2018inc[3];
                   else if (histoName[j].find("Zinc6jet")!= string::npos) SFttbar = ttBarSF_2018inc[4];
                   else if (histoName[j].find("Zinc7jet")!= string::npos) SFttbar = ttBarSF_2018inc[5];
                   else if (histoName[j].find("Zinc8jet")!= string::npos) SFttbar = ttBarSF_2018inc[6];
                   // exclusive jet cut
                   else if (histoName[j].find("Zexc2jet")!= string::npos) SFttbar = ttBarSF_2018exc[0];
                   else if (histoName[j].find("Zexc3jet")!= string::npos) SFttbar = ttBarSF_2018exc[1];
                   else if (histoName[j].find("Zexc4jet")!= string::npos) SFttbar = ttBarSF_2018exc[2];
                   else if (histoName[j].find("Zexc5jet")!= string::npos) SFttbar = ttBarSF_2018exc[3];
                   else if (histoName[j].find("Zexc6jet")!= string::npos) SFttbar = ttBarSF_2018exc[4];
                   else if (histoName[j].find("Zexc7jet")!= string::npos) SFttbar = ttBarSF_2018exc[5];
                   else if (histoName[j].find("Zexc8jet")!= string::npos) SFttbar = ttBarSF_2018exc[6];
                   // if not >= 2 jets cut, set xsec SF to 1
                   else SFttbar = 1.;
               }

               //the two sections below scale the jet multiplicities
               //bin number 3 is 2 jet mult. (histogram bin number 1 is 0 jet bin)
               if (histoName[j].find("ZNGoodJets_Zinc") != string::npos){
                   if (year == 2016){
                       // bin contents
                       hist[i][j]->SetBinContent(3, hist[i][j]->GetBinContent(3)*ttBarSF_2016inc[0]);
                       hist[i][j]->SetBinContent(4, hist[i][j]->GetBinContent(4)*ttBarSF_2016inc[1]);
                       hist[i][j]->SetBinContent(5, hist[i][j]->GetBinContent(5)*ttBarSF_2016inc[2]);
                       hist[i][j]->SetBinContent(6, hist[i][j]->GetBinContent(6)*ttBarSF_2016inc[3]);
                       hist[i][j]->SetBinContent(7, hist[i][j]->GetBinContent(7)*ttBarSF_2016inc[4]);
                       hist[i][j]->SetBinContent(8, hist[i][j]->GetBinContent(8)*ttBarSF_2016inc[5]);
                       hist[i][j]->SetBinContent(9, hist[i][j]->GetBinContent(9)*ttBarSF_2016inc[6]);
                       // bin errors 
                       hist[i][j]->SetBinError(3, hist[i][j]->GetBinError(3)*ttBarSF_2016inc[0]);
                       hist[i][j]->SetBinError(4, hist[i][j]->GetBinError(4)*ttBarSF_2016inc[1]);
                       hist[i][j]->SetBinError(5, hist[i][j]->GetBinError(5)*ttBarSF_2016inc[2]);
                       hist[i][j]->SetBinError(6, hist[i][j]->GetBinError(6)*ttBarSF_2016inc[3]);
                       hist[i][j]->SetBinError(7, hist[i][j]->GetBinError(7)*ttBarSF_2016inc[4]);
                       hist[i][j]->SetBinError(8, hist[i][j]->GetBinError(8)*ttBarSF_2016inc[5]);
                       hist[i][j]->SetBinError(9, hist[i][j]->GetBinError(9)*ttBarSF_2016inc[6]);
                   }
                   else if (year == 2017){
                       // bin contents
                       hist[i][j]->SetBinContent(3, hist[i][j]->GetBinContent(3)*ttBarSF_2017inc[0]);
                       hist[i][j]->SetBinContent(4, hist[i][j]->GetBinContent(4)*ttBarSF_2017inc[1]);
                       hist[i][j]->SetBinContent(5, hist[i][j]->GetBinContent(5)*ttBarSF_2017inc[2]);
                       hist[i][j]->SetBinContent(6, hist[i][j]->GetBinContent(6)*ttBarSF_2017inc[3]);
                       hist[i][j]->SetBinContent(7, hist[i][j]->GetBinContent(7)*ttBarSF_2017inc[4]);
                       hist[i][j]->SetBinContent(8, hist[i][j]->GetBinContent(8)*ttBarSF_2017inc[5]);
                       hist[i][j]->SetBinContent(9, hist[i][j]->GetBinContent(9)*ttBarSF_2017inc[6]);
                       // bin errors 
                       hist[i][j]->SetBinError(3, hist[i][j]->GetBinError(3)*ttBarSF_2017inc[0]);
                       hist[i][j]->SetBinError(4, hist[i][j]->GetBinError(4)*ttBarSF_2017inc[1]);
                       hist[i][j]->SetBinError(5, hist[i][j]->GetBinError(5)*ttBarSF_2017inc[2]);
                       hist[i][j]->SetBinError(6, hist[i][j]->GetBinError(6)*ttBarSF_2017inc[3]);
                       hist[i][j]->SetBinError(7, hist[i][j]->GetBinError(7)*ttBarSF_2017inc[4]);
                       hist[i][j]->SetBinError(8, hist[i][j]->GetBinError(8)*ttBarSF_2017inc[5]);
                       hist[i][j]->SetBinError(9, hist[i][j]->GetBinError(9)*ttBarSF_2017inc[6]);
                   }
                   else {
                       // bin contents
                       hist[i][j]->SetBinContent(3, hist[i][j]->GetBinContent(3)*ttBarSF_2018inc[0]);
                       hist[i][j]->SetBinContent(4, hist[i][j]->GetBinContent(4)*ttBarSF_2018inc[1]);
                       hist[i][j]->SetBinContent(5, hist[i][j]->GetBinContent(5)*ttBarSF_2018inc[2]);
                       hist[i][j]->SetBinContent(6, hist[i][j]->GetBinContent(6)*ttBarSF_2018inc[3]);
                       hist[i][j]->SetBinContent(7, hist[i][j]->GetBinContent(7)*ttBarSF_2018inc[4]);
                       hist[i][j]->SetBinContent(8, hist[i][j]->GetBinContent(8)*ttBarSF_2018inc[5]);
                       hist[i][j]->SetBinContent(9, hist[i][j]->GetBinContent(9)*ttBarSF_2018inc[6]);
                       // bin errors 
                       hist[i][j]->SetBinError(3, hist[i][j]->GetBinError(3)*ttBarSF_2018inc[0]);
                       hist[i][j]->SetBinError(4, hist[i][j]->GetBinError(4)*ttBarSF_2018inc[1]);
                       hist[i][j]->SetBinError(5, hist[i][j]->GetBinError(5)*ttBarSF_2018inc[2]);
                       hist[i][j]->SetBinError(6, hist[i][j]->GetBinError(6)*ttBarSF_2018inc[3]);
                       hist[i][j]->SetBinError(7, hist[i][j]->GetBinError(7)*ttBarSF_2018inc[4]);
                       hist[i][j]->SetBinError(8, hist[i][j]->GetBinError(8)*ttBarSF_2018inc[5]);
                       hist[i][j]->SetBinError(9, hist[i][j]->GetBinError(9)*ttBarSF_2018inc[6]);
                   }
               }
               else if (histoName[j].find("ZNGoodJets_Zexc") != string::npos){
                   if (year == 2016){
                       // bin contents
                       hist[i][j]->SetBinContent(3, hist[i][j]->GetBinContent(3)*ttBarSF_2016exc[0]);
                       hist[i][j]->SetBinContent(4, hist[i][j]->GetBinContent(4)*ttBarSF_2016exc[1]);
                       hist[i][j]->SetBinContent(5, hist[i][j]->GetBinContent(5)*ttBarSF_2016exc[2]);
                       hist[i][j]->SetBinContent(6, hist[i][j]->GetBinContent(6)*ttBarSF_2016exc[3]);
                       hist[i][j]->SetBinContent(7, hist[i][j]->GetBinContent(7)*ttBarSF_2016exc[4]);
                       hist[i][j]->SetBinContent(8, hist[i][j]->GetBinContent(8)*ttBarSF_2016exc[5]);
                       hist[i][j]->SetBinContent(9, hist[i][j]->GetBinContent(9)*ttBarSF_2016exc[6]);
                       // bin errors 
                       hist[i][j]->SetBinError(3, hist[i][j]->GetBinError(3)*ttBarSF_2016exc[0]);
                       hist[i][j]->SetBinError(4, hist[i][j]->GetBinError(4)*ttBarSF_2016exc[1]);
                       hist[i][j]->SetBinError(5, hist[i][j]->GetBinError(5)*ttBarSF_2016exc[2]);
                       hist[i][j]->SetBinError(6, hist[i][j]->GetBinError(6)*ttBarSF_2016exc[3]);
                       hist[i][j]->SetBinError(7, hist[i][j]->GetBinError(7)*ttBarSF_2016exc[4]);
                       hist[i][j]->SetBinError(8, hist[i][j]->GetBinError(8)*ttBarSF_2016exc[5]);
                       hist[i][j]->SetBinError(9, hist[i][j]->GetBinError(9)*ttBarSF_2016exc[6]);
                   }
                   else if (year == 2017){
                       // bin contents
                       hist[i][j]->SetBinContent(3, hist[i][j]->GetBinContent(3)*ttBarSF_2017exc[0]);
                       hist[i][j]->SetBinContent(4, hist[i][j]->GetBinContent(4)*ttBarSF_2017exc[1]);
                       hist[i][j]->SetBinContent(5, hist[i][j]->GetBinContent(5)*ttBarSF_2017exc[2]);
                       hist[i][j]->SetBinContent(6, hist[i][j]->GetBinContent(6)*ttBarSF_2017exc[3]);
                       hist[i][j]->SetBinContent(7, hist[i][j]->GetBinContent(7)*ttBarSF_2017exc[4]);
                       hist[i][j]->SetBinContent(8, hist[i][j]->GetBinContent(8)*ttBarSF_2017exc[5]);
                       hist[i][j]->SetBinContent(9, hist[i][j]->GetBinContent(9)*ttBarSF_2017exc[6]);
                       // bin errors 
                       hist[i][j]->SetBinError(3, hist[i][j]->GetBinError(3)*ttBarSF_2017exc[0]);
                       hist[i][j]->SetBinError(4, hist[i][j]->GetBinError(4)*ttBarSF_2017exc[1]);
                       hist[i][j]->SetBinError(5, hist[i][j]->GetBinError(5)*ttBarSF_2017exc[2]);
                       hist[i][j]->SetBinError(6, hist[i][j]->GetBinError(6)*ttBarSF_2017exc[3]);
                       hist[i][j]->SetBinError(7, hist[i][j]->GetBinError(7)*ttBarSF_2017exc[4]);
                       hist[i][j]->SetBinError(8, hist[i][j]->GetBinError(8)*ttBarSF_2017exc[5]);
                       hist[i][j]->SetBinError(9, hist[i][j]->GetBinError(9)*ttBarSF_2017exc[6]);
                   }
                   else{
                       // bin contents
                       hist[i][j]->SetBinContent(3, hist[i][j]->GetBinContent(3)*ttBarSF_2018exc[0]);
                       hist[i][j]->SetBinContent(4, hist[i][j]->GetBinContent(4)*ttBarSF_2018exc[1]);
                       hist[i][j]->SetBinContent(5, hist[i][j]->GetBinContent(5)*ttBarSF_2018exc[2]);
                       hist[i][j]->SetBinContent(6, hist[i][j]->GetBinContent(6)*ttBarSF_2018exc[3]);
                       hist[i][j]->SetBinContent(7, hist[i][j]->GetBinContent(7)*ttBarSF_2018exc[4]);
                       hist[i][j]->SetBinContent(8, hist[i][j]->GetBinContent(8)*ttBarSF_2018exc[5]);
                       hist[i][j]->SetBinContent(9, hist[i][j]->GetBinContent(9)*ttBarSF_2018exc[6]);
                       // bin errors 
                       hist[i][j]->SetBinError(3, hist[i][j]->GetBinError(3)*ttBarSF_2018exc[0]);
                       hist[i][j]->SetBinError(4, hist[i][j]->GetBinError(4)*ttBarSF_2018exc[1]);
                       hist[i][j]->SetBinError(5, hist[i][j]->GetBinError(5)*ttBarSF_2018exc[2]);
                       hist[i][j]->SetBinError(6, hist[i][j]->GetBinError(6)*ttBarSF_2018exc[3]);
                       hist[i][j]->SetBinError(7, hist[i][j]->GetBinError(7)*ttBarSF_2018exc[4]);
                       hist[i][j]->SetBinError(8, hist[i][j]->GetBinError(8)*ttBarSF_2018exc[5]);
                       hist[i][j]->SetBinError(9, hist[i][j]->GetBinError(9)*ttBarSF_2018exc[6]);
                   }
               }
               else {
                   //it is this line where we scale the ttbar MC for each of the histograms
                   //if the name of the histogram does not contain one of the phrases above, it is simply scaled by 1
                   hist[i][j]->Scale(SFttbar);
               }

	         } // END SCALING OF TTBAR

            // -----------------------------------------------------------------------------------------

            else {
                //if it's any of the other MC, set fill/line colors and legends like normal
                hist[i][j]->SetFillColor(Colors[i]);
                hist[i][j]->SetLineColor(Colors[i]);
                legend[j]->AddEntry(hist[i][j], legendNames[i].c_str(), "f");
            }

        } // end loop over histograms
    } // end loop over files
    
    //Here is where the THStack histSumMC is filled/stacked ---
    //Loop starts with i=1 here because hist[0] refers to the data file (see fileNames.h)
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
    
    ////////////////////////////////////////////////////////////////////////////////////////


    cout << "\n >>> Summed all MC histograms to create stacks! " << endl;

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
        histSumMC[i]->Draw("HIST"); 

        if ( (histoName[i].find("JetEta") != string::npos) || (histoName[i].find("JetAK8Eta") != string::npos) ){
            hist[0][i]->GetXaxis()->SetRangeUser(-2.4,2.4);
            histSumMC[i]->GetXaxis()->SetRangeUser(-2.4,2.4);
        }

        if (histoName[i].find("ZNGoodJets_") != string::npos){
            // nominal
           //hist[0][i]->GetXaxis()->SetRangeUser(0,10);
           //histSumMC[i]->GetXaxis()->SetRangeUser(0,10);
           // up to 8 jets
           hist[0][i]->GetXaxis()->SetRangeUser(0,8);
           histSumMC[i]->GetXaxis()->SetRangeUser(0,8);
            // for ttbar study
            //  hist[0][i]->GetXaxis()->SetRangeUser(2,8);
            //  histSumMC[i]->GetXaxis()->SetRangeUser(2,8);
        }

        hist[0][i]->SetTitle("");
        histSumMC[i]->SetTitle(""); 
        histSumMC[i]->GetYaxis()->SetLabelSize(0.048); //0.04
        histSumMC[i]->GetYaxis()->SetLabelOffset(0.002); 
        histSumMC[i]->GetYaxis()->SetLabelFont(42); 
        histSumMC[i]->GetYaxis()->SetTitle("Event count"); 
        histSumMC[i]->GetYaxis()->SetTitleFont(42); 
        histSumMC[i]->GetYaxis()->SetTitleSize(0.051); //0.04
        histSumMC[i]->GetYaxis()->SetTitleOffset(1.09); //1.2

        histSumMC[i]->SetMinimum(20); // andrew - was the option before
        // histSumMC[i]->SetMinimum(2);

        histSumMC[i]->SetMaximum(1000*histSumMC[i]->GetMaximum()); 

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
        cmsPre[i]->DrawLatex(0.27,0.87, " Work in Progress"); //uncomment later on
       
        if (yearStr == "2016")      intLumi[i]->DrawLatex(0.71,0.955, "35.9 fb^{-1} (13 TeV)");
        else if (yearStr == "2017") intLumi[i]->DrawLatex(0.71,0.955, "41.5 fb^{-1} (13 TeV)");
        else if (yearStr == "2018") intLumi[i]->DrawLatex(0.71,0.955, "59.7 fb^{-1} (13 TeV)");
        else intLumi[i]->DrawLatex(0.71,0.955, "137.16 fb^{-1} (13 TeV)");

        if (histoName[i].find("inc0") == string::npos){
            if (histoName[i].find("AK8") != string::npos){
                ostringstream ptLegend;
                if ( (histoName[i].find("FirstJetAK8Pt_") != string::npos) || (histoName[i].find("SecondJetAK8Pt_") != string::npos) || (histoName[i].find("LepPtPlusLeadingJetAK8Pt_Zinc") != string::npos) ){
                    ptLegend << "p_{T}^{jet} > 30 GeV,  |y^{jet}| < 2.4";
                }
                else{
                    ptLegend << "p_{T}^{jet} > 50 GeV,  |y^{jet}| < 2.4";
                }
                jetCuts[i]->DrawLatex(0.17,0.73, ptLegend.str().c_str());
                jetAlgo[i]->DrawLatex(0.17,0.80, "anti-k_{T} jets,  R = 0.8");
                if (histoName[i].find("LepPtPlusLeadingJetAK8Pt_Zinc") != string::npos){
                    jetAlgo[i]->DrawLatex(0.17,0.66, "p_{T}^{#mu}  > 50 GeV");
                }
                pad1[i]->Draw();
            }
            else{
                ostringstream ptLegend;
                ptLegend << "p_{T}^{jet} > " << JetPtMin << " GeV,  |y^{jet}| < 2.4";
                //ptLegend << "p_{T}^{jet} > 100 GeV, |y^{jet}| < 2.4";  //uncomment for DR plot
                jetCuts[i]->DrawLatex(0.17,0.75, ptLegend.str().c_str());
                jetAlgo[i]->DrawLatex(0.17,0.80, "anti-k_{T} jets,  R = 0.4");
                //jetAlgo[i]->DrawLatex(0.17,0.70, "Leading jet p_{T} > 300 GeV");  //uncomment for DR plot
                pad1[i]->Draw();
            }
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

        hist[0][i]->GetYaxis()->SetRangeUser(0.51, 1.49); // andrew -- the usual y-axis settings for the ratio plot
        // hist[0][i]->GetYaxis()->SetRangeUser(0.01, 1.99);

        hist[0][i]->GetYaxis()->SetNdivisions(5,5,0);
        hist[0][i]->GetYaxis()->SetTitle("Prediction/Data");
        hist[0][i]->GetYaxis()->SetTitleFont(42);
        hist[0][i]->GetYaxis()->SetTitleSize(0.104); //0.1
        hist[0][i]->GetYaxis()->SetTitleOffset(0.48); //0.5
        hist[0][i]->GetYaxis()->CenterTitle();
        hist[0][i]->GetYaxis()->SetLabelSize(0.09);  //0.08
        hist[0][i]->GetYaxis()->SetLabelFont(42); 

        //andrew -- dividing data (hist[0]) by stacked MC (histSumMC)
        hist[0][i]->Divide((TH1D*) histSumMC[i]->GetStack()->Last());
        for (unsigned short j(1); j <= nBins; j++){
     
            //These lines flip the bin content to Sim/Data ---
            //reassign both bincontent and binerrors
            double content(hist[0][i]->GetBinContent(j));
            double contentError(hist[0][i]->GetBinError(j));
            if (content > 0){
                hist[0][i]->SetBinContent(j, 1./content);
                // hist[0][i]->SetBinError(j, contentError*pow(hist[0][i]->GetBinContent(j),2.));
                hist[0][i]->SetBinError(j, contentError/(content*content));
            }
        }

        hist[0][i]->DrawCopy("EP");
        canvas[i]->cd();

        string outputFilePDF = outputFileName + "/" + histoName[i] + ".pdf";
        canvas[i]->Print(outputFilePDF.c_str());
        // outputFile->cd();
        // canvas[i]->Write(); // write TCanvas to the outputFile

        // Next section prints out the same plots on a linear scale -----
    //    histSumMC[i]->SetMaximum(1.5*histSumMC[i]->GetMaximum());
    //    TCanvas *tmpCanvas = (TCanvas*) canvas[i]->Clone();
    //    tmpCanvas->cd();
    //    tmpCanvas->Draw();
    //    TPad *tmpPad = (TPad*) tmpCanvas->GetPrimitive("pad1");
    //    tmpPad->SetLogy(0);
    //    histoName[i] += "_Lin";
    //    tmpCanvas->SetTitle(histoName[i].c_str());
    //    tmpCanvas->SetName(histoName[i].c_str());
    //    string outputFileLinPDF = outputFileName + "/" + histoName[i] + ".pdf";
    //    tmpCanvas->Print(outputFileLinPDF.c_str());
        //outputFile->cd();
        //tmpCanvas->Write();
    }

    // outputFile->cd();
    // outputFile->Close();

    cout << "\nClosing all files..." << endl;
    for (unsigned short i(0); i < nFiles; i++) closeFile(file[i]);

    cout << "\nPlotting finished!" << endl;
}
