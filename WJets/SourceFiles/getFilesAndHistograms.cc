#include <iostream>
#include <sstream>
#include <TFile.h>
#include "fileNames.h"
#include "getFilesAndHistograms.h"
#include <string.h>
#include <algorithm>
using namespace std;

//------------------------------------------------------------
// getEnergy() returns a string, either "7TeV" or "13TeV"
// according to the name of the directory from which the 
// code is being executed.
//------------------------------------------------------------
string getEnergy()
{
    string energy = "";
    ostringstream fileBeingProcessed; fileBeingProcessed << __FILE__;
    if (fileBeingProcessed.str().find("Analysis2012") != string::npos) {
        energy = "13TeV";
    }
    else if (fileBeingProcessed.str().find("Analysis2011") != string::npos) {
        energy = "7TeV";
    }
    else 
    {
        std::cout << "WARNING ! Impossible to retrieve te energy from the current location !" << std::endl;
        energy = "Unknown";
    }
    fileBeingProcessed.str("");

    return energy;
}
//------------------------------------------------------------

TFile* getFile(string histoFilesDirectory, string leptonFlavor, string energy, string Name, int JetPtMin, int JetPtMax, bool doFlat, bool doVarWidth, int doQCD, bool doSSign, bool doInvMassCut, int MET, int doBJets, string closureTest, string syst, bool dodR, bool useUnfoldingFiles)
{
    
    string fileName = histoFilesDirectory; // string to contain the name of the file

    //--- make sure leptonFlavor is short version ---
    if (leptonFlavor == "Muons" || leptonFlavor == "DMu_") leptonFlavor = "DMu";
    else if (leptonFlavor == "Electrons" || leptonFlavor == "DE_") leptonFlavor = "DE";
    else if (leptonFlavor == "Electron" || leptonFlavor == "SE_") leptonFlavor = "SE";
    else if (leptonFlavor == "Muon" || leptonFlavor == "SMu_") leptonFlavor = "SMu";
    else if (leptonFlavor == "MuonElectron" || leptonFlavor == "SMuE_") leptonFlavor = "SMuE";

    fileName += leptonFlavor + "_" + energy + "_" + Name; // update fileName with lepton information
    //-----------------------------------------------


    //--- make string from numbers ---
    ostringstream JetPtMinStr; JetPtMinStr << JetPtMin;
    ostringstream JetPtMaxStr; JetPtMaxStr << JetPtMax;
    ostringstream doQCDStr; doQCDStr << doQCD;
    ostringstream METStr; METStr << MET;
    //--------------------------------

    //--- deal with efficiency correction applied or not ---
    string effiCorr = "1", trigCorr = "0";
    if (Name.find("Data") == 0 || energy == "13TeV") trigCorr = "1"; // trigger correction is applied to data and MC at 13TeV but only to data at 7TeV 
    if (useUnfoldingFiles) { // for cross-section measurement: correct data for efficiencies difference wrt MC
        if (Name.find("Data") == 0) effiCorr = "1";
        else effiCorr = "0";
    }
    else { // for control plots: correct MC for efficiencies difference wrt Data
        if (Name.find("Data") == 0) effiCorr = "0";
        else effiCorr = "1";
    }

    //--- special case for the generator comparison ---
    if (Name.find("Powheg") != string::npos || Name.find("Sherpa") != string::npos) {
        trigCorr = "0"; 
        effiCorr = "0";
    }
    //-------------------------------------------------

    fileName += "_EffiCorr_" + effiCorr + "_TrigCorr_" + trigCorr; // update fileName with efficiency correction information
    //------------------------------------------------------

    //--- update fileName for a bunch of other things ---
    fileName += "_Syst_" + syst + "_JetPtMin_" + JetPtMinStr.str();
    if (JetPtMax != 0 && JetPtMax > JetPtMin) fileName += "_JetPtMax_" + JetPtMaxStr.str();
    if (doFlat && Name.find("Data") == string::npos) fileName += "_Flat";
    if (closureTest != "") fileName += closureTest;
    if (doVarWidth) fileName += "_VarWidth";
    if (dodR) fileName += "_dR5";
    if (doInvMassCut) fileName += "_InvMass";
    if (doSSign) fileName += "_SS";
    if (doBJets > 0) fileName += "_BJets";
    if (doBJets < 0) fileName += "_BVeto";
    if (doQCD > 0) fileName += "_QCD" + doQCDStr.str();
    if (MET > 0) fileName += "_MET" + METStr.str();
    //---------------------------------------------------

    //--- fileName is complete: just add the extension and open it ---
    fileName += ".root";
    TFile *File = new TFile(fileName.c_str(), "READ");
    std::cout << "Opening: " << fileName << "   --->   Opened? " << File->IsOpen() << std::endl;
    return File;
    //----------------------------------------------------------------
}

void getFiles(string histoFilesDirectory, TFile *Files[], string leptonFlavor, string energy, string Name, int JetPtMin, int JetPtMax, bool doFlat, bool doVarWidth, int doQCD, bool doSSign, bool doInvMassCut, int MET, int doBJets, bool useUnfoldingFiles)
{

    //--- make sure leptonFlavor is short version ---
    if (leptonFlavor == "Muons" || leptonFlavor == "DMu_") leptonFlavor = "DMu";
    else if (leptonFlavor == "Electrons" || leptonFlavor == "DE_") leptonFlavor = "DE";
    else if (leptonFlavor == "Electron" || leptonFlavor == "SE_") leptonFlavor = "SE";
    else if (leptonFlavor == "Muon" || leptonFlavor == "SMu_") leptonFlavor = "SMu";
    else if (leptonFlavor == "MuonElectron" || leptonFlavor == "SMuE_") leptonFlavor = "SMuE";
    //-----------------------------------------------

    //--- set the double lepton flag ---
    bool isDoubleLep(0);
    if (leptonFlavor == "DMu" || leptonFlavor == "DE") isDoubleLep = 1;
    //----------------------------------

    vector<string> Syst;
    if ((Name.find("Data") != string::npos) && (Name.find("QCD") == string::npos)) { // for data we have:
        Syst.push_back("0");                 //   0: central
        Syst.push_back("2_Up");              //   2 up: JES up
        Syst.push_back("2_Down");            //   2 down: JES down
        if (leptonFlavor == "DE") {          // additionaly for electron:
            Syst.push_back("5_Up");          //   5 up: scale factor up
            Syst.push_back("5_Down");        //   5 down: scale factor down
        }
    }
    else if (Name.find("UNFOLDING") != string::npos && ((isDoubleLep && Name.find("DYJets") != string::npos) || (!isDoubleLep && Name.find("WJets") != string::npos))) {
        // for DYJets in case of Z+Jets or for WJets in case of W+Jets analysis we have:
        Syst.push_back("0");         // 0: central
        Syst.push_back("1_Up");      // 1 up: PU up
        Syst.push_back("1_Down");    // 1 down: PU down
        Syst.push_back("4_Up");      // 4 up: JER up
    }
    else { // for background we have
        Syst.push_back("0");         // 0: central
        Syst.push_back("1_Up");      // 1 up: PU up
        Syst.push_back("1_Down");    // 1 down: PU down
        Syst.push_back("3_Up");      // 3 up: XSec up
        Syst.push_back("3_Down");    // 3 down: Xsec down
    };

    //--- determnie how many files we have and open them all ---
    int nSyst(Syst.size());
    for (int i(0); i < nSyst; i++) {
        Files[i] = getFile(histoFilesDirectory, leptonFlavor, energy, Name, JetPtMin, JetPtMax, doFlat, doVarWidth, doQCD, doSSign, doInvMassCut, MET, doBJets, "", Syst[i], false, useUnfoldingFiles);
    }
    //----------------------------------------------------------
}

//------------------------------------------------------------
// Close the file if open and delete the pointer
//------------------------------------------------------------
void closeFile(TFile *File)
{
    if (File) {
        if (File->IsOpen()) File->Close();
        std::cout << "Closing: " << File->GetName() << "   --->   Closed ? " << (!(File->IsOpen())) << std::endl;
        //delete File;
        //File = NULL;
    }
}

void closeFiles(TFile *Files[])
{
    if (Files[0]) {
        string fileName = Files[0]->GetName();
        int nFiles;
        //if (fileName.find("Data") != string::npos) {
        if ((fileName.find("Data") != string::npos) && (fileName.find("QCD") == string::npos)){
            nFiles = 3; 
            if (fileName.find("DE") != string::npos) nFiles = 5;
        }
        else if (fileName.find("UNFOLDING") != string::npos) nFiles = 4; 
        else nFiles = 5;

        for (int i(0); i < nFiles; i++){
            closeFile(Files[i]);
        }
    }
}

void closeFiles(TFile *Files[], int nFiles)
{
    string fileName = Files[0]->GetName();
    for (int i(0); i < nFiles; i++){
        closeFile(Files[i]);
        std::cout << "Closing file: " << Files[i]->GetName() << "   --->   Closed ? " << (!(Files[i]->IsOpen())) << std::endl;
    }
}


TH1D* getHisto(TFile *File, const string variable)
{
    TH1D *histo = (TH1D*) File->Get(variable.c_str());
    histo->SetDirectory(0);
    return histo;
}

void getHistos(TH1D *histograms[], TFile *Files[], string variable, bool isDoubleLep)
{
    string fileName = Files[0]->GetName();
    int nFiles;
    //if (fileName.find("Data") != string::npos){
    if ((fileName.find("Data") != string::npos) && (fileName.find("QCD") == string::npos)){
        nFiles = 3; 
        if (fileName.find("DE") != string::npos) nFiles = 5;
    }
    else if (((isDoubleLep && fileName.find("DYJets") != string::npos) || (!isDoubleLep && fileName.find("WJets") != string::npos)) && fileName.find("UNFOLDING") != string::npos) nFiles = 4; 
    else nFiles = 5;

    for (int i(0); i < nFiles; i++){
        histograms[i] = (TH1D*) Files[i]->Get(variable.c_str());
    } 
}







void getStatistics(string leptonFlavor, int JetPtMin, int JetPtMax, bool doFlat, bool doVarWidth, int doQCD, bool doSSign, bool doInvMassCut, int METcut, int doBJets, bool doTTScale)
{

    //variableExc/Inc is the name of the histogram whos bin counts we are recording
    std::string variableExc = "ZNGoodJets_Zexc";
    std::string variableInc = "ZNGoodJets_Zinc";
    
    //string energy = getEnergy();
    //andrew -- just easy fix for now
    string energy = "13TeV";

    std::cout <<"\n>>>>> Getting jet multiplicity statistic!" << std::endl;
    // jet counter
    //int NBins = 11
    int NBins = 8 ;
    double DataEvExc[20][20] = {{0}};
    double DataEvInc[20][20] = {{0}};

    //-- Fetch the data files and histograms --------------
    //andrew -- 2016 w+jets uses NFILESTTBARWJETS
    int usedFiles = NFILESTTBARWJETS ; 
    bool doDY(0) ;
    if (leptonFlavor.find("Muons") != string::npos ||  leptonFlavor.find("Electrons") != string::npos ) {
        usedFiles = NFILESDYJETS;
        doDY = true;
        NBins = 8 ; /// FIXED !!!!
    }
    for (int i(0); i < usedFiles; i++) {
        //"fData" is the name of the file (either data or MC) recording bin counts of
        TFile *fData;
        int sel = i; 
        if (doDY) sel = FilesDYJets[i];
        else if (leptonFlavor.find("SMuE") != string::npos) sel = FilesTTbar[i];
        else sel = FilesTTbarWJets[i];

        //following line skips over any of the QCD control region files
        if ((doQCD > 0 || doInvMassCut || doSSign ) && ProcessInfo[sel].filename.find("QCD") != string::npos) continue;
        
        fData = getFile(FILESDIRECTORY,  leptonFlavor, energy, ProcessInfo[sel].filename, JetPtMin, JetPtMax, doFlat, doVarWidth, doQCD , doSSign,  doInvMassCut, METcut, doBJets,  "","0");
        std::cout << "Opened file #" << sel << ": " << ProcessInfo[sel].filename << std::endl;
        TH1D *hTempExc = getHisto(fData, variableExc);
        TH1D *hTempInc = getHisto(fData, variableInc);
        
        //andrew -- ttbar SFs for 2016 data -- 8 April 2019
        //i==4 is the ttbar file when QCD is turned off (see fileNames.h)
        //i==5 when running the usual distributions with bveto == -1
        //if (i==4){
        if (i==5){
            if (doTTScale){
                std::cout << ">>>>>>> Implementing ttbar SFs!" << std::endl;
                //inclusive jets
                hTempInc->SetBinContent(3, hTempInc->GetBinContent(3)*1.04904374);
                hTempInc->SetBinContent(4, hTempInc->GetBinContent(4)*0.99624372);
                hTempInc->SetBinContent(5, hTempInc->GetBinContent(5)*0.97243212);
                hTempInc->SetBinContent(6, hTempInc->GetBinContent(6)*0.94750793);
                hTempInc->SetBinContent(7, hTempInc->GetBinContent(7)*0.92506162);
                hTempInc->SetBinError(3, hTempInc->GetBinError(3)*1.04904374);
                hTempInc->SetBinError(4, hTempInc->GetBinError(4)*0.99624372);
                hTempInc->SetBinError(5, hTempInc->GetBinError(5)*0.97243212);
                hTempInc->SetBinError(6, hTempInc->GetBinError(6)*0.94750793);
                hTempInc->SetBinError(7, hTempInc->GetBinError(7)*0.92506162);
                
                //exclusive jets
                hTempExc->SetBinContent(3, hTempExc->GetBinContent(3)*1.28532061);
                hTempExc->SetBinContent(4, hTempExc->GetBinContent(4)*1.02880218);
                hTempExc->SetBinContent(5, hTempExc->GetBinContent(5)*0.98994466);
                hTempExc->SetBinContent(6, hTempExc->GetBinContent(6)*0.95915763);
                hTempExc->SetBinContent(7, hTempExc->GetBinContent(7)*0.93695869);
                hTempExc->SetBinError(3, hTempExc->GetBinError(3)*1.28532061);
                hTempExc->SetBinError(4, hTempExc->GetBinError(4)*1.02880218);
                hTempExc->SetBinError(5, hTempExc->GetBinError(5)*0.98994466);
                hTempExc->SetBinError(6, hTempExc->GetBinError(6)*0.95915763);
                hTempExc->SetBinError(7, hTempExc->GetBinError(7)*0.93695869);
            }
            else{
                std::cout << ">>>>>>> No ttbar SFs implemented!" << std::endl;
            }
        }

        for (int j = 1 ; j < NBins + 1 ; j++ ){
            Double_t binContentExc = hTempExc->GetBinContent(j);
            Double_t binContentInc = hTempInc->GetBinContent(j);

            //"DataEv" records the number of bin counts in each of the files
            DataEvExc[i][j] = binContentExc;
            DataEvInc[i][j] = binContentInc;

            //Recording total bin counts of all contributing MC 
            //(The files are indexed from 0 to usedFiles-1 in the loop, so we record sum of all MC in the usedFiles entry)
            //i==0 is the datafile (see fileNames.h)
            if ( i > 0 ) DataEvExc[usedFiles][j]+=int(binContentExc);
            if ( i > 0 ) DataEvInc[usedFiles][j]+=int(binContentInc);
        }
        // close all input root files
        fData->Close();
    }

    std::cout << "\nClosed all files, now creating statistics tables..." << std::endl;
    
    ////////////////////////////////////////////////////////Statistics Tables
    
    ////Create exclusive jets statistics table
    ostringstream nameStrExc;  nameStrExc << "statsTable_" << leptonFlavor << "_" << variableExc << "_JetPtMin_" << JetPtMin;
    if (doInvMassCut) nameStrExc << "_InvMass";
    if (doSSign)   nameStrExc << "_SS";
    if (doQCD > 0) nameStrExc << "_QCD"<<doQCD;
    if (METcut > 0) nameStrExc << "_MET"<< METcut;
    if (doBJets < 0) nameStrExc <<"_BVeto";
    if (doBJets > 0) nameStrExc << "_BJets";
    if (doTTScale) nameStrExc << "_TTSFs";
    if (!doTTScale) nameStrExc << "_noTTSFs";
    nameStrExc << ".tex";
    
    FILE *outFileExc = fopen(nameStrExc.str().c_str(),"w");
    fprintf( outFileExc, "\\footnotesize{\n\\begin{tabular}{l|cccccccc} \n ");
    fprintf( outFileExc, " &  $N_{\\text{jets}} = 0 $ & $N_{\\text{jets}} = 1 $ & $N_{\\text{jets}} = 2 $ & $N_{\\text{jets}} = 3 $ & $N_{\\text{jets}} = 4 $ & $N_{\\text{jets}} = 5 $ & $N_{\\text{jets}} = 6 $ & $N_{\\text{jets}} = 7$ \\\\ \\hline \n ");

    //// Printing statistics of all the MC samples
    for (int i=1; i< usedFiles + 1 ; i++){
        int sel = i ;   
        if ( doDY ) sel = FilesDYJets[i];
        else sel = FilesTTbarWJets[i];
        //Printing information on file
        if (i < usedFiles) fprintf(outFileExc, " %s        & ", ProcessInfo[sel].legend.c_str());
        else {
            fprintf( outFileExc, "\\hline \n");
            fprintf( outFileExc, " MC TOTAL & ");
        }
        //Printing bin contents (rounded to nearest int)
        for (int j = 1 ; j < NBins + 1  ; j++ ) {
            if (j < NBins ) fprintf( outFileExc, "%d & ", int(DataEvExc[i][j]));
            else fprintf( outFileExc, "%d \\\\ \n ", int(DataEvExc[i][j])); ///?????

        }
    }
    
    std::cout << std::endl;

    // print data statistics
    fprintf( outFileExc, "\\hline \n");
    fprintf( outFileExc, " Data          & ");
    for (int j = 1; j< NBins + 1 ; j++){
        if (j < NBins ) fprintf( outFileExc, "%d & ",  int(DataEvExc[0][j]));
        else fprintf( outFileExc, "%d \\\\ \n ",  int(DataEvExc[0][j]));
    }
    // print ratio of MC/Data
    fprintf( outFileExc, " MC/Data Ratio          & ");
    std::cout << "Histogram: " << variableExc << std::endl;
    for (int j=1; j<NBins + 1; j++){
        double temp = DataEvExc[usedFiles][j]/DataEvExc[0][j];
        std:: cout << "Histogram bin number:" << j << "    Total MC: " << DataEvExc[usedFiles][j] << "    Data: " << DataEvExc[0][j] << "    MC/Data Ratio: " << temp << std::endl;
        if (j<NBins) fprintf( outFileExc, "%f & ", float(temp));
        else fprintf( outFileExc, "%f \\\\ \n ",temp);

    }
    fprintf( outFileExc, "\\end{tabular}}");
    fclose(outFileExc);
    
    ////Now create inclusive jets statistics table
    ostringstream nameStrInc;  nameStrInc << "statsTable_" << leptonFlavor << "_" << variableInc << "_JetPtMin_" << JetPtMin;
    if (doInvMassCut) nameStrInc << "_InvMass";
    if (doSSign)   nameStrInc << "_SS";
    if (doQCD > 0) nameStrInc << "_QCD"<<doQCD;
    if (METcut > 0) nameStrInc << "_MET"<< METcut;
    if (doBJets < 0) nameStrInc <<"_BVeto";
    if (doBJets > 0) nameStrInc << "_BJets";
    if (doTTScale) nameStrInc << "_TTSFs";
    if (!doTTScale) nameStrInc << "_noTTSFs";
    nameStrInc << ".tex";
    
    FILE *outFileInc = fopen(nameStrInc.str().c_str(),"w");
    fprintf( outFileInc, "\\footnotesize{\n\\begin{tabular}{l|cccccccc} \n ");
    fprintf( outFileInc, " &  $N_{\\text{jets}} = 0 $ & $N_{\\text{jets}} = 1 $ & $N_{\\text{jets}} = 2 $ & $N_{\\text{jets}} = 3 $ & $N_{\\text{jets}} = 4 $ & $N_{\\text{jets}} = 5 $ & $N_{\\text{jets}} = 6 $ & $N_{\\text{jets}} = 7$ \\\\ \\hline \n ");
    
    //// print statistics of all the MC samples
    for (int i=1; i< usedFiles + 1 ; i++){
        int sel = i ;
        if ( doDY ) sel = FilesDYJets[i];
        else sel = FilesTTbarWJets[i];
        
        if (i < usedFiles) fprintf(outFileInc, " %s        & ", ProcessInfo[sel].legend.c_str());
        else {
            fprintf( outFileInc, "\\hline \n");
            fprintf( outFileInc, " MC TOTAL & ");
        }
        for (int j = 1 ; j < NBins + 1  ; j++ ) {
            if (j < NBins ) fprintf( outFileInc, "%d & ", int(DataEvInc[i][j]));
            else fprintf( outFileInc, "%d \\\\ \n ", int(DataEvInc[i][j]));
            
        }
    }
    
    std::cout << std::endl;
    
    // print data statistics
    fprintf( outFileInc, "\\hline \n");
    fprintf( outFileInc, " Data          & ");
    for (int j = 1; j< NBins + 1 ; j++){
        if (j < NBins ) fprintf( outFileInc, "%d & ",  int(DataEvInc[0][j]));
        else fprintf( outFileInc, "%d \\\\ \n ",  int(DataEvInc[0][j]));
    }
    // print ratio of MC/Data
    fprintf( outFileInc, " MC/Data Ratio          & ");
    std::cout << "Histogram: " << variableInc << std::endl;
    for (int j=1; j<NBins + 1; j++){
        double temp = DataEvInc[usedFiles][j]/DataEvInc[0][j];
        std:: cout << "Histogram bin number:" << j << "    Total MC: " << DataEvInc[usedFiles][j] << "    Data: " << DataEvInc[0][j] << "    MC/Data Ratio: " << temp << std::endl;
        if (j<NBins) fprintf( outFileInc, "%f & ", float(temp));
        else fprintf( outFileInc, "%f \\\\ \n ",temp);
        
    }
    fprintf( outFileInc, "\\end{tabular}}");
    fclose(outFileInc);
}
