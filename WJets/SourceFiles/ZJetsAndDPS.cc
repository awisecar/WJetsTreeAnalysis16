#define PI 3.14159265359
#define DEBUG 0
#define PRINTEVENTINFO 0

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TDatime.h>
#include <TMath.h>
#include <TRandom3.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "standalone_LumiReWeighting.h"
#include "getFilesAndHistograms.h"
#include "functions.h"
#include "HistoSet.h"
#include "ZJetsAndDPS.h"

using namespace std;

ClassImp(ZJetsAndDPS);

void ZJetsAndDPS::Loop(bool hasRecoInfo, bool hasGenInfo, int year, int doQCD, bool doSSign, bool doInvMassCut, int doBJets, int doPUStudy, bool doFlat, bool useRoch, bool doVarWidth){
    std::cout << "\n >>>>>>>>>> ZJetsAndDPS::Loop() >>>>>>>>>> " << std::endl;

    //--- Check weither it is 8 TeV or 13 TeV ---
    //string energy = getEnergy();
    string energy = "13TeV";
    //--------------------------------------
    if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    //--- Counters to check the yields ---
    unsigned int nEvents(0), nEventsIncl0Jets(0);
    unsigned int countEventpassTrig(0), countEvtpassHasMET(0), nEventsPassMETFilter(0), countEventpassLepReq(0), countEventpassBveto(0);
    unsigned int nEventsWithTwoGoodLeptons(0);
    unsigned int nEventsExcl0Jets(0), nEventsExcl1Jets(0), nEventsExcl2Jets(0), nEventsExcl3Jets(0),nEventsIncBJets(0);
    unsigned int GENnEventsIncl0Jets(0), GENnEventsIncl1Jets(0), GENnEventsIncl2Jets(0), GENnEventsIncl3Jets(0);
    double TotalGenWeight(0.), TotalGenWeightPassGEN(0.), TotalGenWeightPassGENPU(0.), TotalGenWeightPassRECO(0.), TotalRecoWeightPassRECO(0.);
    //------------------------------------
    if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    //==========================================================================================================//
    double MTCut(50.);
    double ZMCutLow(71), ZMCutHigh(111);
    //------------------------------------
    bool doZ(true), doW(false), doTT(false), doDR(false);
    if (leptonFlavor == "SingleElectron" || leptonFlavor == "SingleMuon"){
        doW = true; 
        doDR = true;
    }
    if (leptonFlavor == "TTMuE") doTT = true; 
    if (doW || doTT) doZ = false;
    if (fileName.find("_dR_") != string::npos) doDR = true;
    
    // additional muons variables
    double leptonMass(0.00051);  // 
    int LeptonID(11);
    if (leptonFlavor == "Muons" || leptonFlavor == "SingleMuon"){
        leptonMass = 0.105658;
        LeptonID = 13;
    }
    
    // Turn on MET filtering
    bool doMETFiltering = true;

    //==========================================================================================================//
    //         Create output file        //
    //===================================//
    string command = "mkdir -p " + outputDirectory;
    system(command.c_str());
    string outputFileName = CreateOutputFileName(useRoch, doFlat, doPUStudy, doVarWidth, doBJets, doQCD, doSSign , doInvMassCut);
    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    //==========================================================================================================//

    if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    std::cout << std::endl;
    //==========================================================================================================//
    //       Load efficiency tables        //
    //====================================//
    table TableJESunc, LeptIso, LeptID, LeptTrig;

    if (energy == "13TeV"){
        if (year == 2016){
            table TableJESUncertainties("EfficiencyTables/2016_JECUncertainty_Summer16_23Sep2016V4_AK4PF.txt");
            TableJESunc = TableJESUncertainties;
            if (leptonFlavor == "SingleMuon")  {
                table SF_Muon_TightID_ReReco("EfficiencyTables/2016_SMu_SFs_TightId_13TeV_EtaPt.txt");
                table SF_Muon_TightISO_ReReco("EfficiencyTables/2016_SMu_SFs_TightISO_13TeV_EtaPt.txt");
                table SF_Muon_HLTIsoMu24IsoTkMu24_ReReco("EfficiencyTables/2016_SMu_SFs_HLTIsoMu24IsoTkMu24_13TeV_EtaPt.txt");
                LeptID = SF_Muon_TightID_ReReco;
                LeptIso = SF_Muon_TightISO_ReReco;
                LeptTrig = SF_Muon_HLTIsoMu24IsoTkMu24_ReReco;
            }
        }
        if (year == 2017){
            // using 2016 JES uncertainties for now
            table TableJESUncertainties("EfficiencyTables/2016_JECUncertainty_Summer16_23Sep2016V4_AK4PF.txt");
            TableJESunc = TableJESUncertainties;
            if (leptonFlavor == "SingleMuon")  {
                table SF_Muon_TightID_ReReco("EfficiencyTables/2017_SMu_SFs_TightId_13TeV_EtaPt.txt");
                table SF_Muon_TightISO_ReReco("EfficiencyTables/2017_SMu_SFs_TightISO_13TeV_EtaPt.txt");
                table SF_Muon_HLTIsoMu24IsoTkMu24_ReReco("EfficiencyTables/2017_SMu_SFs_HLTIsoMu27_13TeV_EtaPt.txt");
                LeptID = SF_Muon_TightID_ReReco;
                LeptIso = SF_Muon_TightISO_ReReco;
                LeptTrig = SF_Muon_HLTIsoMu24IsoTkMu24_ReReco;
            }
        }
    }
    std::cout << std::endl;
    //==========================================================================================================//
    
    if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    //==========================================================================================================//
    //     Systematics: jec, pu, xsec     //
    //====================================//
    int puYear(0); 
    if (energy == "13TeV") {
        if (year == 2016) puYear = 2016;
        if (year == 2017) puYear = 2017;
    }
    int mode = (systematics == 1) ? direction : 0; // PU
    cout << "Pile Up Distribution: " << puYear << ", Mode: " << mode << endl;
    standalone_LumiReWeighting puWeight(puYear, mode);
   
    int scale(0); // JES
    if (systematics == 2 && direction ==  1) scale =  1;
    if (systematics == 2 && direction == -1) scale = -1;

    double xsec(1.); // XSEC
    if (systematics == 3 && direction ==  1) xsec = 1. + xsecfactor;
    if (systematics == 3 && direction == -1) xsec = 1. - xsecfactor;

    int smearJet(0); // JER
    if (systematics == 4 && direction ==  1) smearJet =  1;
    if (systematics == 4 && direction == -1) smearJet = -1;
    
    int sysLepSF(0); // LepSF
    if (systematics == 5 && direction ==  1) sysLepSF =  1;
    if (systematics == 5 && direction == -1) sysLepSF = -1;
    
    int sysBtagSF(0); // BtagSF
    if (systematics == 6 && direction ==  1) sysBtagSF =  1;
    if (systematics == 6 && direction == -1) sysBtagSF = -1;
    
    double muScale(1.0);
    if (systematics == 7 && direction ==  1) muScale = 1.002;
    if (systematics == 7 && direction == -1) muScale = 0.998;
    
    bool doMer(false); // the number used for MER : 0.006
    double merUncer(0);
    if (systematics == 8 && direction ==  1) doMer = true;
    
    // Wb study
    bool doWbsyst(false);
    double WbSystSF(1.3);
    if (systematics == 9 && direction ==  1) doWbsyst = true;
    
    // bool doRespSyst(false);
    // if (systematics == 10 && direction ==  1) doRespSyst = true;
    
    TRandom3* RandGen = new TRandom3();
    RandGen->SetSeed(12345678);
    if (sysBtagSF != 0) RandGen->SetSeed(23456789);
    
    TRandom3* Rand_MER_Gen = new TRandom3();
    //==========================================================================================================//

    double sumEventW = 0. ;
    cout << "\nMC initial weight :  " << sumEventW <<endl;

    //------------------------------------
    std::cout << "\n-----> Print out variables: " << std::endl;
    std::cout << "---> Switches -- " << std::endl;
    std::cout << "hasRecoInfo: " << hasRecoInfo << std::endl;
    std::cout << "hasGenInfo: " << hasGenInfo << std::endl;
    std::cout << "isData: " << isData << std::endl;
    std::cout << "lumiScale: " << lumiScale << std::endl;
    std::cout << "puScale: " << puScale << std::endl;

    std::cout << "year: " << year << std::endl;

    std::cout << "doW: " << doW << std::endl;
    std::cout << "doZ: " << doZ << std::endl;
    std::cout << "doTT: " << doTT << std::endl;
    std::cout << "LeptonID: " << LeptonID << std::endl;
    std::cout << "leptonFlavor: " << leptonFlavor << std::endl;
    std::cout << "energy: " << energy << std::endl;
    std::cout << "doQCD: " << doQCD << std::endl;
    std::cout << "doBJets: " << doBJets << std::endl;
    std::cout << "doDR: " << doDR << std::endl;
    std::cout << "useEfficiencyCorrection: " << useEfficiencyCorrection << std::endl;
    std::cout << "useTriggerCorrection: " << useTriggerCorrection << std::endl;
    std::cout << "doMETFiltering: " << doMETFiltering << std::endl;
    std::cout << "doVarWidth: " << doVarWidth << std::endl;

    std::cout << "doSSign: " << doSSign << std::endl;
    std::cout << "doPUStudy: " << doPUStudy << std::endl;
    std::cout << "doFlat: " << doFlat << std::endl;
    std::cout << "useRoch: " << useRoch << std::endl;
    std::cout << "doInvMassCut: " << doInvMassCut << std::endl;

    std::cout << "---> Phase Space -- " << std::endl;
    std::cout << "jetPtCutMin: " << jetPtCutMin << std::endl;
    std::cout << "jetPtCutMax: " << jetPtCutMax << std::endl;
    std::cout << "jetEtaCutMin: " << jetEtaCutMin/10. << std::endl;
    std::cout << "jetEtaCutMax: " << jetEtaCutMax/10. << std::endl;
    std::cout << "ZPtCutMin: " << ZPtCutMin << std::endl;
    std::cout << "ZEtaCutMin: " << ZEtaCutMin << std::endl;
    std::cout << "ZEtaCutMax: " << ZEtaCutMax << std::endl;
    std::cout << "ZMCutLow: " << ZMCutLow << std::endl;
    std::cout << "ZMCutHigh: " << ZMCutHigh << std::endl;
    std::cout << "METcut: " << METcut << std::endl;
    std::cout << "MTCut: " << MTCut << std::endl;

    std::cout << "---> Systematic: " << systematics << ", Direction: " << direction << ", XSecFactor: " << xsecfactor << std::endl;

    //==========================================================================================================//
    // Start looping over all the events //
    //===================================//
    printf("\nProcessing : %s    -->   %s \n", fileName.c_str(), outputFileName.c_str()); 

    //--- Initialize the tree branches ---
    Init(hasRecoInfo, hasGenInfo);
    if (fChain == 0) return;
    Long64_t nbytes(0), nb(0);
    Long64_t nentries = fChain->GetEntries();

    // --- Begin Loop All Entries ---
    std::cout << "-----> Begin loop on all entries! " << std::endl;
    std::cout << "-----> Total number of entries: " << nentries << std::endl;
    // eventOfInterest is the event whose content we want to investigate if PRINTEVENTINFO is on
    int eventOfInterest = 4005;
    for (Long64_t jentry(0); jentry < nentries; jentry++){
    //  for (Long64_t jentry(0); jentry < 3000000; jentry++){
        if (PRINTEVENTINFO && jentry == eventOfInterest) cout << "\n" << __LINE__ << " PRINTEVENTINFO: ==================== EVENT INFO for Event # " << eventOfInterest << " ==================== " << endl;

        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;

        if (jentry % 100000 == 0) std::cout << jentry << " of " << nentries << std::endl;
        nb = fChain->GetEntry(jentry);  
        nbytes += nb;
        nEvents++;

        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        //=======================================================================================================//
        //         Computing weight           //
        //====================================//
        //--- weight variable ---
        double weight(1.);
        double genWeight(1.);
        double weightNoSF(1.);
        //-----------------------
        
        weight = weight * lumiScale * xsec;
        
        if (hasRecoInfo && !isData){
            // 0th element of EvtWeights should return
            // main event weight from GenEventInfoProduct
            weight *= EvtWeights->at(0);
            sumEventW += EvtWeights->at(0);
        }

        //--- There is no pile-up so no need to reweight for that ---
        genWeight = weight;
        double genWeightBackup(genWeight);
        TotalGenWeight += genWeightBackup;
        //---
        
        if (hasRecoInfo){

            // old way to do trigger requirements
            // if (energy == "13TeV" && doW) {
            //     if ( (year == 2016) && ((TrigHltMu & 1LL<<11) || (TrigHltMu & 1LL<<16)) ) countEventpassTrig++;
            //     else if ( (year == 2017) && (TrigHltMu & 1LL<<12) ) countEventpassTrig++;
            // }
            if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
            // new way
             if (energy == "13TeV" && doW) {
                 // for WJets 13 TeV 2016, HLT trigger path should be HLTIsoMu24 || HLTIsoTkMu24
                 if ( (year == 2016) && ((MuHltTrgPath1->at(0) == 1) || (MuHltTrgPath2->at(0) == 1)) ) countEventpassTrig++;
                 // for WJets 13 TeV 2017, HLT path of interest is HLT_IsoMu27
                 if ( (year == 2017) && (MuHltTrgPath2->at(0) == 1) ) countEventpassTrig++;
             }
            if (DEBUG) cout << "Stop after line " << __LINE__ << endl;

            if (doW && (METPt->size() > 0)) countEvtpassHasMET++;
        }

        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        //==========================================================================================================//
        //--- MET FILTERING ---
        //--- https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond_2017
        //===================================//
        
        bool passMETFILTER(true);
        if (hasRecoInfo && doMETFiltering){

            // MET filters for data
            if (isData){
                if (year == 2017){
                    passMETFILTER = (
                        (METFilterPath1->at(0) == 1)     //Flag_HBHENoiseFilter
                        && (METFilterPath2->at(0) == 1)  //Flag_HBHENoiseIsoFilter
                        && (METFilterPath3->at(0) == 1)  //Flag_globalSuperTightHalo2016Filter
                        && (METFilterPath4->at(0) == 1)  //Flag_EcalDeadCellTriggerPrimitiveFilter
                        && (METFilterPath5->at(0) == 1)  //Flag_goodVertices
                        && (METFilterPath6->at(0) == 1)  //Flag_eeBadScFilter
                        && (METFilterPath7->at(0) == 1)  //Flag_BadPFMuonFilter
                    );
                    if (passMETFILTER) nEventsPassMETFilter++ ;
                }
            }

            // MET filters for reco MC
            if (!isData){
                if (year == 2017){
                    passMETFILTER = (
                        (METFilterPath1->at(0) == 1)     //Flag_HBHENoiseFilter
                        && (METFilterPath2->at(0) == 1)  //Flag_HBHENoiseIsoFilter
                        && (METFilterPath3->at(0) == 1)  //Flag_globalSuperTightHalo2016Filter
                        && (METFilterPath4->at(0) == 1)  //Flag_EcalDeadCellTriggerPrimitiveFilter
                        && (METFilterPath5->at(0) == 1)  //Flag_goodVertices
                        && (METFilterPath6->at(0) == 1)  //Flag_eeBadScFilter
                        && (METFilterPath7->at(0) == 1)  //Flag_BadPFMuonFilter
                    );
                    if (passMETFILTER) nEventsPassMETFilter++ ;
                }
            }

            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: METFilterPath1->at(0) = " << METFilterPath1->at(0) << endl;
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: METFilterPath2->at(0) = " << METFilterPath2->at(0) << endl;
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: METFilterPath3->at(0) = " << METFilterPath3->at(0) << endl;
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: METFilterPath4->at(0) = " << METFilterPath4->at(0) << endl;
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: METFilterPath5->at(0) = " << METFilterPath5->at(0) << endl;
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: METFilterPath6->at(0) = " << METFilterPath6->at(0) << endl;
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: METFilterPath7->at(0) = " << METFilterPath7->at(0) << endl;
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: passMETFILTER = " <<  passMETFILTER << endl;
        }

		//=======================================================================================================//
        //     Reweighting MC for pileup discrepancy    //
        //==============================================//

        double puWeightFact(1);
        if (hasRecoInfo && !isData){

            // puWeightFact = (double)puWeight.weight(int(EvtPuCntTruth)); // using "true" number of MC PU vertices

            // andrew -- turning off PU reweighting for now!!! - 28 August 2019
            // puWeightFact = (double)puWeight.weight(int(EvtPuCntObs)); // using observed number MC PU vertices

            weight *= puWeightFact;
        }

	    weightNoSF = weight;

        //=======================================================================================================//
        
        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        //=======================================================================================================//
        //         Retrieving leptons          //
        //====================================//
        bool doMuons(leptonFlavor == "Muons" || doW || doTT);
        bool passesLeptonCut(0);
        bool passesLeptonReq(0), passesLeptonAndMT(0), passesBtagReq(1);
        unsigned short nTotLeptons(0), nLeptons(0), nMuons(0), nElectrons(0);
        vector<leptonStruct> leptons, muons, electrons, mets, allLooseLeptons, selLeptons;
        TLorentzVector lep1, lep2, Z;
        leptonStruct lepton1 = {0, 0, 0, 0, 0, 0, 0};
        leptonStruct lepton2 = {0, 0, 0, 0, 0, 0, 0};
        double METphi(0), METpt(0), MT(0);
        int whichMet(0); //  0 - slimmedMETs, 1 - slimmedMETsPUPPI
        
        if (hasRecoInfo) {
            bool eventTrigger = false;

            //--- DO MUONS ---
            // Reco muon cuts here!
            if (doMuons){
                if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
                nTotLeptons = MuEta->size();
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: ==================== RECO MUONS ==================== " << endl;
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: MuEta->size() = " << MuEta->size() << endl;
                
                // Trigger requirement ---
                 if (energy == "13TeV" && doW) {
                     // for WJets 13 TeV 2016, HLT trigger path should be HLTIsoMu24 || HLTIsoTkMu24
                     if ( (year == 2016) && ((MuHltTrgPath1->at(0) == 1) || (MuHltTrgPath2->at(0) == 1)) ) eventTrigger = true;
                     // for WJets 13 TeV 2017, HLT path of interest is HLT_IsoMu27
                     if ( (year == 2017) && (MuHltTrgPath2->at(0) == 1) ) eventTrigger = true;
                 }
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: MuHltTrgPath1->at(0) = " << MuHltTrgPath1->at(0) << endl;
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: MuHltTrgPath2->at(0) = " << MuHltTrgPath2->at(0) << endl;
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: eventTrigger = " << eventTrigger << endl;

                for (unsigned short i(0); i < nTotLeptons; i++) {
                    if (doMer) merUncer = Rand_MER_Gen->Gaus(0, (MuPt->at(i) * 0.006));

                    //grabbing reco muon information as a leptonStruct
                    leptonStruct mu = {(MuPt->at(i) * muScale) + merUncer, MuEta->at(i), MuPhi->at(i), MuE->at(i) * (((MuPt->at(i) * muScale) + merUncer)/MuPt->at(i)), MuCh->at(i), MuPfIso->at(i), 0};
                    if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For mu #" << i << ": pT, eta = " << MuPt->at(i) << ", " << MuEta->at(i) << endl;

                    // pT cut for muon
                    bool muPassesPtCut(false);
                    if ((year == 2016) && (doW && mu.pt >= 26.)) muPassesPtCut = true;
                    if ((year == 2017) && (doW && mu.pt >= 29.)) muPassesPtCut = true;
                    // eta cut for muon
                    bool muPassesEtaLooseCut(fabs(mu.eta) <= 2.4);
                    bool muPassesEtaCut(doW && fabs(mu.eta) <= 2.4);
                    // ID cut for muon
                    bool muPassesIdCut( MuIdTight->at(i) == 1 ); // muon tight ID
                    if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For mu #" << i << " MuIdTight->at(i) = " << MuIdTight->at(i) << endl;
                    // Iso cut for muon
                    bool muPassesIsoCut(doW && MuPfIso->at(i) < 0.15);  
                    if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For mu #" << i << " MuPfIso->at(i) = " << MuPfIso->at(i) << endl;
                    // Iso cut for muon (for doing QCD background)
                    bool muPassesQCDIsoCut(doW && MuPfIso->at(i) >= 0.2); // use 0.15 if you want to cover full Iso space
                    
                    // select the good muons only
                    // all muons that survive eta cut of 2.4 and pT cut of 15
                    if (!doTT && muPassesEtaLooseCut && mu.pt >= 15) muons.push_back(mu);
                    
                    if (PRINTEVENTINFO && jentry == eventOfInterest) {
                        cout << __LINE__ << " PRINTEVENTINFO: Muon analysis cuts --- " << endl;
                        cout << __LINE__ << " PRINTEVENTINFO: For mu #" << i << ", muPassesPtCut = " << muPassesPtCut << endl;
                        cout << __LINE__ << " PRINTEVENTINFO: For mu #" << i << ", muPassesEtaCut = " << muPassesEtaCut << endl;
                        cout << __LINE__ << " PRINTEVENTINFO: For mu #" << i << ", muPassesIdCut = " << muPassesIdCut << endl;
                        cout << __LINE__ << " PRINTEVENTINFO: For mu #" << i << ", eventTrigger = " << eventTrigger << endl;
                        cout << __LINE__ << " PRINTEVENTINFO: For mu #" << i << ", muPassesIsoCut = " << muPassesIsoCut << endl;
                        cout << __LINE__ << " PRINTEVENTINFO: For mu #" << i << ", muPassesQCDIsoCut = " << muPassesQCDIsoCut << endl;
                    }
                    // check if muon passes all required event cuts
                    if (muPassesPtCut && muPassesEtaCut && muPassesIdCut && (!useTriggerCorrection || eventTrigger)){
                        // analysis iso cut
                        if (muPassesIsoCut){  
                            if (doQCD < 2 && leptonFlavor != "SingleElectron") {
                                leptons.push_back(mu); 
                                // if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: Mu #" << i << " passed analysis cuts!" << endl;
                            }
                        } 
                        // QCD isolation cut
                        if (doQCD > 1 && muPassesQCDIsoCut && leptonFlavor != "SingleElectron") leptons.push_back(mu);                        
                    }
                }//End of loop over all the muons
            }
            
            if (DEBUG) cout << "Stop after line " << __LINE__ << endl;

            nMuons = muons.size();
            nElectrons = electrons.size();
            nLeptons = leptons.size();
            
	        // sort leptons by descending pt (these leptons are the muons that pass event selection)
            // this line says, if there's two or more, than we need to pT-sort them
            // if there's only 1 or 0 than there's no need to sort
	        if (nLeptons >= 2) sort(leptons.begin(), leptons.end(), LepDescendingOrder);
			
            vector<leptonStruct> tempVec;
            for ( int iLep = 0 ; iLep < nLeptons ; iLep++){
                tempVec.push_back(leptons[iLep]);
            }
            selLeptons = tempVec ;

            if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        }// end has reco info for leptons

        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;

        //=======================================================================================================//
        //       Retrieving gen leptons        //
        //====================================//
        bool passesGenLeptonCut(0);
        unsigned short nTotGenLeptons(0), nGenLeptons(0), nTotGenPhotons(0);
        vector<leptonStruct> genLeptons;
        vector<int> usedGenPho;
        TLorentzVector genLep1, genLep2, genZ;
        leptonStruct genLepton1, genLepton2;
        double genMT(0.);
        int nuID = 0;

        if (hasGenInfo) {
            
            nTotGenLeptons = GLepBareEta->size();
            nTotGenPhotons = GLepClosePhotEta->size();

            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: GLepBareEta->size() = " << GLepBareEta->size() << endl;
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: GLepClosePhotEta->size() = " << GLepClosePhotEta->size() << endl;

            // getting PDG ID for the neutrino
            if (doW) nuID = 14;
            else if (doW && LeptonID == 11) nuID = 12;
            //-- retriveing generated leptons with status 1
            for (unsigned short i(0); i < nTotGenLeptons; i++) {

                // allow gen lepton to pass only if satisfies the necessary PDG ID
                bool lepSelector( 
                        doW && 
                        (abs(GLepBareId->at(i)) == LeptonID || abs(GLepBareId->at(i)) == nuID)
                    );
                if (!lepSelector) continue ;
                if (!GLepBarePrompt->at(i)) continue ;
                
                // determining gen lepton charge based on PDG ID
                double charge;
                if (abs(GLepBareId->at(i)) == 12 || abs(GLepBareId->at(i)) == 14 || abs(GLepBareId->at(i)) == 16) charge = 0.;
                else if (GLepBareId->at(i) > 0) charge = -1.;
                else charge = 1.;
                
                // grabbing gen lepton info in the form of a leptonStruct
                leptonStruct genLep = {GLepBarePt->at(i), GLepBareEta->at(i), GLepBarePhi->at(i), GLepBareE->at(i), charge, 0., 0.};
                leptonStruct genLepNoFSR = {GLepBarePt->at(i), GLepBareEta->at(i), GLepBarePhi->at(i), GLepBareE->at(i), charge, 0., 0. };
                
                if (lepSelector
                    && GLepBarePrompt->at(i)
                    && GLepBareSt->at(i) == 1
                    && ( abs(GLepBareId->at(i)) == LeptonID || (charge == 0 && doW) ) )
                 {
                    // Only charged lepton(s) will be dressed
                    if( fabs(genLep.charge) > 0 ){
                        TLorentzVector tmpGenLep;
                        tmpGenLep.SetPtEtaPhiM(genLep.pt, genLep.eta, genLep.phi, leptonMass);

                         // loop over all photons
                        for (unsigned short j(0); j < nTotGenPhotons; j++){
                            if( abs(GLepClosePhotSt->at(j)) != 1 || GLepClosePhotPt->at(j) < 0.000001 ) continue;
                            TLorentzVector tmpGenPho;
                            tmpGenPho.SetPtEtaPhiM(GLepClosePhotPt->at(j), GLepClosePhotEta->at(j), GLepClosePhotPhi->at(j), 0.);
                            int used(0);
                            for (unsigned short k(0); k < usedGenPho.size(); k++){
                                if (j == usedGenPho[k]) used = 1;
                            }
                            // dressing leptons using photons within the delR <= 0.1 cone
                            if (deltaR(tmpGenPho.Phi(), tmpGenPho.Eta(), genLepNoFSR.phi, genLepNoFSR.eta) <= 0.1 && !used){
                                tmpGenLep += tmpGenPho;
                                usedGenPho.push_back(j);
                            }
                        }
                        genLep.pt = tmpGenLep.Pt();
                        genLep.eta = tmpGenLep.Eta();
                        genLep.phi = tmpGenLep.Phi();
                        genLep.energy = tmpGenLep.E();
                    }
                    
                    // Event selection on gen lepton
                    if (doZ && genLep.pt >= 20 && fabs(genLep.eta) <= 2.4 && fabs(genLep.charge) > 0){
                        genLeptons.push_back(genLep);
                    }
             
                    // Make gen Muon pT cut match that of reco Muon
                    double genLepPtCut = 25.;
                    if (year == 2016) genLepPtCut = 26.;
                    if (year == 2017) genLepPtCut = 29.;
                    //For WJets, pT and eta cut on the muon, MET cut on the neutrino
                    if (doW && ( (fabs(genLep.charge) > 0 && genLep.pt >= genLepPtCut && fabs(genLep.eta) <= 2.4) || (fabs(genLep.charge) == 0 && genLep.pt >= METcut) ) ){
                        genLeptons.push_back(genLep); 
                        if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: genLep #" << i << " passed analysis cuts!" << endl;
                    }
                }
            } //end loop over gen leptons
            nGenLeptons = genLeptons.size();
            
            if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
            //-- determine if the event passes the leptons requirements
            // for WJets, should have a muon and a neutrino, for ZJets, should have two muons or two electrons
            if (nGenLeptons >= 2){
                
                // sort leptons by descending pt
                sort(genLeptons.begin(), genLeptons.end(), LepDescendingOrder);
                genLepton1 = genLeptons[0];
                genLepton2 = genLeptons[1];
                
                //----- For W+jets -------
                if (doW){
                    //set genLepton1 as the muon, and genLepton2 as the neutrino
                    //assuming that the two highest pT gen leptons are the signal muon and neutrino
                    if (abs(genLeptons[0].charge) > 0 && genLeptons[1].charge == 0){
                        genLepton1 = genLeptons[0];
                        genLepton2 = genLeptons[1];
                        genMT = sqrt(2 *genLepton2.pt * genLepton1.pt * (1 - cos(genLepton2.phi - genLepton1.phi)));
                    }
                    else if (abs(genLeptons[1].charge) > 0 && genLeptons[0].charge == 0){
                        genLepton1 = genLeptons[1];
                        genLepton2 = genLeptons[0];
                        genMT = sqrt(2 *genLepton2.pt * genLepton1.pt * (1 - cos(genLepton2.phi - genLepton1.phi)));
                    }
                    //if the two highest pT leptons aren't exactly one muon and one neutrino, make it fail the MT cut
                    else genMT = -99.0;
                    
                    //MT cut and MET cut on the neutrino
                    if (genMT >= MTCut && genLepton2.pt >= METcut) passesGenLeptonCut = 1;

                    if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: genMT = " << genMT << endl;
                    if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: MET from neutrino = " << genLepton2.pt << endl;
                }

                //----- For Z+jets -------
                // Select the first two leptons with opposite charge
                if (doZ && genLepton1.charge*genLepton2.charge > 0 && nGenLeptons > 2) {
                    //try grabbing the third-highest pT lepton as the second lepton
                    genLepton2 = genLeptons[2];
                }

                // build the TLorentzVectors, the Z candidate and the kinematic
                genLep1.SetPtEtaPhiE(genLepton1.pt, genLepton1.eta, genLepton1.phi, genLepton1.energy);
                genLep2.SetPtEtaPhiE(genLepton2.pt, genLepton2.eta, genLepton2.phi, genLepton2.energy);
                genZ = genLep1 + genLep2; //TLorentzVector for the gen-level W (or Z) boson

                // Apply charge, mass window, pT, and eta cut for Z+jets
                if (doZ && genLepton1.charge*genLepton2.charge < 0 
                        && genZ.M() > ZMCutLow && genZ.M() < ZMCutHigh 
                        && genZ.Eta()*100 > ZEtaCutMin && genZ.Eta()*100 < ZEtaCutMax 
                        && genZ.Pt()>= ZPtCutMin) {
                    passesGenLeptonCut = 1;
                }
            }
        } //end gen lepton if statement
        //=======================================================================================================//

        // andrew - For WJets, passesLeptonReq is true if there is only one muon that passes analysis cuts
        // and no second muon that passes at least 15 GeV and abs(eta) cut of 2.4 (nMuons and nLeptons must be 1 and refer to the same lepton)
        // this secondary muon event veto is done to reduce DYJets BG
        // if (doW && ( (leptonFlavor == "SingleMuon" && nMuons == 1 && nLeptons == 1 && nElectrons == 0) || (leptonFlavor == "SingleElectron" && nMuons == 0 && nLeptons == 1 && nElectrons == 1) )) {
        if (doW && ( leptonFlavor == "SingleMuon" && nMuons == 1 && nLeptons == 1 && nElectrons == 0 )) {
            passesLeptonReq = true;
            countEventpassLepReq++;
        }

        if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: passesLeptonReq = " << passesLeptonReq << endl;

        //=======================================================================================================//
        //          Retrieving AK4 jets           //
        //========================================//
        bool passesJetCut(1); 
        unsigned short nGoodJets(0),nGoodJets_20(0), nTotJets(0), nJetsAdd(0), nJetsNoDRCut(0), nJetsDR02Cut(0), nJetsPt100DR04(0);
        double jetsHT(0); 
        double jetsHT_20(0);
        double XMETscale(0.), YMETscale(0.);            // for calculating METscale
        double  XMETpt(0.), YMETpt(0.) ; // for calculating METscale

        vector<jetStruct> jets, jetsAdditional, jetsPuMva, jets_20, jetsDR02, jetsNoDRCut, jetsPt100DR04; // additional jet collection with pt threshold of 20 GeV
        TLorentzVector leadJ, secondJ, jet1Plus2, jet1Minus2;
        
        TLorentzVector newLeadJ, newSecondJ, newThirdJ, newFourthJ;
        double ForwardJetRapidity(0), BackwardJetRapidity(0);
        vector<TLorentzVector> vJetYOrdered;
        
        int countBJets(0);
	    int countDR02CutBJets(0), countDR04CutBJets(0);
        int countWbBjets(0); // Wb study

        if (hasRecoInfo) {
            int countNJetsVSBeta[10] = {0};
            nTotJets = JetAk04Eta->size();
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: ==================== AK4 JETS ==================== " << endl;
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: JetAk04Eta->size() = " << JetAk04Eta->size() << endl;
            
            if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
            
            //--- loop over all the jets ----------
            for (unsigned short i(0); i < nTotJets; i++) {
                double jetPtTemp(0.); // for calculating METscale
                bool passBJets(0);
                
                //btagging criterion
                //nominal btagging criterion
                if ( (year == 2016) && (JetAk04BDiscCisvV2->at(i) >= 0.8484)) passBJets = true; // Moriond17 for 13 TeV Btag POG recommended discriminator with medium wp cut
                if ( (year == 2017) && (JetAk04BDiscCisvV2->at(i) >= 0.8838)) passBJets = true; // btag medium wp cut for CombinedSecondaryVertexv2 tagger
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", JetAk04BDiscCisvV2->at(i) = " << JetAk04BDiscCisvV2->at(i) << endl;

                //************************* B-tag Veto Correction *******************************//
                
                // Get pt, eta of the i-th jet in the jet loop
                float pt = JetAk04Pt->at(i);
                float eta = JetAk04Eta->at(i);
                
                // Get a uniformly distributed float between 0 and 1
                // This value should change every time we look at a new jet in the jet loop
                float this_rand = RandGen->Rndm(); 

                //************************* Begin emulating btagging efficiencies in MC (MC only)********//

                // Correct for data/MC b-tag eff. discrepancy by updating b-tag score on jet-by-jet basis

                // do if (!isData) instead???
                if (isData == false){
                    if (year == 2016){
                    
                        bool passBJets_SFB_sys_up = passBJets;     // Initialize the systematic_up as the central value
                        bool passBJets_SFB_sys_down = passBJets; // Initialize the systematic_down as the central value
                        
                        // Get jet flavor (hadron definition)
                        int jetflavour= JetAk04HadFlav->at(i);
                        
                        // ---------------- For Truth-Level B-jets --------------- //
                        if (abs(jetflavour)==5){
                            float effb = 0.67298;
                            float effb_corr = 1;
                            if (false) {
                            }
                            else{
                                if (pt < 30.)                 effb = 0.67298  * effb_corr;
                                if (pt >= 30.  && pt < 50.)   effb = 0.67298  * effb_corr;
                                if (pt >= 50.  && pt < 70.)   effb = 0.70732  * effb_corr;
                                if (pt >= 70.  && pt < 100.)  effb = 0.715349 * effb_corr;
                                if (pt >= 100. && pt < 140.)  effb = 0.712663 * effb_corr;
                                if (pt >= 140. && pt < 200.)  effb = 0.694916 * effb_corr;
                                if (pt >= 200. && pt < 300.)  effb = 0.667261 * effb_corr;
                                if (pt >= 300. && pt < 600.)  effb = 0.607964 * effb_corr;
                                if (pt >= 600. && pt < 1000.) effb = 0.489927 * effb_corr;
                                if (pt >= 1000.)              effb = 0.489927 * effb_corr;
                            }
                            
                            //--- CSVv2_Moriond17_B_H.csv values (run period independent)
                            float           SFb = 0.561694*((1.+(0.31439*pt))/(1.+(0.17756*pt)));
                            if (pt < 20.)   SFb = 0.561694*((1.+(0.31439*20))/(1.+(0.17756*20)));
                            if (pt > 1000.) SFb = 0.561694*((1.+(0.31439*1000))/(1.+(0.17756*1000)));
                            
                            float SFb_error = 0.0;
                            if (pt < 20.)                 SFb_error = 0.040213499218225479*2.;
                            if (pt >= 20.  && pt < 30.)   SFb_error = 0.040213499218225479;
                            if (pt >= 30.  && pt < 50.)   SFb_error = 0.014046305790543556;
                            if (pt >= 50.  && pt < 70.)   SFb_error = 0.012372690252959728;
                            if (pt >= 70.  && pt < 100.)  SFb_error = 0.012274007312953472;
                            if (pt >= 100. && pt < 140.)  SFb_error = 0.011465910822153091;
                            if (pt >= 140. && pt < 200.)  SFb_error = 0.012079551815986633;
                            if (pt >= 200. && pt < 300.)  SFb_error = 0.014995276927947998;
                            if (pt >= 300. && pt < 600.)  SFb_error = 0.021414462476968765;
                            if (pt >= 600. && pt < 1000.) SFb_error = 0.032377112656831741;
                            if (pt >= 1000.)              SFb_error = 0.032377112656831741*2.;
                            
                            float SFb_up = SFb + SFb_error;
                            float SFb_down = SFb - SFb_error;
                            
                            // F values for rand comparison
                            float f = 0.0;
                            float f_up = 0.0;
                            float f_down = 0.0;
                            
                            if (SFb <1.0) f = (1.0 - SFb);
                            if (SFb_up <1.0) f_up = (1.0 - SFb_up);
                            if (SFb_down <1.0) f_down = (1.0 - SFb_down);
                            
                            if (SFb > 1.0) f = (1.0 - SFb)/(1.0 - 1.0/effb);
                            if (SFb_up > 1.0) f_up = (1.0 - SFb_up)/(1.0 - 1.0/effb);
                            if (SFb_down > 1.0) f_down = (1.0 - SFb_down)/(1.0 - 1.0/effb);
                            
                            passBJets_SFB_sys_up = passBJets;     // Initialize the systematic_up as the central value
                            passBJets_SFB_sys_down = passBJets; // Initialize the systematic_down as the central value
                            
                            // Untag a tagged jet
                            if ((passBJets==true) && (SFb<1.0) && (this_rand < f)) passBJets = false; // for central value
                            if ((passBJets_SFB_sys_up==true)   && (SFb_up<1.0) && (this_rand < f_up))   passBJets_SFB_sys_up = false; // for systematic_up
                            if ((passBJets_SFB_sys_down==true) && (SFb_down<1.0) && (this_rand < f_down)) passBJets_SFB_sys_down = false; // for sytematic_down
                            
                            
                            // Tag an untagged jet
                            if ((passBJets==false) && (SFb>1.0) && (this_rand < f)) passBJets = true; // for central value
                            if ((passBJets_SFB_sys_up==false)   && (SFb_up>1.0) && (this_rand < f_up))   passBJets_SFB_sys_up = true; // for systematic_up
                            if ((passBJets_SFB_sys_down==false) && (SFb_down>1.0) && (this_rand < f_down)) passBJets_SFB_sys_down = true; // for sytematic_down
                            
                        } // end b-jet section
                        
                        // ---------------- For Truth-Level C-jets--------------- //
                        if (abs(jetflavour)==4){
                            float effc = 0.185922;
                            float effc_corr = 1;
                            if (false) {
                            }
                            else{
                                if (pt < 30.)                 effc = 0.185922 * effc_corr;
                                if (pt >= 30.  && pt < 50.)   effc = 0.185922 * effc_corr;
                                if (pt >= 50.  && pt < 70.)   effc = 0.182578 * effc_corr;
                                if (pt >= 70.  && pt < 100.)  effc = 0.18893 * effc_corr;
                                if (pt >= 100. && pt < 140.)  effc = 0.189251 * effc_corr;
                                if (pt >= 140. && pt < 200.)  effc = 0.186575 * effc_corr;
                                if (pt >= 200. && pt < 300.)  effc = 0.194419 * effc_corr;
                                if (pt >= 300. && pt < 600.)  effc = 0.171999 * effc_corr;
                                if (pt >= 600. && pt < 1000.) effc = 0.143033 * effc_corr;
                                if (pt >= 1000.)              effc = 0.143033 * effc_corr;
                            }
                            
                            //--- CSVv2_Moriond17_B_H.csv values (run period independent)
                            float           SFc = 0.561694*((1.+(0.31439*pt))/(1.+(0.17756*pt)));
                            if (pt < 20.)   SFc = 0.561694*((1.+(0.31439*20))/(1.+(0.17756*20)));
                            if (pt > 1000.) SFc = 0.561694*((1.+(0.31439*1000))/(1.+(0.17756*1000)));
                            
                            float SFc_error = 0.0;
                            if (pt < 20.)                 SFc_error = 0.12064050137996674*2.;
                            if (pt >= 20.  && pt < 30.)   SFc_error = 0.12064050137996674;
                            if (pt >= 30.  && pt < 50.)   SFc_error = 0.042138919234275818;
                            if (pt >= 50.  && pt < 70.)   SFc_error = 0.03711806982755661;
                            if (pt >= 70.  && pt < 100.)  SFc_error = 0.036822021007537842;
                            if (pt >= 100. && pt < 140.)  SFc_error = 0.034397732466459274;
                            if (pt >= 140. && pt < 200.)  SFc_error = 0.0362386554479599;
                            if (pt >= 200. && pt < 300.)  SFc_error = 0.044985830783843994;
                            if (pt >= 300. && pt < 600.)  SFc_error = 0.064243391156196594;
                            if (pt >= 600. && pt < 1000.) SFc_error = 0.097131341695785522;
                            if (pt >= 1000.)              SFc_error = 0.097131341695785522*2.;
                            
                            float SFc_up = SFc + SFc_error;
                            float SFc_down = SFc - SFc_error;
                            
                            // F values for rand comparison
                            float f = 0.0;
                            float f_up = 0.0;
                            float f_down = 0.0;
                            
                            if (SFc <1.0) f = (1.0 - SFc);
                            if (SFc_up <1.0) f_up = (1.0 - SFc_up);
                            if (SFc_down <1.0) f_down = (1.0 - SFc_down);
                            
                            if (SFc > 1.0) f = (1.0 - SFc)/(1.0 - 1.0/effc);
                            if (SFc_up > 1.0) f_up = (1.0 - SFc_up)/(1.0 - 1.0/effc);
                            if (SFc_down > 1.0) f_down = (1.0 - SFc_down)/(1.0 - 1.0/effc);
                            
                            passBJets_SFB_sys_up = passBJets;     // Initialize the systematic_up as the central value
                            passBJets_SFB_sys_down = passBJets; // Initialize the systematic_down as the central value
                            
                            // Untag a tagged jet
                            if ((passBJets==true) && (SFc<1.0) && (this_rand < f)) passBJets = false; // for central value
                            if ((passBJets_SFB_sys_up==true)   && (SFc_up<1.0) && (this_rand < f_up))   passBJets_SFB_sys_up = false; // for systematic_up
                            if ((passBJets_SFB_sys_down==true) && (SFc_down<1.0) && (this_rand < f_down)) passBJets_SFB_sys_down = false; // for sytematic_down
                            
                            // Tag an untagged jet
                            if ((passBJets==false) && (SFc>1.0) && (this_rand < f)) passBJets = true; // for central value
                            if ((passBJets_SFB_sys_up==false)   && (SFc_up>1.0) && (this_rand < f_up))   passBJets_SFB_sys_up = true; // for systematic_up
                            if ((passBJets_SFB_sys_down==false) && (SFc_down>1.0) && (this_rand < f_down)) passBJets_SFB_sys_down = true; // for sytematic_down
                            
                        } // end c-jet section
                        
                        // ---------------- For Truth-Level Light-jets --------------- //
                        if (abs(jetflavour)<4){
                            float eff_l = 0.0200484;
                            float eff_l_corr = 1;

                            if (false) {
                            }
                            else{
                                if (pt < 30.)                 eff_l = 0.0200484 * eff_l_corr;
                                if (pt >= 30.  && pt < 50.)   eff_l = 0.0200484 * eff_l_corr;
                                if (pt >= 50.  && pt < 70.)   eff_l = 0.0171176 * eff_l_corr;
                                if (pt >= 70.  && pt < 100.)  eff_l = 0.0158081 * eff_l_corr;
                                if (pt >= 100. && pt < 140.)  eff_l = 0.0155987 * eff_l_corr;
                                if (pt >= 140. && pt < 200.)  eff_l = 0.0211239 * eff_l_corr;
                                if (pt >= 200. && pt < 300.)  eff_l = 0.0223319 * eff_l_corr;
                                if (pt >= 300. && pt < 600.)  eff_l = 0.030051  * eff_l_corr;
                                if (pt >= 600. && pt < 1000.) eff_l = 0.0427809 * eff_l_corr;
                                if (pt >= 1000.)              eff_l = 0.0427809 * eff_l_corr;
                            }
                            
                            //--- CSVv2_Moriond17_B_H.csv values (run period independent)
                            float           SFlight = 1.0589+(0.000382569*pt)+(-2.4252e-07*(pt*pt))+(2.20966e-10*(pt*pt*pt));
                            if (pt < 20.)   SFlight = 1.0589+(0.000382569*20)+(-2.4252e-07*(20*20))+(2.20966e-10*(20*20*20));
                            if (pt > 1000.) SFlight = 1.0589+(0.000382569*1000)+(-2.4252e-07*(1000*1000))+(2.20966e-10*(1000*1000*1000));
                            
                            float           SFlight_up = SFlight * (1+(0.100485+(3.95509e-05*pt)+(-4.90326e-08*(pt*pt))));
                            if (pt < 20.)   SFlight_up = SFlight * (1+(2*(0.100485+(3.95509e-05*20)+(-4.90326e-08*(20*20)))));
                            if (pt > 1000.) SFlight_up = SFlight * (1+(2*(0.100485+(3.95509e-05*1000)+(-4.90326e-08*(1000*1000)))));
                            
                            float           SFlight_down = SFlight * (1-(0.100485+(3.95509e-05*pt)+(-4.90326e-08*(pt*pt))));
                            if (pt < 20.)   SFlight_down = SFlight * (1-(2*(0.100485+(3.95509e-05*20)+(-4.90326e-08*(20*20)))));
                            if (pt > 1000.) SFlight_down = SFlight * (1-(2*(0.100485+(3.95509e-05*1000)+(-4.90326e-08*(1000*1000)))));
                            
                            // F values for rand comparison
                            float f = 0.0;
                            float f_up = 0.0;
                            float f_down = 0.0;
                            
                            if (SFlight <1.0) f = (1.0 - SFlight);
                            if (SFlight_up <1.0) f_up = (1.0 - SFlight_up);
                            if (SFlight_down <1.0) f_down = (1.0 - SFlight_down);
                            
                            if (SFlight > 1.0) f = (1.0 - SFlight)/(1.0 - 1.0/eff_l);
                            if (SFlight_up > 1.0) f_up = (1.0 - SFlight_up)/(1.0 - 1.0/eff_l);
                            if (SFlight_down > 1.0) f_down = (1.0 - SFlight_down)/(1.0 - 1.0/eff_l);
                            
                            passBJets_SFB_sys_up = passBJets;     // Initialize the systematic_up as the central value
                            passBJets_SFB_sys_down = passBJets; // Initialize the systematic_down as the central value
                            
                            // Untag a tagged jet
                            if ((passBJets==true) && (SFlight<1.0) && (this_rand < f)) passBJets = false; // for central value
                            if ((passBJets_SFB_sys_up==true)   && (SFlight_up<1.0) && (this_rand < f_up))   passBJets_SFB_sys_up = false; // for systematic_up
                            if ((passBJets_SFB_sys_down==true) && (SFlight_down<1.0) && (this_rand < f_down)) passBJets_SFB_sys_down = false; // for sytematic_down
                            
                            // Tag an untagged jet
                            if ((passBJets==false) && (SFlight>1.0) && (this_rand < f)) passBJets = true; // for central value
                            if ((passBJets_SFB_sys_up==false)   && (SFlight_up>1.0) && (this_rand < f_up))   passBJets_SFB_sys_up = true; // for systematic_up
                            if ((passBJets_SFB_sys_down==false) && (SFlight_down>1.0) && (this_rand < f_down)) passBJets_SFB_sys_down = true; // for sytematic_down
                        }   // end light jet section
                        
                        // for btagging SF systematics
                        if (sysBtagSF ==  1) passBJets = passBJets_SFB_sys_up;
                        if (sysBtagSF == -1) passBJets = passBJets_SFB_sys_down;
                        
                        // Wb study
                        if (abs(jetflavour)==5) countWbBjets++;

                    } // end b-tag efficiency SFs for 2016 MC

                    if (year == 2017){
                        
                        // Get jet flavor (hadron definition)
                        int jetflavour= JetAk04HadFlav->at(i);

                        // Initialize the systematic variations as the central value
                        bool passBJets_SFB_sys_up = passBJets; 
                        bool passBJets_SFB_sys_down = passBJets;  
                        
                        // ---------------- For Truth-Level B-jets --------------- //
                        if (abs(jetflavour)==5){

                            // Set the b-tag eff's determined in MC for truth b jets
                            float effb = 1.;
                            if (pt < 30.)                 effb = 0.658148;
                            if (pt >= 30.  && pt < 50.)   effb = 0.658148;
                            if (pt >= 50.  && pt < 70.)   effb = 0.692716;
                            if (pt >= 70.  && pt < 100.)  effb = 0.702473;
                            if (pt >= 100. && pt < 140.)  effb = 0.694655;
                            if (pt >= 140. && pt < 200.)  effb = 0.669370;
                            if (pt >= 200. && pt < 300.)  effb = 0.617494;
                            if (pt >= 300. && pt < 600.)  effb = 0.542055;
                            if (pt >= 600. && pt < 1000.) effb = 0.438289;
                            if (pt >= 1000.)              effb = 0.438289;
                            
                            // Get the central SF -----
                            // From: CSVv2_94XSF_WP_V2_B_F.csv (run period independent), CSVv2 medium WP, "comb" values 
                            // Assuming that the SF is pT-dependent (?)
                            float           SFb = 1.09079*((1.+(0.180764*pt))/(1.+(0.216797*pt)));
                            if (pt < 20.)   SFb = 1.09079*((1.+(0.180764*20.))/(1.+(0.216797*20.)));
                            if (pt > 1000.) SFb = 1.09079*((1.+(0.180764*1000.))/(1.+(0.216797*1000.)));
                            
                            // Grab SF errors for use in systematic variations
                            // Errors listed in the .csv file are symmetric
                            float SFb_error = 0.;
                            if (pt < 20.)                 SFb_error = 0.048631865531206131 * 2.;
                            if (pt >= 20.  && pt < 30.)   SFb_error = 0.048631865531206131;
                            if (pt >= 30.  && pt < 50.)   SFb_error = 0.014063102193176746;
                            if (pt >= 50.  && pt < 70.)   SFb_error = 0.013459851965308189;
                            if (pt >= 70.  && pt < 100.)  SFb_error = 0.012704650871455669;
                            if (pt >= 100. && pt < 140.)  SFb_error = 0.014372736215591431;
                            if (pt >= 140. && pt < 200.)  SFb_error = 0.015085947699844837;
                            if (pt >= 200. && pt < 300.)  SFb_error = 0.033626105636358261;
                            if (pt >= 300. && pt < 600.)  SFb_error = 0.045323032885789871;
                            if (pt >= 600. && pt < 1000.) SFb_error = 0.058395344763994217;
                            if (pt >= 1000.)              SFb_error = 0.058395344763994217 * 2.;
                            float SFb_up = SFb + SFb_error;
                            float SFb_down = SFb - SFb_error;
                            
                            // f values (jet fractions) for comparison to rand
                            // method of computation depends on wheter SF < 1, or SF > 1
                            float f = 0.;
                            float f_up = 0.;
                            float f_down = 0.;
                        
                            // If SF < 1, randomly untag a tagged jet ------
                            // Compute the untag fraction
                            if (SFb < 1.0)           f = (1.0 - SFb);
                            if (SFb_up < 1.0)     f_up = (1.0 - SFb_up);
                            if (SFb_down < 1.0) f_down = (1.0 - SFb_down);
                            // We untag a fraction f of the tagged jets
                            if ((SFb < 1.0)      && (passBJets==true)              && (this_rand < f))       passBJets = false; // for central value
                            if ((SFb_up < 1.0)   && (passBJets_SFB_sys_up==true)   && (this_rand < f_up))    passBJets_SFB_sys_up = false; // for systematic_up
                            if ((SFb_down < 1.0) && (passBJets_SFB_sys_down==true) && (this_rand < f_down))  passBJets_SFB_sys_down = false; // for sytematic_down
                            
                            // If SF > 1, randomly tag an untagged jet ------
                            // Compute the re-tag fraction
                            if (SFb > 1.0)           f = (1.0 - SFb)/(1.0 - 1.0/effb);
                            if (SFb_up > 1.0)     f_up = (1.0 - SFb_up)/(1.0 - 1.0/effb);
                            if (SFb_down > 1.0) f_down = (1.0 - SFb_down)/(1.0 - 1.0/effb);
                            // We re-tag a fraction f of the untagged jets
                            if ((SFb > 1.0)      && (passBJets==false)              && (this_rand < f))      passBJets = true; // for central value
                            if ((SFb_up > 1.0)   && (passBJets_SFB_sys_up==false)   && (this_rand < f_up))   passBJets_SFB_sys_up = true; // for systematic_up
                            if ((SFb_down > 1.0) && (passBJets_SFB_sys_down==false) && (this_rand < f_down)) passBJets_SFB_sys_down = true; // for sytematic_down
                            
                        } // end b-jet section

                        // ---------------- For Truth-Level C-jets--------------- //
                        if (abs(jetflavour)==4){
                        } // end c-jet section

                        // ---------------- For Truth-Level Light-jets --------------- //
                        if (abs(jetflavour)<4){
                        } // end light jet section

                        // for btagging SF systematics
                        if (sysBtagSF ==  1) passBJets = passBJets_SFB_sys_up;
                        if (sysBtagSF == -1) passBJets = passBJets_SFB_sys_down;
                        
                        // Wb study
                        if (abs(jetflavour)==5) countWbBjets++;
                        
                    } // end b-tag efficiency SFs for 2017 MC

                } // --------- End MC-only b-tag efficiency SF section

                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", passBJets = " << passBJets << endl;

                //************************* End B-tag Veto Correction ***********************************//
               
                // grabbing reco jet information in the form of a jetStruct
                jetStruct jet = {JetAk04Pt->at(i), JetAk04Eta->at(i), JetAk04Phi->at(i), JetAk04E->at(i), i, passBJets, 0, 0};
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ": pT, eta = " << JetAk04Pt->at(i) << ", " << JetAk04Eta->at(i) << endl;

                //-- apply jet energy scale uncertainty (need to change the scale when initiating the object)
                double jetEnergyCorr = 0.; 
                bool jetPassesPtCut(jet.pt >= 20); // for MET uncertainty should the cut be before or aftes adding unc.?????
                jetEnergyCorr = TableJESunc.getEfficiency(jet.pt, jet.eta);

                jetPtTemp = jet.pt; // for calculating METscale
                // scale is 0 unless JES uncertainty turned on
                // if JES uncertainty, vary pT and energy by a pT, eta-dependent uncertainty
                jet.pt *= (1 + scale * jetEnergyCorr);
                jet.energy *= (1 + scale * jetEnergyCorr);
				
				// for MET scale (JES uncertainty)
                // if "scale" is non-zero then we need to re-calculate MET to take the JES uncertainties into account
				if (fabs(scale) > 0.){
                    // summed up for all of the jets looped over in the event
					XMETscale += ( scale * jetEnergyCorr * jetPtTemp * cos(jet.phi) ) ;
					YMETscale += ( scale * jetEnergyCorr * jetPtTemp * sin(jet.phi) ) ;
				}
				
                // Now do reco jet event selection cuts
                TLorentzVector jetr;
                jetr.SetPtEtaPhiE(jet.pt, jet.eta, jet.phi, jet.energy);
                // abs(rap) cut
                bool jetPassesEtaCut( (jetr.Rapidity() >= jetEtaCutMin/10.) && (jetr.Rapidity() <= jetEtaCutMax/10.) );
                // PF jet ID: 2017 jet ID is by default tight ID
                bool jetPassesIdCut(JetAk04Id->at(i) > 0);
                // if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", jetPassesIdCut = " << jetPassesIdCut << endl;
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", JetAk04Id->at(i) = " << JetAk04Id->at(i) << endl;
                // PU ID cut to help mitigate pileup jets
                // Loose ID used here
                bool jetPassesPuIdCut(JetAk04PuIdLoose->at(i) == 1);
                // if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", jetPassesPuIdCut = " << jetPassesPuIdCut << endl;
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", JetAk04PuIdLoose->at(i) = " << JetAk04PuIdLoose->at(i) << endl;
                // PU MVA cut (only jets with pT below 100 GeV have to explicitly pass)
                double tempMVA = JetAk04PuMva->at(i);
                bool jetPassesMVACut(0);
                if (energy == "13TeV") {
                    //  for 22Jan rereco, we use simple loose PU ID : -1  - does not pass, 1 passes
                    //  for 13TeV 76x we use the cut for jets having pT < 100 GeV
                    if(jet.pt >= 100.) jetPassesMVACut = true ;
                    else { 
                        if (tempMVA > -0.3) jetPassesMVACut = true ; 
                    }
                }

                if (DEBUG) cout << "Stop after line " << __LINE__ << endl;

                // jet dR cut (wrt to muons)
                bool jetPassesdRCut(1), jetPassesdR02Cut(1);
                unsigned short nRemovedLep = min(int(nLeptons), doW ? 1:2);
                for (unsigned short j(0); j < nRemovedLep; j++) {
                    // determine if passes dRCut
                    // wrt to the leptons that have passed the event selection cuts (should just be the one muon for W+jets)
                    // if the jet and muon are within a dR cone of < 0.4, then jetPassesdRCut turned to 0
                    if ( (doDR) && deltaR(jet.phi, jet.eta, selLeptons[j].phi, selLeptons[j].eta) < 0.4 ) jetPassesdRCut = 0;
					if ( (doDR) && deltaR(jet.phi, jet.eta, selLeptons[j].phi, selLeptons[j].eta) < 0.2 ) jetPassesdR02Cut = 0;
                }
				if (jetPassesdRCut)   jet.passDR04 = true;
				if (jetPassesdR02Cut) jet.passDR02 = true;
				
                if (DEBUG) cout << "Stop after line " << __LINE__ << endl;

                if (PRINTEVENTINFO && jentry == eventOfInterest) {
                    cout << __LINE__ << " PRINTEVENTINFO: Jet analysis cuts --- " << endl;
                    cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", jetPassesPtCut = " << jetPassesPtCut << endl;
                    cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", jetPassesEtaCut = " << jetPassesEtaCut << endl;
                    cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", jetPassesIdCut = " << jetPassesIdCut << endl;
                    cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", jetPassesPuIdCut = " << jetPassesPuIdCut << endl;
                    cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", jetPassesMVACut = " << jetPassesMVACut << endl;
                    cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", jetPassesdRCut = " << jetPassesdRCut << endl;
                }
                // rapidity, PF ID, PU ID, and Pt Cut (lower than analysis cut at 10 GeV)
                if ( jetPassesEtaCut && jetPassesIdCut && jetPassesPuIdCut && jetPassesPtCut) {
                    // getting information about Bjets if doBJets is not 0
                    // passBJets is the marker for if a jet is btagged or not
					if ( fabs(doBJets) > 0 && passBJets == true) countBJets++ ;// count BJets, used for BVeto
					if ( fabs(doBJets) > 0 && passBJets == true && jetPassesdRCut)   countDR04CutBJets++ ;
					if ( fabs(doBJets) > 0 && passBJets == true && jetPassesdR02Cut) countDR02CutBJets++ ;
				
                    // PU MVA cut
					if (energy == "13TeV" && doPUStudy < 0  && jetPassesMVACut) jetsNoDRCut.push_back(jet);
                    // PU MVA cut and dR 0.4 cut
					if (energy == "13TeV" && doPUStudy < 0  && jetPassesMVACut && jetPassesdRCut) jets.push_back(jet);
                }
                // for jetsPuMva
                // rapidity, PF ID, PU ID, dR 0.4 cut, analysis pT cut (30 GeV)
                // w/ no PU MVA ID cut yet
                if (jetPassesEtaCut && jetPassesIdCut && jetPassesPuIdCut && jetPassesdRCut && jet.pt >= jetPtCutMin) {
					jetsPuMva.push_back(jet);
				}
                if (jetPassesEtaCut && jetPassesIdCut && jetPassesPuIdCut && jetPassesdRCut && jet.pt >=  15.) jetsAdditional.push_back(jet);
				
                if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
            } //--- End of loop over all the jets ---
            
            for ( int k = 0 ; k < 10 ; k++){
                ZNGoodJetsBeta_Zexc->Fill(countNJetsVSBeta[k], k  , weight);
            }

            // nGoodJets/jets is what we use for standard analysis histos
            nGoodJets = jets.size();
            nJetsAdd = jetsAdditional.size();
            nJetsNoDRCut = jetsNoDRCut.size();
            
        }  // END IF HAS RECO
        //=======================================================================================================//

        //=======================================================================================================//
        //        Retrieving AK4 gen jets         //
        //========================================//
        bool passesGenJetCut(1);
        unsigned short nGoodGenJets(0), nGenJetsAdd(0), nGoodGenJets_20(0), nTotGenJets(0), nGenJetsNoDRCut(0), nGenJetsDR02Cut(0), nGenJetsPt100DR04(0);
        double genJetsHT(0);
        double genJetsHT_20(0);
        vector<jetStruct> genJets, genJetsAdditional, genJets_20, genJetsDR02, genJetsNoDRCut, genJetsPt100DR04;
        TLorentzVector genLeadJ, genSecondJ, genJet1Plus2, genJet1Minus2;
        
        TLorentzVector genNewLeadJ, genNewSecondJ, genNewThirdJ, genNewFourthJ;
        double genForwardJetRapidity(0), genBackwardJetRapidity(0);
        vector<TLorentzVector> genvJetYOrdered;

        if (hasGenInfo){
            nTotGenJets = GJetAk04Eta->size();

            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: GJetAk04Eta->size() = " << GJetAk04Eta->size() << endl;

            //-- retrieving generated jets
            for (unsigned short i(0); i < nTotGenJets; i++){
                // getting getJet info in the form of a jetStruct
                jetStruct genJet = {GJetAk04Pt->at(i), GJetAk04Eta->at(i), GJetAk04Phi->at(i), GJetAk04E->at(i), i, 0, 0, 0};
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For genJet #" << i << ": pT, eta = " << GJetAk04Pt->at(i) << ", " << GJetAk04Eta->at(i) << endl;
                // genJet dR information (wrt lepton)
                bool genJetPassesdRCut(1), genJetPassesdR02Cut(1);
                for (unsigned short j(0); j < nGenLeptons; j++){ 
                    //andrew - counting both the muon and the neutrino that are in genLeptons here? or just muon?
                    //genLeptons should contain the muon and the neutrino that pass analysis cuts
                    //regardless, our topology suggests that the jets kick off of the W-system, so they should be well seprated from all gen leptonic activity
                    if (genJet.pt >=  jetPtCutMin){
                        gendeltaRjetMu->Fill(deltaR(genJet.phi, genJet.eta, genLeptons[j].phi, genLeptons[j].eta), genWeight); //couting both the dR from muon and neutrino?
                    }
					if ( (genLeptons[j].charge != 0) && (doDR) ){
                        //ah, it's here that we demand our genLepton to be a muon
						if (deltaR(genJet.phi, genJet.eta, genLeptons[j].phi, genLeptons[j].eta) < 0.4) genJetPassesdRCut = 0;
						if (deltaR(genJet.phi, genJet.eta, genLeptons[j].phi, genLeptons[j].eta) < 0.2) genJetPassesdR02Cut = 0;
                    }
                }
				
				if (genJetPassesdRCut)   genJet.passDR04 = true;
				if (genJetPassesdR02Cut) genJet.passDR02 = true;

                if (genJetPassesdRCut && genJet.pt >= 10 && fabs(genJet.eta) <= 4.7){ // Apichart Z+jets uses 5.0
                    genJets.push_back(genJet);                  
                    if (genJet.pt >=  15.) genJetsAdditional.push_back(genJet);
                }		
				if (genJet.pt >= 10 && fabs(genJet.eta) <= 4.7){
					genJetsNoDRCut.push_back(genJet);
				}
            }
            nGoodGenJets = genJets.size();
            nGenJetsAdd = genJetsAdditional.size();
			nGenJetsNoDRCut = genJetsNoDRCut.size();

        } //end if hasGenInfo for gen AK4 jets

        //=======================================================================================================//

        //=======================================================================================================//
        //          Retrieving AK8 jets           //
        //========================================//
        // bool passesAK8JetCut(1); 
        unsigned short nTotJetsAK8(0), nJetsAK8NoDRCut(0), nGoodJetsAK8(0);
        double jetsAK8HT(0); 

        vector<jetStruct> jetsAK8, jetsAK8NoDRCut; 

        int countBJetsAK8(0);
        int countDR04CutBJetsAK8(0);

        if (hasRecoInfo) {
            nTotJetsAK8 = JetAk08Eta->size();
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: ==================== AK8 JETS ==================== " << endl;
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: JetAk08Eta->size() = " << JetAk08Eta->size() << endl;
            if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
            
            //--- loop over all the jets ----------
            for (unsigned short i(0); i < nTotJetsAK8; i++) {
                
                //begin btagging section -------
                //nominal btagging criterion
                bool passBJetsAK8(0);
                if ( (year == 2017) && (JetAk08BDiscCisvV2->at(i) >= 0.8838)) passBJetsAK8 = true; // btag medium wp cut for CombinedSecondaryVertexv2 tagger
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", JetAk08BDiscCisvV2->at(i) = " << JetAk08BDiscCisvV2->at(i) << endl;
                
                // ~~~~~ B-TAGGING SCALE FACTOR CORRECTIONS ~~~~~
                // if (isData == false){
                // }

                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", passBJetsAK8 = " << passBJetsAK8 << endl;
                //end btagging section -------

                jetStruct jetAK8 = {JetAk08Pt->at(i), JetAk08Eta->at(i), JetAk08Phi->at(i), JetAk08E->at(i), i, passBJetsAK8, 0, 0};
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ": pT, eta = " << JetAk08Pt->at(i) << ", " << JetAk08Eta->at(i) << endl;

                // ------------------------------------------------
                //
                // ~~~~~ Jet Energy Scale Uncertainty section ~~~~~
                //
                // ------------------------------------------------


                // Now do reco AK8 jet event selection cuts -----
                // pt cut
                bool jetAK8PassesPtCut(jetAK8.pt >= 200.);
                // abs(rap) cut
                TLorentzVector jetAK8r;
                jetAK8r.SetPtEtaPhiE(jetAK8.pt, jetAK8.eta, jetAK8.phi, jetAK8.energy);
                bool jetAK8PassesEtaCut( (jetAK8r.Rapidity() >= jetEtaCutMin/10.) && (jetAK8r.Rapidity() <= jetEtaCutMax/10.) );
                // PF jet ID: 2017 jet ID is by default tight ID
                bool jetAK8PassesIdCut(JetAk08Id->at(i) > 0);
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", JetAk08Id->at(i) = " << JetAk08Id->at(i) << endl;
                // jet dR cut (wrt to muons)
                bool jetAK8PassesdRCut(1);
                unsigned short nRemovedLep = min(int(nLeptons), doW ? 1:2);
                for (unsigned short j(0); j < nRemovedLep; j++) {
                    // determine if passes dRCut
                    // wrt to the leptons that have passed the event selection cuts (should just be the one muon for W+jets)
                    // if the jet and muon are within a dR cone of < 0.4, then jetPassesdRCut turned to 0
                    if ( (doDR) && deltaR(jetAK8.phi, jetAK8.eta, selLeptons[j].phi, selLeptons[j].eta) < 0.8 ) jetAK8PassesdRCut = 0;
                }
                // the dR(mu, j) cut for AK8 jets is 0.8
				if (jetAK8PassesdRCut) jetAK8.passDR04 = true;

                if (PRINTEVENTINFO && jentry == eventOfInterest) {
                    cout << __LINE__ << " PRINTEVENTINFO: AK8 jet analysis cuts --- " << endl;
                    cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", jetAK8PassesPtCut = " << jetAK8PassesPtCut << endl;
                    cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", jetAK8PassesEtaCut = " << jetAK8PassesEtaCut << endl;
                    cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", jetAK8PassesIdCut = " << jetAK8PassesIdCut << endl;
                    cout << __LINE__ << " PRINTEVENTINFO: For jet #" << i << ", jetAK8PassesdRCut = " << jetAK8PassesdRCut << endl;
                }

                // Rapidity, PF ID, and Pt Cut
                if ( jetAK8PassesPtCut && jetAK8PassesEtaCut && jetAK8PassesIdCut ) {
                    // getting information about Bjets if doBJets is not 0
                    // passBJetsAK8 is the marker for if a jet is btagged or not
					if ( fabs(doBJets) > 0 && passBJetsAK8 == true) countBJetsAK8++ ;// count BJets, used for BVeto
					if ( fabs(doBJets) > 0 && passBJetsAK8 == true && jetAK8PassesdRCut) countDR04CutBJetsAK8++ ;
				
					jetsAK8NoDRCut.push_back(jetAK8);
                }
            }

            nJetsAK8NoDRCut = jetsAK8NoDRCut.size();

        } // END IF HAS RECO FOR AK8 JETS

        //=======================================================================================================//

        //=======================================================================================================//
        //          Retrieving AK8 gen jets           //
        //============================================//
        // bool passesGenJetCut(1);
        unsigned short nTotGenJetsAK8(0), nGenJetsAK8NoDRCut(0), nGoodGenJetsAK8(0);
        double genJetsAK8HT(0);
        vector<jetStruct> genJetsAK8, genJetsAK8NoDRCut;

        if (hasGenInfo) {
            nTotGenJetsAK8 = GJetAk08Eta->size();
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: GJetAk08Eta->size() = " << GJetAk08Eta->size() << endl;

            //-- retrieving generated AK8 jets
            for (unsigned short i(0); i < nTotGenJetsAK8; i++){

                jetStruct genJetAK8 = {GJetAk08Pt->at(i), GJetAk08Eta->at(i), GJetAk08Phi->at(i), GJetAk08E->at(i), i, 0, 0, 0};
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For genJet #" << i << ": pT, eta = " << GJetAk08Pt->at(i) << ", " << GJetAk08Eta->at(i) << endl;

                // genJet dR information (wrt lepton)
                bool genJetAK8PassesdRCut(1);
                for (unsigned short j(0); j < nGenLeptons; j++){ 
                    //genLeptons should contain the muon and the neutrino that pass analysis cuts
                    //regardless, our topology suggests that the jets kick off of the W-system, so they should be well seprated from all gen leptonic activity
					if ( (genLeptons[j].charge != 0) && (doDR) ){
                        //ah, it's here that we demand our genLepton to be a muon
						if (deltaR(genJetAK8.phi, genJetAK8.eta, genLeptons[j].phi, genLeptons[j].eta) < 0.8) genJetAK8PassesdRCut = 0;
                    }
                }
                // the dR(mu, j) cut for AK8 jets is 0.8
				if (genJetAK8PassesdRCut) genJetAK8.passDR04 = true;
	
				if (genJetAK8.pt >= 170. && fabs(genJetAK8.eta) <= 4.7){
					genJetsAK8NoDRCut.push_back(genJetAK8);
				}
            }
			nGenJetsAK8NoDRCut = genJetsAK8NoDRCut.size();

        } // END IF HAS GEN FOR AK8 JETS

        // if (DEBUG) cout << "Stop after line " << __LINE__ << "   " << hasGenInfo <<"    gen Wgh = " << genWeight << "  pass gen cuts = " << passesGenLeptonCut <<"  nGenJets = " << nGoodGenJets <<  endl;
        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;

        //=======================================================================================================//
        //     Matching gen and reco jets (for jet pT smearing)    //
        //=========================================================//
        double XMetSmear(0), YMetSmear(0);
        vector<int> genJetsIndex(nGenJetsNoDRCut, 0);
        vector<vector<int> > matchingTable(nJetsNoDRCut, genJetsIndex);

        vector<int> genJetsAK8Index(nGenJetsAK8NoDRCut, 0);
        vector<vector<int> > matchingTableAK8(nJetsAK8NoDRCut, genJetsAK8Index);
       
        // Only signal MC (W+jets) is hasRecoInfo and hasGenInfo
        if (hasRecoInfo && hasGenInfo){

            // AK4 jets -------------------------------------------------------------------------------------------------------------------
            for (unsigned short i(0); i < nJetsNoDRCut; i++){
                int index(-1);
                double mindR(0.2); // dR of matching cone is half of that used to cluster the jets (for AK4 it's 0.4)
                double dR(9999);

                //Searching for the closest gen-reco jet match within the dR 0.5 cone drawn around the reco jet
                for (unsigned short j(0); j < nGenJetsNoDRCut; j++){
                    dR = deltaR(genJetsNoDRCut[j].phi, genJetsNoDRCut[j].eta, jetsNoDRCut[i].phi, jetsNoDRCut[i].eta);
                    if (dR <= mindR){
                        mindR = dR;
                        index = j;
                    }
                } 

                if (index > -1 ){
                    // so for gen-reco matched jets, use the "scaling method" approach

                    // if a gen-reco jet match found
                    // matchingTable keeps track of which reco-gen jet pairs were matched
                    matchingTable[i][index] = 1;
                    
                    // Smearing of reco MC jets to match jet energy resolution seen in data
                    // andrew -- no JER smearing factors yet! -- 28 august 2019
                    double oldJetPt = jetsNoDRCut[i].pt;
                    // double newJetPt = SmearJetPt(oldJetPt, genJetsNoDRCut[index].pt, jetsNoDRCut[i].eta, smearJet, year);
                    double newJetPt = oldJetPt;

                    // Smearing done for both reco jet pT and energy
                    jetsNoDRCut[i].pt = newJetPt;
                    jetsNoDRCut[i].energy = jetsNoDRCut[i].energy * (newJetPt / oldJetPt);
                    
                    // for recalculating MET
                    XMetSmear += (newJetPt - oldJetPt) * cos(jetsNoDRCut[i].phi);
                    YMetSmear += (newJetPt - oldJetPt) * sin(jetsNoDRCut[i].phi);
                    
                    puMVA_JetsMatchGenJets->Fill(JetAk04PuMva->at(jetsNoDRCut[i].patIndex), weight);
                    jetsEta_JetsMatchGenJets->Fill(JetAk04Eta->at(jetsNoDRCut[i].patIndex), weight);
                }

                else {
                    // if no match found within dR 0.2 cone, then index == -1
                    puMVA_JetsNoMatchGenJets->Fill(JetAk04PuMva->at(jetsNoDRCut[i].patIndex), weight);
                    jetsEta_JetsNoMatchGenJets->Fill(JetAk04Eta->at(jetsNoDRCut[i].patIndex), weight);
                    
                    // now do stochastic smearing method for jets that are not gen-matched!

                    

                    
                }

            }

            // AK8 jets -------------------------------------------------------------------------------------------------------------------
            for (unsigned short i(0); i < nJetsAK8NoDRCut; i++){
                int index(-1);
                double mindR(0.4); // dR of matching cone is half of that used to cluster the jets (for AK8 it's 0.8)
                double dR(9999);

                //Searching for the closest gen-reco jet match within the dR 0.5 cone drawn around the reco jet
                for (unsigned short j(0); j < nGenJetsAK8NoDRCut; j++){
                    dR = deltaR(genJetsAK8NoDRCut[j].phi, genJetsAK8NoDRCut[j].eta, jetsAK8NoDRCut[i].phi, jetsAK8NoDRCut[i].eta);
                    if (dR <= mindR){
                        mindR = dR;
                        index = j;
                    }
                } 

                if (index > -1 ){
                    // if a gen-reco jet match found
                    // matchingTableAK8 keeps track of which reco-gen jet pairs were matched
                    matchingTableAK8[i][index] = 1;

                    // Smearing of reco MC jets to match jet energy resolution seen in data
                    // andrew -- no JER smearing factors yet! -- 28 august 2019
                    double oldJetPt = jetsAK8NoDRCut[i].pt;
                    // double newJetPt = SmearJetPt(oldJetPt, genJetsAK8NoDRCut[index].pt, jetsAK8NoDRCut[i].eta, smearJet, year);
                    double newJetPt = oldJetPt;

                    // Smearing done for both reco jet pT and energy
                    jetsAK8NoDRCut[i].pt = newJetPt;
                    jetsAK8NoDRCut[i].energy = jetsAK8NoDRCut[i].energy * (newJetPt / oldJetPt);
                    
                    // jetsAK8Eta_AK8JetsMatchGenJets->Fill(JetAk08Eta->at(jetsAK8NoDRCut[i].patIndex), weight);
                }

            }
        } //end if hasRecoInfo and hasGenInfo for reco jet smearing
        
        //=======================================================================================================//
        //          Retrieving MET and final selection on leptonic activity             //
        //==============================================================================//
        if (hasRecoInfo){
            if (doW && !(METPt->size() > 0)) continue; // this has no effect, all event has METPt->size() = 2 (76x)
            
            METpt = METPt->at(whichMet);
            TLorentzVector tmpVecMet;
            tmpVecMet.SetPxPyPzE(METPx->at(whichMet), METPy->at(whichMet), 0., METE->at(whichMet));
            //METphi = tmpVecMet.Phi();
            METphi = METPhi->at(whichMet); //grabbing MET Phi value directly from tree rather than calculating it

            if (passesLeptonReq) {
                // recalculate METpt and METphi
                // basically either in the case for JES uncertainties (data), for JER (done for W+jets signal)
                if ( (fabs(scale) > 0.) || (hasRecoInfo && hasGenInfo) ){
                    // get reco MET from branches
                    XMETpt = METPx->at(whichMet);
                    YMETpt = METPy->at(whichMet);
                    // if do JES systematic variations (Data only), we recalculate MET
                    if (fabs(scale) > 0.){
                        XMETpt -= XMETscale;
                        YMETpt -= YMETscale;
                    }
                    // Recalculate MET because Smearing jets
                    if (hasRecoInfo && hasGenInfo){
                        XMETpt -= XMetSmear;
                        YMETpt -= YMetSmear;
                    }
                    // use this recalculated MET
                    TVector2 METvec;
                    METvec.Set(XMETpt, YMETpt);
                    // assign new values for METpt and METphi
                    METpt  = METvec.Mod();
                    METphi = METvec.Phi_mpi_pi(METvec.Phi());
                }
                // lepton1 is the muon that passed event selection
                // should only be this one muon in "leptons" if passesLeptonReq is True
                lepton1 = leptons[0];
                leptonStruct tempMet = {METpt, 0., METphi, METpt, 0, 0, 0};
                lepton2 = tempMet;
                
                MT = sqrt(2 * METpt * lepton1.pt * (1 - cos(METphi - lepton1.phi)));

                // build the TLorentzVectors, the Z candidate and the kinematic
                lep1.SetPtEtaPhiM(lepton1.pt, lepton1.eta, lepton1.phi, leptonMass);
                lep2.SetPtEtaPhiM(METpt, 0, METphi, 0);
                Z = lep1 + lep2; //reco Z or W boson
                
                // correct for identification and isolation efficiencies if required by useEfficiencyCorrection
                // and then trigger efficiencies if useTriggerCorrection
                // Applying scale factors only on reco MC
                if (useEfficiencyCorrection) {
                    double effWeight = 1.;
                    if (leptonFlavor == "SingleMuon") {
                        effWeight = LeptID.getEfficiency(lepton1.pt, fabs(lepton1.eta), sysLepSF);
                        effWeight *= LeptIso.getEfficiency(lepton1.pt, fabs(lepton1.eta), sysLepSF);
                        if (useTriggerCorrection) effWeight *= LeptTrig.getEfficiency(fabs(lepton1.pt), fabs(lepton1.eta), sysLepSF);
                    }
                    if (isData) weight /= effWeight;
                    // The following line should give the emulation of the ID, Iso, Trig efficiencies seen in data to the MC
                    else weight *= effWeight; 
                }

                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: METpt = " << METpt << endl;
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: MT = " << MT << endl;

                // Apply transverse mass and MET cut
                if (METpt >= METcut && passMETFILTER && (((doQCD % 2) == 0 && MT >= MTCut) || ((doQCD % 2) == 1 && MT < MTCut))) {
                    passesLeptonCut = true;
                    passesLeptonAndMT = true;
                    nEventsWithTwoGoodLeptons++;
                }
            }
            
        } // END IF RECO FOR MET
        
        if (passesGenLeptonCut) {
            TotalGenWeightPassGEN += genWeightBackup;
            TotalGenWeightPassGENPU += weight;
        }
        //=======================================================================================================//
        
        
        // Re-analyze the jets collections and Cut on the Pt
        // Needed to smear the jets to account for JER

	    double mindRjmu(99999), mindRj04Pt100mu(99999);
        double  mindRdijet04Pt100mu(99999);
	    int IndClosestDR02Jet(-1), IndEvtClosestDR04Jet(-1), IndEvtClosestDR04DiJet(-1);
	    double genMindRjmu(99999), genMindRj04Pt100mu(99999); 
        double genMindRdijet04Pt100mu(99999);
	    int genIndClosestDR02Jet(-1), genIndEvtClosestDR04Jet(-1), genIndEvtClosestDR04DiJet(-1);

        if (hasRecoInfo){
            
            // AK4 jets -------------------------------------------------------------------------------------------------------------------
            vector<jetStruct> tmpJets;
            for (unsigned short i(0); i < nJetsNoDRCut; i++){
                // the jetsNoDRCut collection are the jets that have been smeared for JER
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: Final reco AK4 jet selections --- " << endl;
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jetsNoDRCut #" << i << ", jetsNoDRCut[i].pt = " << jetsNoDRCut[i].pt << endl;
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jetsNoDRCut #" << i << ", jetsNoDRCut[i].passDR04 = " << jetsNoDRCut[i].passDR04 << endl;
                if (jetsNoDRCut[i].pt >= jetPtCutMin && jetsNoDRCut[i].passDR04) tmpJets.push_back(jetsNoDRCut[i]);
                if (jetsNoDRCut[i].pt >= 20.         && jetsNoDRCut[i].passDR04) jets_20.push_back(jetsNoDRCut[i]);
				if (jetsNoDRCut[i].pt >= jetPtCutMin && jetsNoDRCut[i].passDR02) jetsDR02.push_back(jetsNoDRCut[i]);
				if (jetsNoDRCut[i].pt >= 100.        && jetsNoDRCut[i].passDR04) jetsPt100DR04.push_back(jetsNoDRCut[i]);
            }
            jets.clear(); 
            jets = tmpJets; 
            tmpJets.clear(); 

            // jets is the main collection to pass event selection
            // jets_20 is the collection used for W+jets alpha-s
            // we lower the jet pT cut before unfolding to account for bin migrations up
            // in the lower part of the spectrum
            nGoodJets = jets.size();
            nGoodJets_20 = jets_20.size();

            nJetsDR02Cut = jetsDR02.size();
			nJetsPt100DR04 = jetsPt100DR04.size();

            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: nGoodJets_20 = " << nGoodJets_20 << endl;
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: nGoodJets = " << nGoodJets << endl;

            // AK8 jets -------------------------------------------------------------------------------------------------------------------
            vector<jetStruct> tmpJetsAK8;
            for (unsigned short i(0); i < nJetsAK8NoDRCut; i++){
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: Final reco AK8 jet selections --- " << endl;
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jetsAK8NoDRCut #" << i << ", jetsAK8NoDRCut[i].pt = " << jetsAK8NoDRCut[i].pt << endl;
                if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: For jetsAK8NoDRCut #" << i << ", jetsAK8NoDRCut[i].passDR04 = " << jetsAK8NoDRCut[i].passDR04 << endl;
                if (jetsAK8NoDRCut[i].pt >= 200. && jetsAK8NoDRCut[i].passDR04) tmpJetsAK8.push_back(jetsAK8NoDRCut[i]);
            }
            jetsAK8.clear(); 
            jetsAK8 = tmpJetsAK8; 
            tmpJetsAK8.clear();             

            nGoodJetsAK8 = jetsAK8.size();

            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: nGoodJetsAK8 = " << nGoodJetsAK8 << endl;
			
			//--- Sort jets by decreasing pT ----------------------------------------------------------------------------------------------
            // AK4
			if (nGoodJets    >= 1) sort(jets.begin(), jets.end(), JetDescendingOrder);
            if (nGoodJets_20 >= 1) sort(jets_20.begin(), jets_20.end(), JetDescendingOrder);
			if (nJetsAdd     >= 1) sort(jetsAdditional.begin(), jetsAdditional.end(), JetDescendingOrder);
			if (nJetsDR02Cut >= 1) sort(jetsDR02.begin(), jetsDR02.end(), JetDescendingOrder);
			if (nJetsPt100DR04 >= 1) sort(jetsPt100DR04.begin(), jetsPt100DR04.end(), JetDescendingOrder);
            // AK8
            if (nGoodJetsAK8 >= 1) sort(jetsAK8.begin(), jetsAK8.end(), JetDescendingOrder);
            // -----------------------------------------------------------------------------------------------------------------------------
            
            if (nGoodJets >= 1){
                leadJ.SetPtEtaPhiE(jets[0].pt, jets[0].eta, jets[0].phi, jets[0].energy);               
                newLeadJ.SetPtEtaPhiE(jets[0].pt, jets[0].eta, jets[0].phi, jets[0].energy);
                vector<double> vJetRapidity;
                for (unsigned short i(0); i < nGoodJets; i++) {
                    TLorentzVector LVJet;
                    LVJet.SetPtEtaPhiE(jets[i].pt, jets[i].eta, jets[i].phi, jets[i].energy);
                    vJetRapidity.push_back(LVJet.Rapidity());
                    vJetYOrdered.push_back(LVJet);
                }
                ForwardJetRapidity = *max_element(vJetRapidity.begin(), vJetRapidity.end());
                BackwardJetRapidity = *min_element(vJetRapidity.begin(), vJetRapidity.end());
                sort(vJetYOrdered.begin(), vJetYOrdered.end(), JetYDescendingOrder);
                // jet pT max cut (if desired)
                if (jetPtCutMax > jetPtCutMin){
                    passesJetCut = jets[0].pt < jetPtCutMax;
                }
            }
            if (nGoodJets >= 2){
                secondJ.SetPtEtaPhiE(jets[1].pt, jets[1].eta, jets[1].phi, jets[1].energy);               
                newSecondJ.SetPtEtaPhiE(jets[1].pt, jets[1].eta, jets[1].phi, jets[1].energy);
		        jet1Plus2 = leadJ + secondJ;
                jet1Minus2 = leadJ - secondJ;
            }
            if (nGoodJets >= 3){
                newThirdJ.SetPtEtaPhiE(jets[2].pt, jets[2].eta, jets[2].phi, jets[2].energy);
            }
            if (nGoodJets >= 4){
                newFourthJ.SetPtEtaPhiE(jets[3].pt, jets[3].eta, jets[3].phi, jets[3].energy);
            }

            // Next three little sections grab the jet H_T's for the event
            // AK4 jets
	        jetsHT = 0;
            for (unsigned short i(0); i < nGoodJets; i++){
                jetsHT += jets[i].pt;  
            }
	        jetsHT_20 = 0;
            for (unsigned short i(0); i < nGoodJets_20; i++){
                jetsHT_20 += jets[i].pt;  
            }
            // AK8 jets
            jetsAK8HT = 0.;
            for (unsigned short i(0); i < nGoodJetsAK8; i++){
                jetsAK8HT += jetsAK8[i].pt;
            }

			// calculate the closest jet w.r.t to the muon for jetsDR02 collection
			if (nJetsDR02Cut >= 1 && passesLeptonReq){
				for (unsigned short i(0); i < nJetsDR02Cut; i++){
					double dRjmu(9999);
					dRjmu = fabs(deltaR(jetsDR02[i].phi, jetsDR02[i].eta, lepton1.phi, lepton1.eta));
					if (dRjmu < mindRjmu){
						mindRjmu = dRjmu;
						IndClosestDR02Jet = i;
					}
				}
			}
			// calculate the closest jet w.r.t to the muon for jetsPt100DR04 collection
			if (nJetsPt100DR04 >= 1 && passesLeptonReq){
				for (unsigned short i(0); i < nJetsPt100DR04; i++){
					double dRjmu(9999);
					dRjmu = fabs(deltaR(jetsPt100DR04[i].phi, jetsPt100DR04[i].eta, lepton1.phi, lepton1.eta));
					if (dRjmu < mindRj04Pt100mu){
						mindRj04Pt100mu = dRjmu;
						IndEvtClosestDR04Jet = i;
					}
				}
			}

			// calculate the closest jet in dijet w.r.t to the muon for jetsPt100DR04 collection
			if (nJetsPt100DR04 >= 2 && passesLeptonReq){
                double dRj1mu(9999);
                double dRj2mu(9999);
                dRj1mu = fabs(deltaR(jetsPt100DR04[0].phi, jetsPt100DR04[0].eta, lepton1.phi, lepton1.eta));
                dRj2mu = fabs(deltaR(jetsPt100DR04[1].phi, jetsPt100DR04[1].eta, lepton1.phi, lepton1.eta));

                if (dRj1mu < dRj2mu){
                    mindRdijet04Pt100mu = dRj1mu;
                    IndEvtClosestDR04DiJet = 0;
                }
                if (dRj2mu < dRj1mu){
                    mindRdijet04Pt100mu = dRj2mu;
                    IndEvtClosestDR04DiJet = 1;
                }
            }
        } //end if hasRecoInfo

        if (hasGenInfo){

            // AK4 jets -------------------------------------------------------------------------------------------------------------------
            vector< jetStruct> tmpJets;
            for (unsigned short i(0); i < nGenJetsNoDRCut; i++){
                TLorentzVector gjetr;
                gjetr.SetPtEtaPhiE(genJetsNoDRCut[i].pt, genJetsNoDRCut[i].eta, genJetsNoDRCut[i].phi, genJetsNoDRCut[i].energy);
                if (genJetsNoDRCut[i].pt >= jetPtCutMin && fabs(gjetr.Rapidity()) <= 0.1*jetEtaCutMax && genJetsNoDRCut[i].passDR04) tmpJets.push_back(genJetsNoDRCut[i]);
                if (genJetsNoDRCut[i].pt >= 20          && fabs(gjetr.Rapidity()) <= 0.1*jetEtaCutMax && genJetsNoDRCut[i].passDR04) genJets_20.push_back(genJetsNoDRCut[i]);
				if (genJetsNoDRCut[i].pt >= jetPtCutMin && fabs(gjetr.Rapidity()) <= 0.1*jetEtaCutMax && genJetsNoDRCut[i].passDR02) genJetsDR02.push_back(genJetsNoDRCut[i]);
				if (genJetsNoDRCut[i].pt >= 100			&& fabs(gjetr.Rapidity()) <= 0.1*jetEtaCutMax && genJetsNoDRCut[i].passDR04) genJetsPt100DR04.push_back(genJetsNoDRCut[i]);
            }
            genJets.clear();
            genJets = tmpJets; 
            tmpJets.clear(); 

            // genJets is the main collection to pass event selection
            // genJets_20 is the collection used for W+jets alpha-s
            // we lower the jet pT cut before unfolding to account for bin migrations up
            // in the lower part of the spectrum
            nGoodGenJets = genJets.size();
            nGoodGenJets_20 = genJets_20.size();

            nGenJetsDR02Cut = genJetsDR02.size();
			nGenJetsPt100DR04 = genJetsPt100DR04.size();

            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: nGoodGenJets_20 = " << nGoodGenJets_20 << endl;
            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: nGoodGenJets = " << nGoodGenJets << endl;

            // AK8 jets -------------------------------------------------------------------------------------------------------------------
            vector< jetStruct> tmpJetsAK8;
            for (unsigned short i(0); i < nGenJetsAK8NoDRCut; i++){
                TLorentzVector gjetAK8r;
                gjetAK8r.SetPtEtaPhiE(genJetsAK8NoDRCut[i].pt, genJetsAK8NoDRCut[i].eta, genJetsAK8NoDRCut[i].phi, genJetsAK8NoDRCut[i].energy);
                if (genJetsAK8NoDRCut[i].pt >= 200. && fabs(gjetAK8r.Rapidity()) <= jetEtaCutMax/10. && genJetsAK8NoDRCut[i].passDR04) tmpJetsAK8.push_back(genJetsAK8NoDRCut[i]);
            }
            genJetsAK8.clear();
            genJetsAK8 = tmpJetsAK8; 
            tmpJetsAK8.clear(); 

            nGoodGenJetsAK8 = genJetsAK8.size();

            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: nGoodGenJetsAK8 = " << nGoodGenJetsAK8 << endl;

			//--- Sort jets by decreasing pT ----------------------------------------------------------------------------------------------
            // AK4
			if (nGoodGenJets >= 1) sort(genJets.begin(), genJets.end(), JetDescendingOrder);
			if (nGoodGenJets_20 >= 1) sort(genJets_20.begin(), genJets_20.end(), JetDescendingOrder);
			if (nGenJetsDR02Cut >= 1) sort(genJetsDR02.begin(), genJetsDR02.end(), JetDescendingOrder);
			if (nGenJetsPt100DR04 >= 1) sort(genJetsPt100DR04.begin(), genJetsPt100DR04.end(), JetDescendingOrder);
            // AK8
            if (nGoodGenJetsAK8 >= 1) sort(genJetsAK8.begin(), genJetsAK8.end(), JetDescendingOrder);
            // -----------------------------------------------------------------------------------------------------------------------------

            if (nGoodGenJets >= 1){
                genLeadJ.SetPtEtaPhiE(genJets[0].pt, genJets[0].eta, genJets[0].phi, genJets[0].energy);
                genNewLeadJ.SetPtEtaPhiE(genJets[0].pt, genJets[0].eta, genJets[0].phi, genJets[0].energy);
                vector<double> vJetRapidity;
                for (unsigned short i(0); i < nGoodGenJets; i++) {
                    TLorentzVector LVJet;
                    LVJet.SetPtEtaPhiE(genJets[i].pt, genJets[i].eta, genJets[i].phi, genJets[i].energy);
                    vJetRapidity.push_back(LVJet.Rapidity());
                    genvJetYOrdered.push_back(LVJet);
                }
                genForwardJetRapidity = *max_element(vJetRapidity.begin(), vJetRapidity.end());
                genBackwardJetRapidity = *min_element(vJetRapidity.begin(), vJetRapidity.end());
                sort(genvJetYOrdered.begin(), genvJetYOrdered.end(), JetYDescendingOrder);
                // jet pT max cut (if desired)
                if (jetPtCutMax > jetPtCutMin){
                    passesGenJetCut = genJets[0].pt < jetPtCutMax;
                }
            }
            if (nGoodGenJets >= 2){
                genSecondJ.SetPtEtaPhiE(genJets[1].pt, genJets[1].eta, genJets[1].phi, genJets[1].energy);
                genNewSecondJ.SetPtEtaPhiE(genJets[1].pt, genJets[1].eta, genJets[1].phi, genJets[1].energy);
                genJet1Plus2 = genLeadJ + genSecondJ;
                genJet1Minus2 = genLeadJ - genSecondJ;
            }
            if (nGoodGenJets >= 3){
                genNewThirdJ.SetPtEtaPhiE(genJets[2].pt, genJets[2].eta, genJets[2].phi, genJets[2].energy);
            }
            if (nGoodGenJets >= 4){
                genNewFourthJ.SetPtEtaPhiE(genJets[3].pt, genJets[3].eta, genJets[3].phi, genJets[3].energy);
            }

            // Next three little sections grab the jet H_T's for the event
            genJetsHT = 0.;
            for (unsigned short i(0); i < nGoodGenJets; i++){
                genJetsHT += genJets[i].pt;  
            }
            genJetsHT_20 = 0.;
            for (unsigned short i(0); i < nGoodGenJets_20; i++){
                genJetsHT_20 += genJets[i].pt;  
            }
            genJetsAK8HT = 0.;
            for (unsigned short i(0); i < nGoodGenJetsAK8; i++){
                genJetsAK8HT += genJetsAK8[i].pt;  
            }

			// calculate the closest gen jet w.r.t to the muon for jetsDR02 collection
			if (nGenJetsDR02Cut >= 1 && passesGenLeptonCut){
				for (unsigned short i(0); i < nGenJetsDR02Cut; i++){
					double dRGenjmu(9999);
					dRGenjmu = fabs(deltaR(genJetsDR02[i].phi, genJetsDR02[i].eta, genLepton1.phi, genLepton1.eta));
					if (dRGenjmu < genMindRjmu){
						genMindRjmu = dRGenjmu;
						genIndClosestDR02Jet = i;
					}
				}
			}
			// calculate the closest gen jet w.r.t to the muon for jetsDR02 collection
			if (nGenJetsPt100DR04 >= 1 && passesGenLeptonCut){
				for (unsigned short i(0); i < nGenJetsPt100DR04; i++){
					double dRGenjmu(9999);
					dRGenjmu = fabs(deltaR(genJetsPt100DR04[i].phi, genJetsPt100DR04[i].eta, genLepton1.phi, genLepton1.eta));
					if (dRGenjmu < genMindRj04Pt100mu){
						genMindRj04Pt100mu = dRGenjmu;
						genIndEvtClosestDR04Jet = i;
					}
				}
			}
			// calculate the closest gen jet in dijet w.r.t to the muon for genJetsPt100DR04 collection
			if (nGenJetsPt100DR04 >= 2 && passesGenLeptonCut){
                double dRGenj1mu(9999);
                double dRGenj2mu(9999);
                dRGenj1mu = fabs(deltaR(genJetsPt100DR04[0].phi, genJetsPt100DR04[0].eta, genLepton1.phi, genLepton1.eta));
                dRGenj2mu = fabs(deltaR(genJetsPt100DR04[1].phi, genJetsPt100DR04[1].eta, genLepton1.phi, genLepton1.eta));

                if (dRGenj1mu < dRGenj2mu){
                    genMindRdijet04Pt100mu = dRGenj1mu;
                    genIndEvtClosestDR04DiJet = 0;
                }
                if (dRGenj2mu < dRGenj1mu){
                    genMindRdijet04Pt100mu = dRGenj2mu;
                    genIndEvtClosestDR04DiJet = 1;
                }
            }
        } //end if hasGenInfo
        //=======================================================================================================//

        // Wb study
        if ( doWbsyst && countWbBjets == 1 && hasGenInfo) {
            weight = weight * WbSystSF;
            genWeight = genWeight * WbSystSF;
        }
        
        // Filling puMVA, btag eff histograms before 
        // final event selection for filling analysis histograms
        if (hasRecoInfo) {

            //--- Fill puMVA ---
            for (unsigned short i(0); i < jetsPuMva.size() ; i++){
                if (passesLeptonCut) puMVA->Fill(JetAk04PuMva->at(jetsPuMva[i].patIndex), weight);
            }
          
            //--- For calculating b-tagging efficiency-------------------------
            float btagWP(1.);
            if (year == 2016) btagWP = 0.8484;
            if (year == 2017) btagWP = 0.8838;

            // advantage of using analysis-cuts jets is that we
            // cut out many jets from pileup
            for (unsigned short i(0); i < nGoodJets; i++){
                int jet_ind = jets[i].patIndex;

                // b-flavor jet
                if(fabs(JetAk04HadFlav->at(jet_ind)) == 5){
                    h_pt_eta_b->Fill(jets[i].pt, jets[i].eta, weightNoSF);
                    h_pt_b->Fill(jets[i].pt, weightNoSF);
                    // track which reco b-jets pass the btag score requirement
                    if(JetAk04BDiscCisvV2->at(jet_ind) >= btagWP){
                        h_pt_eta_b_tagged->Fill(jets[i].pt, jets[i].eta, weightNoSF);
                        h_pt_b_tagged->Fill(jets[i].pt, weightNoSF);
                    }
                }

                // c-flavor jet
                else if(fabs(JetAk04HadFlav->at(jet_ind)) == 4){
                    h_pt_eta_c->Fill(jets[i].pt, jets[i].eta, weightNoSF);
                    h_pt_c->Fill(jets[i].pt, weightNoSF);
                    // track which reco c-jets pass the btag score requirement
                    if(JetAk04BDiscCisvV2->at(jet_ind) >= btagWP){
                        h_pt_eta_c_tagged->Fill(jets[i].pt, jets[i].eta, weightNoSF);
                        h_pt_c_tagged->Fill(jets[i].pt, weightNoSF);
                    }
                }

                // light jet (up, down, strange, gluon)
                else {
                    h_pt_eta_udsg->Fill(jets[i].pt, jets[i].eta, weightNoSF);
                    h_pt_udsg->Fill(jets[i].pt, weightNoSF);
                    // track which reco light jets pass the btag score requirement
                    if(JetAk04BDiscCisvV2->at(jet_ind) >= btagWP){
                        h_pt_eta_udsg_tagged->Fill(jets[i].pt, jets[i].eta, weightNoSF);
                        h_pt_udsg_tagged->Fill(jets[i].pt, weightNoSF);
                    }
                }
            }
            //--- End for calculating b-tagging efficiency-------------------------
        }
        //---
        
        //======= Final Selections: =======
        // B-tagging veto/acceptance for event selection
        if (hasRecoInfo){

            // Note: if doBJets == 0, then we're agnostic to btags in terms of event selection
            // doBJets < 0 is a btag veto
            // if doBJets == -1, then veto event if one or more bjets
            // if doBJets == -2, then veto if two or more
            if (doBJets < 0 && countBJets >= fabs(doBJets)) {
                passesLeptonCut = 0;
                passesBtagReq = false;
                nEventsIncBJets++;
            }

            // doBJets > 0 is a btag requirement
            // if doBJets == 2, then veto event if there's not at least two btags
            // used for moving to the ttbar control region for estimating ttbar SFs
            if (doBJets > 0 && doBJets < 99 && countBJets < fabs(doBJets) ){
                passesLeptonCut = 0;
                passesBtagReq = false;
            }
            if (doBJets == 101 && countBJets != 1){
                passesLeptonCut = 0;
                passesBtagReq = false;
            }

            if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: passesBtagReq = " << passesBtagReq << endl;
        }
        //=================================

        if (PRINTEVENTINFO && jentry == eventOfInterest) cout << __LINE__ << " PRINTEVENTINFO: ==================== HISTOGRAM FILLING ==================== " << endl;

        // if (DEBUG) cout << "Stop after line " << __LINE__ << "   " << hasGenInfo <<"    gen Wgh = " << genWeight << "  pass gen cuts = " << passesGenLeptonCut <<"  nGenJets = " << nGoodGenJets <<  endl;
        //=======================================================================================================//
        //   Filling gen and parton histos    //
        //====================================//
        if (hasGenInfo && PRINTEVENTINFO && jentry == eventOfInterest) {
            cout << __LINE__ << " PRINTEVENTINFO: Requirements for filling gen histos --- " << endl;
            cout << __LINE__ << " PRINTEVENTINFO: hasGenInfo = " << hasGenInfo << endl;
            cout << __LINE__ << " PRINTEVENTINFO: passesGenLeptonCut = " << passesGenLeptonCut << endl;
            cout << __LINE__ << " PRINTEVENTINFO: passesGenJetCut = " << passesGenJetCut << endl;
        }
        if (hasGenInfo && passesGenLeptonCut && passesGenJetCut){
            GENnEventsIncl0Jets++;
            genZNGoodJets_Zexc->Fill(nGoodGenJets, genWeight);
            genZNGoodJetsFull_Zexc->Fill(nGoodGenJets, genWeight);
            genZNGoodJets_Zinc->Fill(0., genWeight);
            genZNGoodJetsFull_Zinc->Fill(0., genWeight);
            genZMass_Zinc0jet->Fill(genZ.M(), genWeight);
            genZPt_Zinc0jet->Fill(genZ.Pt(), genWeight);
            genZRapidity_Zinc0jet->Fill(genZ.Rapidity(), genWeight);
            genZEta_Zinc0jet->Fill(genZ.Eta(), genWeight);
            GLepBarePtZinc0jet->Fill(genLep1.Pt(), genWeight);
            GLepBareEtaZinc0jet->Fill(genLep1.Eta(), genWeight);
            if (doZ || doTT){
                GLepBarePtZinc0jet->Fill(genLep2.Pt(), genWeight);
                GLepBareEtaZinc0jet->Fill(genLep2.Eta(), genWeight);
            }
            
            if (nGoodGenJets_20 >= 1){
                genFirstJetPt_Zinc1jet  ->Fill(genJets_20[0].pt, genWeight);
                genFirstJetPt_1_Zinc1jet->Fill(genJets_20[0].pt, genWeight);
                genFirstJetPt_2_Zinc1jet->Fill(genJets_20[0].pt, genWeight);
                genJetsHT_20_Zinc1jet->Fill(genJetsHT_20, genWeight);

                //andrew - alpha-s
                genLeadingJetPt_Zinc1jet->Fill(genJets_20[0].pt, genWeight);
                // genLeadingJetPt_1_Zinc1jet->Fill(genJets_20[0].pt, genWeight);
                genLeadingJetPt_2_Zinc1jet->Fill(genJets_20[0].pt, genWeight);
                genLepPtPlusLeadingJetPt_Zinc1jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                // genLepPtPlusLeadingJetPt_1_Zinc1jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                genLepPtPlusLeadingJetPt_2_Zinc1jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                genZPt_Zinc1jet->Fill(genZ.Pt(), genWeight);
                // genZPt_1_Zinc1jet->Fill(genZ.Pt(), genWeight);
                genZPt_2_Zinc1jet->Fill(genZ.Pt(), genWeight);
                genZPtPlusLeadingJetPt_Zinc1jet->Fill(genZ.Pt()+genJets_20[0].pt, genWeight);
                // genZPtPlusLeadingJetPt_1_Zinc1jet->Fill(genZ.Pt()+genJets_20[0].pt, genWeight);
                genZPtPlusLeadingJetPt_2_Zinc1jet->Fill(genZ.Pt()+genJets_20[0].pt, genWeight);

                if (nGoodGenJets_20 == 1){
                    genLeadingJetPt_Zexc1jet->Fill(genJets_20[0].pt, genWeight);
                    genLeadingJetPt_2_Zexc1jet->Fill(genJets_20[0].pt, genWeight);
                    genLepPtPlusLeadingJetPt_Zexc1jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                    genLepPtPlusLeadingJetPt_2_Zexc1jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                }
            }
            if (nGoodGenJets_20 >= 2){
                genSecondJetPt_Zinc2jet  ->Fill(genJets_20[1].pt, genWeight);
                genSecondJetPt_1_Zinc2jet->Fill(genJets_20[1].pt, genWeight);
                genSecondJetPt_2_Zinc2jet->Fill(genJets_20[1].pt, genWeight);
                genJetsHT_20_Zinc2jet->Fill(genJetsHT_20, genWeight);

                //andrew - alpha-s
                genLeadingJetPt_Zinc2jet->Fill(genJets_20[0].pt, genWeight);
                // genLeadingJetPt_1_Zinc2jet->Fill(genJets_20[0].pt, genWeight);
                genLeadingJetPt_2_Zinc2jet->Fill(genJets_20[0].pt, genWeight);
                genHTover2_Zinc2jet->Fill((genJets_20[0].pt + genJets_20[1].pt)/2. , genWeight);
                // genHTover2_1_Zinc2jet->Fill((genJets_20[0].pt + genJets_20[1].pt)/2. , genWeight);
                genHTover2_2_Zinc2jet->Fill((genJets_20[0].pt + genJets_20[1].pt)/2. , genWeight);
                genLepPtPlusLeadingJetPt_Zinc2jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                // genLepPtPlusLeadingJetPt_1_Zinc2jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                genLepPtPlusLeadingJetPt_2_Zinc2jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                genLepPtPlusHTover2_Zinc2jet->Fill(genLep1.Pt() + (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                // genLepPtPlusHTover2_1_Zinc2jet->Fill(genLep1.Pt() + (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                genLepPtPlusHTover2_2_Zinc2jet->Fill(genLep1.Pt() + (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                genZPt_Zinc2jet->Fill(genZ.Pt(), genWeight);
                // genZPt_1_Zinc2jet->Fill(genZ.Pt(), genWeight);
                genZPt_2_Zinc2jet->Fill(genZ.Pt(), genWeight);
                genZPtPlusLeadingJetPt_Zinc2jet->Fill(genZ.Pt()+genJets_20[0].pt, genWeight);
                // genZPtPlusLeadingJetPt_1_Zinc2jet->Fill(genZ.Pt()+genJets_20[0].pt, genWeight);
                genZPtPlusLeadingJetPt_2_Zinc2jet->Fill(genZ.Pt()+genJets_20[0].pt, genWeight);
                genZPtPlusHTover2_Zinc2jet->Fill(genZ.Pt()+ (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                // genZPtPlusHTover2_1_Zinc2jet->Fill(genZ.Pt()+ (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                genZPtPlusHTover2_2_Zinc2jet->Fill(genZ.Pt()+ (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);

                if (nGoodGenJets_20 == 2){
                    genLeadingJetPt_Zexc2jet->Fill(genJets_20[0].pt, genWeight);
                    genLeadingJetPt_2_Zexc2jet->Fill(genJets_20[0].pt, genWeight);
                    genLepPtPlusLeadingJetPt_Zexc2jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                    genLepPtPlusLeadingJetPt_2_Zexc2jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                }
            }
            if (nGoodGenJets_20 >= 3){
                genThirdJetPt_Zinc3jet  ->Fill(genJets_20[2].pt, genWeight);
                genThirdJetPt_1_Zinc3jet->Fill(genJets_20[2].pt, genWeight);
                genThirdJetPt_2_Zinc3jet->Fill(genJets_20[2].pt, genWeight);
                genJetsHT_20_Zinc3jet->Fill(genJetsHT_20, genWeight);

                //andrew - alpha-s
                genLeadingJetPt_Zinc3jet->Fill(genJets_20[0].pt, genWeight);
                // genLeadingJetPt_1_Zinc3jet->Fill(genJets_20[0].pt, genWeight);
                genLeadingJetPt_2_Zinc3jet->Fill(genJets_20[0].pt, genWeight);
                genHTover2_Zinc3jet->Fill((genJets_20[0].pt + genJets_20[1].pt)/2. , genWeight);
                // genHTover2_1_Zinc3jet->Fill((genJets_20[0].pt + genJets_20[1].pt)/2. , genWeight);
                genHTover2_2_Zinc3jet->Fill((genJets_20[0].pt + genJets_20[1].pt)/2. , genWeight);
                genLepPtPlusLeadingJetPt_Zinc3jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                // genLepPtPlusLeadingJetPt_1_Zinc3jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                genLepPtPlusLeadingJetPt_2_Zinc3jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                genLepPtPlusHTover2_Zinc3jet->Fill(genLep1.Pt() + (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                // genLepPtPlusHTover2_1_Zinc3jet->Fill(genLep1.Pt() + (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                genLepPtPlusHTover2_2_Zinc3jet->Fill(genLep1.Pt() + (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                genZPt_Zinc3jet->Fill(genZ.Pt(), genWeight);
                // genZPt_1_Zinc3jet->Fill(genZ.Pt(), genWeight);
                genZPt_2_Zinc3jet->Fill(genZ.Pt(), genWeight);
                genZPtPlusLeadingJetPt_Zinc3jet->Fill(genZ.Pt()+genJets_20[0].pt, genWeight);
                // genZPtPlusLeadingJetPt_1_Zinc3jet->Fill(genZ.Pt()+genJets_20[0].pt, genWeight);
                genZPtPlusLeadingJetPt_2_Zinc3jet->Fill(genZ.Pt()+genJets_20[0].pt, genWeight);
                genZPtPlusHTover2_Zinc3jet->Fill(genZ.Pt()+ (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                // genZPtPlusHTover2_1_Zinc3jet->Fill(genZ.Pt()+ (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                genZPtPlusHTover2_2_Zinc3jet->Fill(genZ.Pt()+ (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);

                if (nGoodGenJets_20 == 3){
                    genLeadingJetPt_Zexc3jet->Fill(genJets_20[0].pt, genWeight);
                    genLeadingJetPt_2_Zexc3jet->Fill(genJets_20[0].pt, genWeight);
                    genLepPtPlusLeadingJetPt_Zexc3jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                    genLepPtPlusLeadingJetPt_2_Zexc3jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                }
            }
            if (nGoodGenJets_20 >= 4){
                genFourthJetPt_Zinc4jet  ->Fill(genJets_20[3].pt, genWeight);
                genFourthJetPt_1_Zinc4jet->Fill(genJets_20[3].pt, genWeight);
                genFourthJetPt_2_Zinc4jet->Fill(genJets_20[3].pt, genWeight);

                //andrew - alpha-s
                genLeadingJetPt_Zinc4jet->Fill(genJets_20[0].pt, genWeight);
                // genLeadingJetPt_1_Zinc4jet->Fill(genJets_20[0].pt, genWeight);
                genLeadingJetPt_2_Zinc4jet->Fill(genJets_20[0].pt, genWeight);
                genHTover2_Zinc4jet->Fill((genJets_20[0].pt + genJets_20[1].pt)/2. , genWeight);
                // genHTover2_1_Zinc4jet->Fill((genJets_20[0].pt + genJets_20[1].pt)/2. , genWeight);
                genHTover2_2_Zinc4jet->Fill((genJets_20[0].pt + genJets_20[1].pt)/2. , genWeight);
                genLepPtPlusLeadingJetPt_Zinc4jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                // genLepPtPlusLeadingJetPt_1_Zinc4jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                genLepPtPlusLeadingJetPt_2_Zinc4jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                genLepPtPlusHTover2_Zinc4jet->Fill(genLep1.Pt()+(genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                // genLepPtPlusHTover2_1_Zinc4jet->Fill(genLep1.Pt()+(genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                genLepPtPlusHTover2_2_Zinc4jet->Fill(genLep1.Pt()+(genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                genZPt_Zinc4jet->Fill(genZ.Pt(), genWeight);
                // genZPt_1_Zinc4jet->Fill(genZ.Pt(), genWeight);
                genZPt_2_Zinc4jet->Fill(genZ.Pt(), genWeight);
                genZPtPlusLeadingJetPt_Zinc4jet->Fill(genZ.Pt()+genJets_20[0].pt, genWeight);
                // genZPtPlusLeadingJetPt_1_Zinc4jet->Fill(genZ.Pt()+genJets_20[0].pt, genWeight);
                genZPtPlusLeadingJetPt_2_Zinc4jet->Fill(genZ.Pt()+genJets_20[0].pt, genWeight);
                genZPtPlusHTover2_Zinc4jet->Fill(genZ.Pt()+ (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                // genZPtPlusHTover2_1_Zinc4jet->Fill(genZ.Pt()+ (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);
                genZPtPlusHTover2_2_Zinc4jet->Fill(genZ.Pt()+ (genJets_20[0].pt + genJets_20[1].pt)/2., genWeight);

                // if (nGoodGenJets_20 == 4){
                //     genLeadingJetPt_Zexc4jet->Fill(genJets_20[0].pt, genWeight);
                //     genLeadingJetPt_2_Zexc4jet->Fill(genJets_20[0].pt, genWeight);
                //     genLepPtPlusLeadingJetPt_Zexc4jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                //     genLepPtPlusLeadingJetPt_2_Zexc4jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);
                // }
            }
            if (nGoodGenJets_20 >= 5){
                genFifthJetPt_Zinc5jet  ->Fill(genJets_20[4].pt, genWeight);
                genFifthJetPt_1_Zinc5jet->Fill(genJets_20[4].pt, genWeight);
                genFifthJetPt_2_Zinc5jet->Fill(genJets_20[4].pt, genWeight);
            }
            if (nGoodGenJets_20 >= 6){
                genSixthJetPt_Zinc6jet  ->Fill(genJets_20[5].pt, genWeight);
                genSixthJetPt_1_Zinc6jet->Fill(genJets_20[5].pt, genWeight);
            }
            if (nGenJetsDR02Cut >= 1){
                gendRLepCloseJet_Zinc1jet->Fill(genMindRjmu, genWeight);
                gendRLepCloseJet_2_Zinc1jet->Fill(genMindRjmu, genWeight);
                if(genJetsDR02[0].pt > 500){
                    gendRLepCloseJetCo_Zinc1jet->Fill(genMindRjmu, genWeight);
                    gendRLepCloseJetCo_2_Zinc1jet->Fill(genMindRjmu, genWeight);
                }
            }
            if (nGenJetsPt100DR04 >= 1 && genJetsPt100DR04[0].pt > 300){
                gendRptmin100LepCloseJetCo300dR04_Zinc1jet->Fill(genMindRj04Pt100mu, genWeight);
                gendRptmin100LepCloseJetCo300dR04_2_Zinc1jet->Fill(genMindRj04Pt100mu, genWeight);
            }
            if (nGenJetsPt100DR04 >= 2 && genJetsPt100DR04[0].pt > 300 && deltaPhi(genLeadJ, genSecondJ) > (PI - 0.3)){
                gendRptmin100LepCloseDiJetCo300dR04_Zinc2jet->Fill(genMindRdijet04Pt100mu, genWeight);
                gendRptmin100LepCloseDiJetCo300dR04_2_Zinc2jet->Fill(genMindRdijet04Pt100mu, genWeight);
            }
            if (nGoodGenJets >= 1){
                GENnEventsIncl1Jets++;
                genZNGoodJets_Zinc->Fill(1., genWeight);
                genZNGoodJetsFull_Zinc->Fill(1., genWeight);
                genZMass_Zinc1jet->Fill(genZ.M(), genWeight);
                GLepBarePtZinc1jet->Fill(genLep1.Pt(), genWeight);
                GLepBareEtaZinc1jet->Fill(genLep1.Eta(), genWeight);
                if (doZ || doTT){
                    GLepBarePtZinc1jet->Fill(genLep2.Pt(), genWeight);
                    GLepBareEtaZinc1jet->Fill(genLep2.Eta(), genWeight);
                }
                // genZPt_Zinc1jet->Fill(genZ.Pt(), genWeight);
                genZRapidity_Zinc1jet->Fill(genZ.Rapidity(), genWeight);
                genZEta_Zinc1jet->Fill(genZ.Eta(), genWeight);
                
                genFirstJetEta_Zinc1jet->Fill(fabs(genLeadJ.Eta()), genWeight);
                genFirstJetEta_2_Zinc1jet->Fill(fabs(genLeadJ.Eta()), genWeight);
                genFirstJetAbsRapidity_Zinc1jet->Fill(fabs(genNewLeadJ.Rapidity()), genWeight);
                genFirstJetAbsRapidity_2_Zinc1jet->Fill(fabs(genNewLeadJ.Rapidity()), genWeight);
                genFirstJetRapidityFull_Zinc1jet->Fill(genNewLeadJ.Rapidity(), genWeight);
                genMeanNJetsHT_1D_Zinc1jet->Fill(genJetsHT, genWeight*nGoodGenJets);
                genMeanNJetsHT_Zinc1jet->Fill(genJetsHT, nGoodGenJets, genWeight);
                genFirstJetPtEta_Zinc1jet->Fill(genLeadJ.Pt(), fabs(genLeadJ.Eta()), genWeight);
                genFirstHighestJetPt_Zinc1jet->Fill(genLeadJ.Pt(), genWeight);
                genJetsHT_Zinc1jet->Fill(genJetsHT, genWeight);
                genJetsHT_20_30_Zinc1jet->Fill(genJetsHT_20, genWeight);
                genJetsHT_1_Zinc1jet->Fill(genJetsHT, genWeight);
                genJetsHT_2_Zinc1jet->Fill(genJetsHT, genWeight);
                
                gendPhiLepJet1_Zinc1jet->Fill(deltaPhi(genLep1, genNewLeadJ), genWeight);
                gendPhiLepJet1_2_Zinc1jet->Fill(deltaPhi(genLep1, genNewLeadJ), genWeight);

                genMT_Zinc1jet->Fill(genMT, genWeight);

                for ( int i =0 ; i < NbinsEta2D - 1 ; i++){
                    if ( fabs(genLeadJ.Eta()) >= j_Y_range[i] &&  fabs(genLeadJ.Eta()) < j_Y_range[i+1] ) genFirstJetPt_Zinc1jet_Eta[i]->Fill(fabs(genLeadJ.Pt()), genWeight);
                }
                if ( doW ) gendEtaBosonJet_Zinc1jet->Fill(fabs(genLeadJ.Eta() - genLep1.Eta()), genWeight);
                else gendEtaBosonJet_Zinc1jet->Fill(fabs(genLeadJ.Eta()-genZ.Eta()), genWeight);
                if (nGoodGenJets == 1){
                    genFirstJetPt_Zexc1jet->Fill(genLeadJ.Pt(), genWeight);
                    //    gendEtaBosonJet_Zexc1jet->Fill(fabs(genLeadJ.Eta()), genWeight);
                    if ( doW ) gendEtaBosonJet_Zexc1jet->Fill(fabs(genLeadJ.Eta() - genLep1.Eta()), genWeight);
                    else gendEtaBosonJet_Zexc1jet->Fill(fabs(genLeadJ.Eta()-genZ.Eta()), genWeight);

                }
            }
            if (nGoodGenJetsAK8 >= 1){
                genLeadingJetAK8Pt_Zinc1jet->Fill(genJetsAK8[0].pt, genWeight);
                genLeadingJetAK8Pt_2_Zinc1jet->Fill(genJetsAK8[0].pt, genWeight);
                genLepPtPlusLeadingJetAK8Pt_Zinc1jet->Fill(genLep1.Pt()+genJetsAK8[0].pt, genWeight);
                genLepPtPlusLeadingJetAK8Pt_2_Zinc1jet->Fill(genLep1.Pt()+genJetsAK8[0].pt, genWeight);

                if (nGoodGenJetsAK8 == 1){
                    genLeadingJetAK8Pt_Zexc1jet->Fill(genJetsAK8[0].pt, genWeight);
                    genLeadingJetAK8Pt_2_Zexc1jet->Fill(genJetsAK8[0].pt, genWeight);
                    genLepPtPlusLeadingJetAK8Pt_Zexc1jet->Fill(genLep1.Pt()+genJetsAK8[0].pt, genWeight);
                    genLepPtPlusLeadingJetAK8Pt_2_Zexc1jet->Fill(genLep1.Pt()+genJetsAK8[0].pt, genWeight);

                }
            }
            if (nGoodGenJets >= 2){
                TLorentzVector genJet1Plus2PlusZ = genJet1Plus2 + genZ;
                GENnEventsIncl2Jets++;
                genZNGoodJets_Zinc->Fill(2., genWeight);
                genZNGoodJetsFull_Zinc->Fill(2., genWeight);
                genZMass_Zinc2jet->Fill(genZ.M(), genWeight);
                genTwoJetsPtDiff_Zinc2jet->Fill(genJet1Minus2.Pt(), genWeight);
                // genBestTwoJetsPtDiff_Zinc2jet->Fill(genBestJet1Minus2.Pt(), genWeight);
                genJetsMass_Zinc2jet->Fill(genJet1Plus2.M(), genWeight);
                GLepBarePtZinc2jet->Fill(genLep1.Pt(), genWeight);
                GLepBareEtaZinc2jet->Fill(genLep1.Eta(), genWeight);
                if (doZ || doTT){
                    GLepBarePtZinc2jet->Fill(genLep2.Pt(), genWeight);
                    GLepBareEtaZinc2jet->Fill(genLep2.Eta(), genWeight);
                }
                // genZPt_Zinc2jet->Fill(genZ.Pt(), genWeight);
                genZRapidity_Zinc2jet->Fill(genZ.Rapidity(), genWeight);
                genZEta_Zinc2jet->Fill(genZ.Eta(), genWeight);
                genFirstHighestJetPt_Zinc2jet->Fill(genLeadJ.Pt(), genWeight);
                genSecondHighestJetPt_Zinc2jet->Fill(genSecondJ.Pt(), genWeight);


                genSecondJetEta_Zinc2jet->Fill(fabs(genSecondJ.Eta()), genWeight);
                genSecondJetEta_2_Zinc2jet->Fill(fabs(genSecondJ.Eta()), genWeight);
                genSecondJetAbsRapidity_Zinc2jet->Fill(fabs(genNewSecondJ.Rapidity()), genWeight);
                genSecondJetAbsRapidity_2_Zinc2jet->Fill(fabs(genNewSecondJ.Rapidity()), genWeight);
                genSecondJetRapidityFull_Zinc2jet->Fill(genNewSecondJ.Rapidity(), genWeight);
                gendRapidityJets_Zinc2jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), genWeight);
                gendRapidityJets_2_Zinc2jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), genWeight);
                
                gendRapidityJetsFB_Zinc2jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, genWeight);
                gendRapidityJetsFB_2_Zinc2jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, genWeight);
                
                gendiJetMass_Zinc2jet->Fill(genJet1Plus2.M(), genWeight);
                gendiJetMass_2_Zinc2jet->Fill(genJet1Plus2.M(), genWeight);
                
                gendPhiJetsFB_Zinc2jet->Fill(deltaPhi(genvJetYOrdered[0], genvJetYOrdered[genvJetYOrdered.size() - 1] ), genWeight);
                gendPhiJetsFB_2_Zinc2jet->Fill(deltaPhi(genvJetYOrdered[0], genvJetYOrdered[genvJetYOrdered.size() - 1] ), genWeight);
                
                gendRJets_Zinc2jet->Fill(deltaRYPhi(genNewLeadJ, genNewSecondJ), genWeight);
                gendRJets_2_Zinc2jet->Fill(deltaRYPhi(genNewLeadJ, genNewSecondJ), genWeight);
                
                gendiJetPt_Zinc2jet->Fill(genJet1Plus2.Pt(), genWeight);
                gendiJetPt_2_Zinc2jet->Fill(genJet1Plus2.Pt(), genWeight);
                
                gendPhiLepJet2_Zinc2jet->Fill(deltaPhi(genLep1, genNewSecondJ), genWeight);
                gendPhiLepJet2_2_Zinc2jet->Fill(deltaPhi(genLep1, genNewSecondJ), genWeight);
                
                genMeanNJetsHT_1D_Zinc2jet->Fill(genJetsHT, genWeight*nGoodGenJets);
                genMeanNJetsdRapidity_1D_Zinc2jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), genWeight*nGoodGenJets);
                genMeanNJetsdRapidityFB_1D_Zinc2jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, genWeight*nGoodGenJets);
                
                genMeanNJetsHT_Zinc2jet->Fill(genJetsHT, nGoodGenJets, genWeight);
                genMeanNJetsdRapidity_Zinc2jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), nGoodGenJets, genWeight);
                genMeanNJetsdRapidityFB_Zinc2jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, nGoodGenJets, genWeight);
                genSecondJetPtEta_Zinc2jet->Fill(genSecondJ.Pt(), fabs(genSecondJ.Eta()), genWeight);
                genRatioJetPt21_Zinc2jet->Fill(genJets[1].pt/genJets[0].pt, genWeight);
                genJetsHT_Zinc2jet->Fill(genJetsHT, genWeight);
                genJetsHT_20_30_Zinc2jet->Fill(genJetsHT_20, genWeight);
                genJetsHT_1_Zinc2jet->Fill(genJetsHT, genWeight);
                genJetsHT_2_Zinc2jet->Fill(genJetsHT, genWeight);
                
                genptBal_Zinc2jet->Fill(genJet1Plus2PlusZ.Pt(), genWeight);
                
                gendPhiJets_Zinc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                gendPhiJets_2_Zinc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                
                // genBestdPhiJets_Zinc2jet->Fill(deltaPhi(genBestTwoJets.first, genBestTwoJets.second), genWeight);
                gendEtaJets_Zinc2jet->Fill(genLeadJ.Eta() - genSecondJ.Eta(), genWeight);
                gendEtaFirstJetZ_Zinc2jet->Fill(genLeadJ.Eta() - genZ.Eta(), genWeight);
                gendEtaSecondJetZ_Zinc2jet->Fill(genSecondJ.Eta() - genZ.Eta(), genWeight);
                gendEtaJet1Plus2Z_Zinc2jet->Fill(genJet1Plus2.Eta() - genZ.Eta(), genWeight);
                genPHI_Zinc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                // genBestPHI_Zinc2jet->Fill(PHI(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                genPHI_T_Zinc2jet->Fill(PHI_T(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                // genBestPHI_T_Zinc2jet->Fill(PHI_T(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                genSpT_Zinc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                // genBestSpT_Zinc2jet->Fill(SpT(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                genSpTJets_Zinc2jet->Fill(SpTsub(genLeadJ, genSecondJ), genWeight);
                // genBestSpTJets_Zinc2jet->Fill(SpTsub(genBestTwoJets.first, genBestTwoJets.second), genWeight);
                genSpTLeptons_Zinc2jet->Fill(SpTsub(genLep1, genLep2), genWeight);
                genSPhi_Zinc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                // genBestSPhi_Zinc2jet->Fill(SPhi(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                for ( int i =0 ; i < NbinsEta2D - 1 ; i++){
                    if ( fabs(genSecondJ.Eta()) >= j_Y_range[i] &&  fabs(genSecondJ.Eta()) < j_Y_range[i+1] ) genSecondJetPt_Zinc2jet_Eta[i]->Fill(fabs(genSecondJ.Pt()), genWeight);
                }
                
                if (genZ.Pt() < 25){
                    genptBal_LowPt_Zinc2jet->Fill(genJet1Plus2PlusZ.Pt(), genWeight);
                    gendPhiJets_LowPt_Zinc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                    // genBestdPhiJets_LowPt_Zinc2jet->Fill(deltaPhi(genBestTwoJets.first, genBestTwoJets.second), genWeight);
                    gendPhiLeptons_LowPt_Zinc2jet->Fill(deltaPhi(genLep1, genLep2), genWeight);
                    genPHI_T_LowPt_Zinc2jet->Fill(PHI_T(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    // genBestPHI_T_LowPt_Zinc2jet->Fill(PHI_T(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                    genPHI_LowPt_Zinc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    // genBestPHI_LowPt_Zinc2jet->Fill(PHI(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                    genSpTJets_LowPt_Zinc2jet->Fill(SpTsub(genLeadJ, genSecondJ), genWeight);
                    // genBestSpTJets_LowPt_Zinc2jet->Fill(SpTsub(genBestTwoJets.first, genBestTwoJets.second), genWeight);
                    genSpTLeptons_LowPt_Zinc2jet->Fill(SpTsub(genLep1, genLep2), genWeight);
                    genSpT_LowPt_Zinc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    // genBestSpT_LowPt_Zinc2jet->Fill(SpT(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                    genSPhi_LowPt_Zinc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    // genBestSPhi_LowPt_Zinc2jet->Fill(SPhi(genLep1, genLep2, genBestTwoJets.first, genBestTwoJets.second), genWeight);
                    if (SpT(genLep1, genLep2, genLeadJ, genSecondJ) < 0.5){ 
                        genPHI_LowSpT_LowPt_Zinc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genSPhi_LowSpT_LowPt_Zinc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    }
                    else {
                        genPHI_HighSpT_LowPt_Zinc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genSPhi_HighSpT_LowPt_Zinc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    }
                    if (SPhi(genLep1, genLep2, genLeadJ, genSecondJ) < 0.5){
                        genSpT_LowSPhi_LowPt_Zinc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    }
                    else {
                        genSpT_HighSPhi_LowPt_Zinc2jet ->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    }
                }
                else {
                    genptBal_HighPt_Zinc2jet->Fill(genJet1Plus2PlusZ.Pt(), genWeight);
                    gendPhiJets_HighPt_Zinc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                    gendPhiLeptons_HighPt_Zinc2jet->Fill(deltaPhi(genLep1, genLep2), genWeight);
                    genPHI_HighPt_Zinc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    genPHI_T_HighPt_Zinc2jet->Fill(PHI_T(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    genSpTJets_HighPt_Zinc2jet->Fill(SpTsub(genLeadJ, genSecondJ), genWeight);
                    genSpTLeptons_HighPt_Zinc2jet->Fill(SpTsub(genLep1, genLep2), genWeight);
                    genSpT_HighPt_Zinc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    genSPhi_HighPt_Zinc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                }
                if (nGoodGenJets == 2){
                    genTwoJetsPtDiff_Zexc2jet->Fill(genJet1Minus2.Pt(), genWeight);
                    genJetsMass_Zexc2jet->Fill(genJet1Plus2.M(), genWeight);
                    genSecondJetPt_Zexc2jet->Fill(genSecondJ.Pt(), genWeight);
                    gendPhiJets_Zexc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                    genPHI_Zexc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    genPHI_T_Zexc2jet->Fill(PHI_T(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    gendEtaJets_Zexc2jet->Fill(genLeadJ.Eta() - genSecondJ.Eta(), genWeight);
                    gendEtaFirstJetZ_Zexc2jet->Fill(genLeadJ.Eta() - genZ.Eta(), genWeight);
                    gendEtaSecondJetZ_Zexc2jet->Fill(genSecondJ.Eta() - genZ.Eta(), genWeight);
                    gendEtaJet1Plus2Z_Zexc2jet->Fill(genJet1Plus2.Eta() - genZ.Eta(), genWeight);
                    genSpT_Zexc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    genSpTJets_Zexc2jet->Fill(SpTsub(genLeadJ, genSecondJ), genWeight);
                    genSpTLeptons_Zexc2jet->Fill(SpTsub(genLep1, genLep2), genWeight);
                    genSPhi_Zexc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    genptBal_Zexc2jet->Fill(genJet1Plus2PlusZ.Pt(), genWeight);


                    if (genZ.Pt() < 25){
                        genptBal_LowPt_Zexc2jet->Fill(genJet1Plus2PlusZ.Pt(), genWeight);
                        gendPhiJets_LowPt_Zexc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                        gendPhiLeptons_LowPt_Zexc2jet->Fill(deltaPhi(genLep1, genLep2), genWeight);
                        genPHI_T_LowPt_Zexc2jet->Fill(PHI_T(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genPHI_LowPt_Zexc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genSpTJets_LowPt_Zexc2jet->Fill(SpTsub(genLeadJ, genSecondJ), genWeight);
                        genSpTLeptons_LowPt_Zexc2jet->Fill(SpTsub(genLep1, genLep2), genWeight);
                        genSpT_LowPt_Zexc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genSPhi_LowPt_Zexc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        if (SpT(genLep1, genLep2, genLeadJ, genSecondJ) < 0.5) { 
                            genPHI_LowSpT_LowPt_Zexc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            genSPhi_LowSpT_LowPt_Zexc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        }
                        else {
                            genPHI_HighSpT_LowPt_Zexc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                            genSPhi_HighSpT_LowPt_Zexc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        }
                        if (SPhi(genLep1, genLep2, genLeadJ, genSecondJ) < 0.5) {
                            genSpT_LowSPhi_LowPt_Zexc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        }
                        else {
                            genSpT_HighSPhi_LowPt_Zexc2jet ->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        }
                    }
                    else {
                        genptBal_HighPt_Zexc2jet->Fill(genJet1Plus2PlusZ.Pt(), genWeight);
                        gendPhiJets_HighPt_Zexc2jet->Fill(deltaPhi(genLeadJ, genSecondJ), genWeight);
                        gendPhiLeptons_HighPt_Zexc2jet->Fill(deltaPhi(genLep1, genLep2), genWeight);
                        genPHI_HighPt_Zexc2jet->Fill(PHI(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genPHI_T_HighPt_Zexc2jet->Fill(PHI_T(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genSpTJets_HighPt_Zexc2jet->Fill(SpTsub(genLeadJ, genSecondJ), genWeight);
                        genSpTLeptons_HighPt_Zexc2jet->Fill(SpTsub(genLep1, genLep2), genWeight);
                        genSpT_HighPt_Zexc2jet->Fill(SpT(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                        genSPhi_HighPt_Zexc2jet->Fill(SPhi(genLep1, genLep2, genLeadJ, genSecondJ), genWeight);
                    }
                }

            }
            if (nGoodGenJetsAK8 >= 2){
                genLeadingJetAK8Pt_Zinc2jet->Fill(genJetsAK8[0].pt, genWeight);
                genLeadingJetAK8Pt_2_Zinc2jet->Fill(genJetsAK8[0].pt, genWeight);
                genLepPtPlusLeadingJetAK8Pt_Zinc2jet->Fill(genLep1.Pt()+genJetsAK8[0].pt, genWeight);
                genLepPtPlusLeadingJetAK8Pt_2_Zinc2jet->Fill(genLep1.Pt()+genJetsAK8[0].pt, genWeight);

                if (nGoodGenJetsAK8 == 2){
                    genLeadingJetAK8Pt_Zexc2jet->Fill(genJetsAK8[0].pt, genWeight);
                    genLeadingJetAK8Pt_2_Zexc2jet->Fill(genJetsAK8[0].pt, genWeight);
                    genLepPtPlusLeadingJetAK8Pt_Zexc2jet->Fill(genLep1.Pt()+genJetsAK8[0].pt, genWeight);
                    genLepPtPlusLeadingJetAK8Pt_2_Zexc2jet->Fill(genLep1.Pt()+genJetsAK8[0].pt, genWeight);
                }
                
            }
            if (nGoodGenJets >= 3){
                GENnEventsIncl3Jets++;
                genZNGoodJets_Zinc->Fill(3., genWeight);
                genZNGoodJetsFull_Zinc->Fill(3., genWeight);
                genZMass_Zinc3jet->Fill(genZ.M(), genWeight);
                //  genZPt_Zinc3jet->Fill(genZ.Pt(), genWeight);
                genZRapidity_Zinc3jet->Fill(genZ.Rapidity(), genWeight);
                genZEta_Zinc3jet->Fill(genZ.Eta(), genWeight);

                genThirdJetEta_Zinc3jet->Fill(fabs(genJets[2].eta), genWeight);
                genThirdJetEta_2_Zinc3jet->Fill(fabs(genJets[2].eta), genWeight);
                genThirdJetAbsRapidity_Zinc3jet->Fill(fabs(genNewThirdJ.Rapidity()), genWeight);
                genThirdJetAbsRapidity_2_Zinc3jet->Fill(fabs(genNewThirdJ.Rapidity()), genWeight);
                genThirdJetRapidityFull_Zinc3jet->Fill(genNewThirdJ.Rapidity(), genWeight);
                gendRapidityJets_Zinc3jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), genWeight);
                gendRapidityJets_2_Zinc3jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), genWeight);
                
                gendRapidityJets_First_Third_Zinc3jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewThirdJ.Rapidity()), genWeight);
                gendRapidityJets_2_First_Third_Zinc3jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewThirdJ.Rapidity()), genWeight);
                
                gendRapidityJets_Second_Third_Zinc3jet->Fill(fabs(genNewSecondJ.Rapidity() - genNewThirdJ.Rapidity()), genWeight);
                gendRapidityJets_2_Second_Third_Zinc3jet->Fill(fabs(genNewSecondJ.Rapidity() - genNewThirdJ.Rapidity()), genWeight);
                
                gendRapidityJetsFB_Zinc3jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, genWeight);
                gendRapidityJetsFB_2_Zinc3jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, genWeight);
                
                gendiJetMass_Zinc3jet->Fill(genJet1Plus2.M(), genWeight);
                gendiJetMass_2_Zinc3jet->Fill(genJet1Plus2.M(), genWeight);
                
                gendiJetPt_Zinc3jet->Fill(genJet1Plus2.Pt(), genWeight);
                gendiJetPt_2_Zinc3jet->Fill(genJet1Plus2.Pt(), genWeight);

                gendPhiLepJet3_Zinc3jet->Fill(deltaPhi(genLep1, genNewThirdJ), genWeight);
                gendPhiLepJet3_2_Zinc3jet->Fill(deltaPhi(genLep1, genNewThirdJ), genWeight);
                genThirdJetPtEta_Zinc3jet->Fill(genJets[2].pt, fabs(genJets[2].eta), genWeight);
                genRatioJetPt32_Zinc3jet->Fill(genJets[2].pt/genJets[1].pt, genWeight);
                genFirstHighestJetPt_Zinc3jet->Fill(genJets[0].pt, genWeight);
                genSecondHighestJetPt_Zinc3jet->Fill(genJets[1].pt, genWeight);
                genThirdHighestJetPt_Zinc3jet->Fill(genJets[2].pt, genWeight);
                genJetsHT_Zinc3jet->Fill(genJetsHT, genWeight);
                genJetsHT_20_30_Zinc3jet->Fill(genJetsHT_20, genWeight);
                genJetsHT_1_Zinc3jet->Fill(genJetsHT, genWeight);
                genJetsHT_2_Zinc3jet->Fill(genJetsHT, genWeight);
                
            }
            if (nGoodGenJetsAK8 >= 3){
                genLeadingJetAK8Pt_Zinc3jet->Fill(genJetsAK8[0].pt, genWeight);
                genLeadingJetAK8Pt_2_Zinc3jet->Fill(genJetsAK8[0].pt, genWeight);
                genLepPtPlusLeadingJetAK8Pt_Zinc3jet->Fill(genLep1.Pt()+genJetsAK8[0].pt, genWeight);
                genLepPtPlusLeadingJetAK8Pt_2_Zinc3jet->Fill(genLep1.Pt()+genJetsAK8[0].pt, genWeight);

                if (nGoodGenJetsAK8 == 3){
                    genLeadingJetAK8Pt_Zexc3jet->Fill(genJetsAK8[0].pt, genWeight);
                    genLeadingJetAK8Pt_2_Zexc3jet->Fill(genJetsAK8[0].pt, genWeight);
                    genLepPtPlusLeadingJetAK8Pt_Zexc3jet->Fill(genLep1.Pt()+genJetsAK8[0].pt, genWeight);
                    genLepPtPlusLeadingJetAK8Pt_2_Zexc3jet->Fill(genLep1.Pt()+genJetsAK8[0].pt, genWeight);
                }
        
            }
            if (nGoodGenJets >= 4){
                genZNGoodJets_Zinc->Fill(4., genWeight);
                genZNGoodJetsFull_Zinc->Fill(4., genWeight);
                genZMass_Zinc4jet->Fill(genZ.M(), genWeight);
                //  genZPt_Zinc4jet->Fill(genZ.Pt(), genWeight);
                genZRapidity_Zinc4jet->Fill(genZ.Rapidity(), genWeight);
                genZEta_Zinc4jet->Fill(genZ.Eta(), genWeight);
                
                genFourthJetEta_Zinc4jet->Fill(fabs(genJets[3].eta), genWeight);
                genFourthJetEta_2_Zinc4jet->Fill(fabs(genJets[3].eta), genWeight);
                genFourthJetAbsRapidity_Zinc4jet->Fill(fabs(genNewFourthJ.Rapidity()), genWeight);
                genFourthJetAbsRapidity_2_Zinc4jet->Fill(fabs(genNewFourthJ.Rapidity()), genWeight);
                genFourthJetRapidityFull_Zinc4jet->Fill(genNewFourthJ.Rapidity(), genWeight);
                gendRapidityJets_Zinc4jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), genWeight);
                gendRapidityJets_2_Zinc4jet->Fill(fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), genWeight);
                
                gendRapidityJetsFB_Zinc4jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, genWeight);
                gendRapidityJetsFB_2_Zinc4jet->Fill(genForwardJetRapidity - genBackwardJetRapidity, genWeight);
                
                gendiJetMass_Zinc4jet->Fill(genJet1Plus2.M(), genWeight);
                gendiJetMass_2_Zinc4jet->Fill(genJet1Plus2.M(), genWeight);
                
                gendiJetPt_Zinc4jet->Fill(genJet1Plus2.Pt(), genWeight);
                gendiJetPt_2_Zinc4jet->Fill(genJet1Plus2.Pt(), genWeight);

                gendPhiLepJet4_Zinc4jet->Fill(deltaPhi(genLep1, genNewFourthJ), genWeight);
                gendPhiLepJet4_2_Zinc4jet->Fill(deltaPhi(genLep1, genNewFourthJ), genWeight);

                genFourthJetPtEta_Zinc4jet->Fill(genJets[3].pt, fabs(genJets[3].eta), genWeight);
                genFirstHighestJetPt_Zinc4jet->Fill(genJets[0].pt, genWeight);
                genSecondHighestJetPt_Zinc4jet->Fill(genJets[1].pt, genWeight);
                genThirdHighestJetPt_Zinc4jet->Fill(genJets[2].pt, genWeight);
                genJetsHT_Zinc4jet->Fill(genJetsHT, genWeight);
                genJetsHT_1_Zinc4jet->Fill(genJetsHT, genWeight);
                genJetsHT_2_Zinc4jet->Fill(genJetsHT, genWeight);


            }
            if (nGoodGenJets >= 5){
                genZNGoodJets_Zinc->Fill(5., genWeight);
                genZNGoodJetsFull_Zinc->Fill(5., genWeight);
                genZMass_Zinc5jet->Fill(genZ.M(), genWeight);
                genZPt_Zinc5jet->Fill(genZ.Pt(), genWeight);
                genZRapidity_Zinc5jet->Fill(genZ.Rapidity(), genWeight);
                genZEta_Zinc5jet->Fill(genZ.Eta(), genWeight);
                
                genFifthJetEta_Zinc5jet->Fill(fabs(genJets[4].eta), genWeight);
                genFifthJetEta_2_Zinc5jet->Fill(fabs(genJets[4].eta), genWeight);
                genFifthJetPtEta_Zinc5jet->Fill(genJets[4].pt, fabs(genJets[4].eta), genWeight);
                genFirstHighestJetPt_Zinc5jet->Fill(genJets[0].pt, genWeight);
                genSecondHighestJetPt_Zinc5jet->Fill(genJets[1].pt, genWeight);
                genThirdHighestJetPt_Zinc5jet->Fill(genJets[2].pt, genWeight);
                genJetsHT_Zinc5jet->Fill(genJetsHT, genWeight);
                genJetsHT_1_Zinc5jet->Fill(genJetsHT, genWeight);
                genJetsHT_2_Zinc5jet->Fill(genJetsHT, genWeight);

                //gendPhiLepJet5_Zinc5jet->Fill(deltaPhi(genLep1, genNewFifthJ), genWeight);
                //gendPhiLepJet5_2_Zinc5jet->Fill(deltaPhi(genLep1, genNewFifthJ), genWeight);
            }
            if (nGoodGenJets >= 6){
                genZNGoodJets_Zinc->Fill(6., genWeight);
                genZNGoodJetsFull_Zinc->Fill(6., genWeight);
                genZMass_Zinc6jet->Fill(genZ.M(), genWeight);
                genZPt_Zinc6jet->Fill(genZ.Pt(), genWeight);
                genZRapidity_Zinc6jet->Fill(genZ.Rapidity(), genWeight);
                genZEta_Zinc6jet->Fill(genZ.Eta(), genWeight);
                
                genSixthJetEta_Zinc6jet->Fill(fabs(genJets[5].eta), genWeight);
                genSixthJetPtEta_Zinc6jet->Fill(genJets[5].pt, fabs(genJets[5].eta), genWeight);
                genFirstHighestJetPt_Zinc6jet->Fill(genJets[0].pt, genWeight);
                genSecondHighestJetPt_Zinc6jet->Fill(genJets[1].pt, genWeight);
                genThirdHighestJetPt_Zinc6jet->Fill(genJets[2].pt, genWeight);
                genJetsHT_Zinc6jet->Fill(genJetsHT, genWeight);
                genJetsHT_1_Zinc6jet->Fill(genJetsHT, genWeight);
            }
            if (nGoodGenJets >= 7){
                genZNGoodJets_Zinc->Fill(7., genWeight);
                genZNGoodJetsFull_Zinc->Fill(7., genWeight);
            }
            if (nGoodGenJets >= 8) genZNGoodJets_Zinc->Fill(8., genWeight);
            if (nGoodGenJets >= 9) genZNGoodJets_Zinc->Fill(9., genWeight);
            if (nGoodGenJets >= 10) genZNGoodJets_Zinc->Fill(10., genWeight);
            
        } //end filling gen histos
        
        //=======================================================================================================//

        //=======================================================================================================//
        //      Selection for Reco Histos      //
        //====================================//
        if (hasRecoInfo && PRINTEVENTINFO && jentry == eventOfInterest) {
            cout << __LINE__ << " PRINTEVENTINFO: Requirements for filling reco histos --- " << endl;
            cout << __LINE__ << " PRINTEVENTINFO: hasRecoInfo = " << hasRecoInfo << endl;
            cout << __LINE__ << " PRINTEVENTINFO: passesLeptonCut = " << passesLeptonCut << endl;
            cout << __LINE__ << " PRINTEVENTINFO: passesJetCut = " << passesJetCut << "\n" << endl;
        }
        if (hasRecoInfo && passesLeptonCut && passesJetCut) {
            //=======================================================================================================//
            
            if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
            //=======================================================================================================//
            //      Start filling histograms      //
            //====================================//
            
            countEventpassBveto++;
            TotalRecoWeightPassRECO+=weight;
            TotalGenWeightPassRECO+=genWeightBackup;
            NVtx->Fill(EvtVtxCnt, weight);

            PUWeight->Fill(puWeight.weight(int(EvtPuCntObs)), 1);
            if (nGoodJets == 0){
                PUWeight0->Fill(puWeight.weight(int(EvtPuCntObs)), 1);
            }
            else {
                PUWeight1->Fill(puWeight.weight(int(EvtPuCntObs)), 1);
            }
            
            if (lepton1.charge > 0){
                MuPlusPt->Fill(lepton1.pt, weight);
                MuPlusEta->Fill(lepton1.eta, weight);
                if (doZ || doTT){
                    MuMinusPt->Fill(lepton2.pt, weight);
                    MuMinusEta->Fill(lepton2.eta, weight);
                }
            }
            else {
                MuMinusPt->Fill(lepton1.pt, weight);
                MuMinusEta->Fill(lepton1.eta, weight);
                if (doZ || doTT){
                    MuPlusPt->Fill(lepton2.pt, weight);
                    MuPlusEta->Fill(lepton2.eta, weight);
                }
            }
            nEventsIncl0Jets++;
            ZNGoodJetsNVtx_Zexc->Fill(nGoodJets, EvtVtxCnt, weight);
            ZNGoodJets_Zinc->Fill(0., weight);
            ZNGoodJetsFull_Zinc->Fill(0., weight);
            ZNGoodJets_Zexc->Fill(nGoodJets, weight);
            ZNGoodJetsFull_Zexc->Fill(nGoodJets, weight);
            ZNGoodJets_Zinc_NoWeight->Fill(0.);
            ZMass_Zinc0jet->Fill(Z.M(), weight);
            MET_Zinc0jet->Fill(METpt, weight);
            //MET_1_Zinc0jet->Fill(METpt, weight);
            METphi_Zinc0jet->Fill(METphi, weight);
            MT_Zinc0jet->Fill(MT, weight);
            
            MuPFIso_Zinc0jet->Fill(lepton1.iso, weight);
            MuPFIso_2ndZinc0jet->Fill(lepton1.iso, weight);
            MuPFIso_3rdZinc0jet->Fill(lepton1.iso, weight);
            
            ZPt_Zinc0jet->Fill(Z.Pt(), weight);
            ZRapidity_Zinc0jet->Fill(Z.Rapidity(), weight);
            ZEta_Zinc0jet->Fill(Z.Eta(), weight);
            lepPt_Zinc0jet->Fill(lepton1.pt, weight);
            lepEta_Zinc0jet->Fill(lepton1.eta, weight);
            lepPhi_Zinc0jet->Fill(lepton1.phi, weight);
            if (doZ || doTT){
                lepPt_Zinc0jet->Fill(lepton2.pt, weight);
                lepEta_Zinc0jet->Fill(lepton2.eta, weight);
                lepPhi_Zinc0jet->Fill(lepton2.phi, weight);
            }
            lepEtaEta_Zinc0jet->Fill(lepton1.eta, lepton2.eta, weight);
            dPhiLeptons_Zinc0jet->Fill(deltaPhi(lep1, lep2), weight);
            dEtaLeptons_Zinc0jet->Fill(lepton1.eta - lepton2.eta, weight);
            dRLeptons_Zinc0jet->Fill(deltaR(lepton1.phi, lepton1.eta, lepton2.phi, lepton2.eta), weight);
            SpTLeptons_Zinc0jet->Fill(SpTsub(lep1, lep2), weight);
            if (nGoodJets == 0){
                nEventsExcl0Jets++;
                ZNGoodJets_Zexc_NoWeight->Fill(0.);
                ZMass_Zexc0jet->Fill(Z.M(), weight);
                ZPt_Zexc0jet->Fill(Z.Pt(), weight);
                ZRapidity_Zexc0jet->Fill(Z.Rapidity(), weight);
                ZEta_Zexc0jet->Fill(Z.Eta(), weight);
                lepPt_Zexc0jet->Fill(lepton1.pt, weight);
                lepEta_Zexc0jet->Fill(lepton1.eta, weight);
                if ( doZ || doTT){
                    lepPt_Zexc0jet->Fill(lepton2.pt, weight);
                    lepEta_Zexc0jet->Fill(lepton2.eta, weight);
                }
                dPhiLeptons_Zexc0jet->Fill(deltaPhi(lep1, lep2), weight);
                dEtaLeptons_Zexc0jet->Fill(lepton1.eta - lepton2.eta, weight);
                SpTLeptons_Zexc0jet->Fill(SpTsub(lep1, lep2), weight);
            }
            
            if (nGoodJets_20 >= 1){
                FirstJetPt_Zinc1jet  ->Fill(jets_20[0].pt, weight);
                FirstJetPt_1_Zinc1jet->Fill(jets_20[0].pt, weight);
                FirstJetPt_2_Zinc1jet->Fill(jets_20[0].pt, weight);
                JetsHT_20_Zinc1jet->Fill(jetsHT_20, weight);

                //andrew - alpha-s
                LeadingJetPt_Zinc1jet->Fill(jets_20[0].pt, weight);
                // LeadingJetPt_1_Zinc1jet->Fill(jets_20[0].pt, weight);
                LeadingJetPt_2_Zinc1jet->Fill(jets_20[0].pt, weight);
                LepPtPlusLeadingJetPt_Zinc1jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                // LepPtPlusLeadingJetPt_1_Zinc1jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                LepPtPlusLeadingJetPt_2_Zinc1jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                ZPt_Zinc1jet->Fill(Z.Pt(), weight);
                // ZPt_1_Zinc1jet->Fill(Z.Pt(), weight);
                ZPt_2_Zinc1jet->Fill(Z.Pt(), weight);
                ZPtPlusLeadingJetPt_Zinc1jet->Fill(Z.Pt()+jets_20[0].pt, weight);
                // ZPtPlusLeadingJetPt_1_Zinc1jet->Fill(Z.Pt()+jets_20[0].pt, weight);
                ZPtPlusLeadingJetPt_2_Zinc1jet->Fill(Z.Pt()+jets_20[0].pt, weight);

                if (nGoodJets_20 == 1){
                    LeadingJetPt_Zexc1jet->Fill(jets_20[0].pt, weight);
                    LeadingJetPt_2_Zexc1jet->Fill(jets_20[0].pt, weight);
                    LepPtPlusLeadingJetPt_Zexc1jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                    LepPtPlusLeadingJetPt_2_Zexc1jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                }

            }
            if (nGoodJets_20 >= 2){
                SecondJetPt_Zinc2jet  ->Fill(jets_20[1].pt, weight);
                SecondJetPt_1_Zinc2jet->Fill(jets_20[1].pt, weight);
                SecondJetPt_2_Zinc2jet->Fill(jets_20[1].pt, weight);
                JetsHT_20_Zinc2jet->Fill(jetsHT_20, weight);

                //andrew - alpha-s
                LeadingJetPt_Zinc2jet->Fill(jets_20[0].pt, weight);
                // LeadingJetPt_1_Zinc2jet->Fill(jets_20[0].pt, weight);
                LeadingJetPt_2_Zinc2jet->Fill(jets_20[0].pt, weight);
                HTover2_Zinc2jet->Fill((jets_20[0].pt + jets_20[1].pt)/2. , weight);
                // HTover2_1_Zinc2jet->Fill((jets_20[0].pt + jets_20[1].pt)/2. , weight);
                HTover2_2_Zinc2jet->Fill((jets_20[0].pt + jets_20[1].pt)/2. , weight);
                LepPtPlusLeadingJetPt_Zinc2jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                // LepPtPlusLeadingJetPt_1_Zinc2jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                LepPtPlusLeadingJetPt_2_Zinc2jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                LepPtPlusHTover2_Zinc2jet->Fill(lepton1.pt + (jets_20[0].pt + jets_20[1].pt)/2., weight);
                // LepPtPlusHTover2_1_Zinc2jet->Fill(lepton1.pt + (jets_20[0].pt + jets_20[1].pt)/2., weight);
                LepPtPlusHTover2_2_Zinc2jet->Fill(lepton1.pt + (jets_20[0].pt + jets_20[1].pt)/2., weight);
                ZPt_Zinc2jet->Fill(Z.Pt(), weight);
                // ZPt_1_Zinc2jet->Fill(Z.Pt(), weight);
                ZPt_2_Zinc2jet->Fill(Z.Pt(), weight);
                ZPtPlusLeadingJetPt_Zinc2jet->Fill(Z.Pt()+jets_20[0].pt, weight);
                // ZPtPlusLeadingJetPt_1_Zinc2jet->Fill(Z.Pt()+jets_20[0].pt, weight);
                ZPtPlusLeadingJetPt_2_Zinc2jet->Fill(Z.Pt()+jets_20[0].pt, weight);
                ZPtPlusHTover2_Zinc2jet->Fill(Z.Pt()+ (jets_20[0].pt + jets_20[1].pt)/2., weight);
                // ZPtPlusHTover2_1_Zinc2jet->Fill(Z.Pt()+ (jets_20[0].pt + jets_20[1].pt)/2., weight);
                ZPtPlusHTover2_2_Zinc2jet->Fill(Z.Pt()+ (jets_20[0].pt + jets_20[1].pt)/2., weight);

                if (nGoodJets_20 == 2){
                    LeadingJetPt_Zexc2jet->Fill(jets_20[0].pt, weight);
                    LeadingJetPt_2_Zexc2jet->Fill(jets_20[0].pt, weight);
                    LepPtPlusLeadingJetPt_Zexc2jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                    LepPtPlusLeadingJetPt_2_Zexc2jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                }
            }
            if (nGoodJets_20 >= 3){
                ThirdJetPt_Zinc3jet  ->Fill(jets_20[2].pt, weight);
                ThirdJetPt_1_Zinc3jet->Fill(jets_20[2].pt, weight);
                ThirdJetPt_2_Zinc3jet->Fill(jets_20[2].pt, weight);
                JetsHT_20_Zinc3jet->Fill(jetsHT_20, weight);

                //andrew - alpha-s
                LeadingJetPt_Zinc3jet->Fill(jets_20[0].pt, weight);
                // LeadingJetPt_1_Zinc3jet->Fill(jets_20[0].pt, weight);
                LeadingJetPt_2_Zinc3jet->Fill(jets_20[0].pt, weight);
                HTover2_Zinc3jet->Fill((jets_20[0].pt+jets_20[1].pt)/2. , weight);
                // HTover2_1_Zinc3jet->Fill((jets_20[0].pt+jets_20[1].pt)/2. , weight);
                HTover2_2_Zinc3jet->Fill((jets_20[0].pt+jets_20[1].pt)/2. , weight);
                LepPtPlusLeadingJetPt_Zinc3jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                // LepPtPlusLeadingJetPt_1_Zinc3jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                LepPtPlusLeadingJetPt_2_Zinc3jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                LepPtPlusHTover2_Zinc3jet->Fill(lepton1.pt + (jets_20[0].pt+jets_20[1].pt)/2., weight);
                // LepPtPlusHTover2_1_Zinc3jet->Fill(lepton1.pt + (jets_20[0].pt+jets_20[1].pt)/2., weight);
                LepPtPlusHTover2_2_Zinc3jet->Fill(lepton1.pt + (jets_20[0].pt+jets_20[1].pt)/2., weight);
                ZPt_Zinc3jet->Fill(Z.Pt(), weight);
                // ZPt_1_Zinc3jet->Fill(Z.Pt(), weight);
                ZPt_2_Zinc3jet->Fill(Z.Pt(), weight);
                ZPtPlusLeadingJetPt_Zinc3jet->Fill(Z.Pt()+jets_20[0].pt, weight);
                // ZPtPlusLeadingJetPt_1_Zinc3jet->Fill(Z.Pt()+jets_20[0].pt, weight);
                ZPtPlusLeadingJetPt_2_Zinc3jet->Fill(Z.Pt()+jets_20[0].pt, weight);
                ZPtPlusHTover2_Zinc3jet->Fill(Z.Pt()+ (jets_20[0].pt+jets_20[1].pt)/2., weight);
                // ZPtPlusHTover2_1_Zinc3jet->Fill(Z.Pt()+ (jets_20[0].pt+jets_20[1].pt)/2., weight);
                ZPtPlusHTover2_2_Zinc3jet->Fill(Z.Pt()+ (jets_20[0].pt+jets_20[1].pt)/2., weight);

                if (nGoodJets_20 == 3){
                    LeadingJetPt_Zexc3jet->Fill(jets_20[0].pt, weight);
                    LeadingJetPt_2_Zexc3jet->Fill(jets_20[0].pt, weight);
                    LepPtPlusLeadingJetPt_Zexc3jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                    LepPtPlusLeadingJetPt_2_Zexc3jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                }
            }
            if (nGoodJets_20 >= 4){
                FourthJetPt_Zinc4jet  ->Fill(jets_20[3].pt, weight);
                FourthJetPt_1_Zinc4jet->Fill(jets_20[3].pt, weight);
                FourthJetPt_2_Zinc4jet->Fill(jets_20[3].pt, weight);

                //andrew - alpha-s
                LeadingJetPt_Zinc4jet->Fill(jets_20[0].pt, weight);
                // LeadingJetPt_1_Zinc4jet->Fill(jets_20[0].pt, weight);
                LeadingJetPt_2_Zinc4jet->Fill(jets_20[0].pt, weight);
                HTover2_Zinc4jet->Fill((jets_20[0].pt + jets_20[1].pt)/2. , weight);
                // HTover2_1_Zinc4jet->Fill((jets_20[0].pt + jets_20[1].pt)/2. , weight);
                HTover2_2_Zinc4jet->Fill((jets_20[0].pt + jets_20[1].pt)/2. , weight);
                LepPtPlusLeadingJetPt_Zinc4jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                // LepPtPlusLeadingJetPt_1_Zinc4jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                LepPtPlusLeadingJetPt_2_Zinc4jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                LepPtPlusHTover2_Zinc4jet->Fill(lepton1.pt + (jets_20[0].pt + jets_20[1].pt)/2., weight);
                // LepPtPlusHTover2_1_Zinc4jet->Fill(lepton1.pt + (jets_20[0].pt + jets_20[1].pt)/2., weight);
                LepPtPlusHTover2_2_Zinc4jet->Fill(lepton1.pt + (jets_20[0].pt + jets_20[1].pt)/2., weight);
                ZPt_Zinc4jet->Fill(Z.Pt(), weight);
                // ZPt_1_Zinc4jet->Fill(Z.Pt(), weight);
                ZPt_2_Zinc4jet->Fill(Z.Pt(), weight);
                ZPtPlusLeadingJetPt_Zinc4jet->Fill(Z.Pt()+jets_20[0].pt, weight);
                // ZPtPlusLeadingJetPt_1_Zinc4jet->Fill(Z.Pt()+jets_20[0].pt, weight);
                ZPtPlusLeadingJetPt_2_Zinc4jet->Fill(Z.Pt()+jets_20[0].pt, weight);
                ZPtPlusHTover2_Zinc4jet->Fill(Z.Pt()+ (jets_20[0].pt + jets_20[1].pt)/2., weight);
                // ZPtPlusHTover2_1_Zinc4jet->Fill(Z.Pt()+ (jets_20[0].pt + jets_20[1].pt)/2., weight);
                ZPtPlusHTover2_2_Zinc4jet->Fill(Z.Pt()+ (jets_20[0].pt + jets_20[1].pt)/2., weight);

                // if (nGoodJets_20 == 4){
                //     LeadingJetPt_Zexc4jet->Fill(jets_20[0].pt, weight);
                //     LeadingJetPt_2_Zexc4jet->Fill(jets_20[0].pt, weight);
                //     LepPtPlusLeadingJetPt_Zexc4jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                //     LepPtPlusLeadingJetPt_2_Zexc4jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                // }
            }
            if (nGoodJets_20 >= 5){
                FifthJetPt_Zinc5jet  ->Fill(jets_20[4].pt, weight);
                FifthJetPt_1_Zinc5jet->Fill(jets_20[4].pt, weight);
                FifthJetPt_2_Zinc5jet->Fill(jets_20[4].pt, weight);
            }
            if (nGoodJets_20 >= 6){
                SixthJetPt_Zinc6jet  ->Fill(jets_20[5].pt, weight);
                SixthJetPt_1_Zinc6jet->Fill(jets_20[5].pt, weight);
            }	
			if (nJetsDR02Cut >= 1){
				dRLepCloseJet_Zinc1jet->Fill(mindRjmu, weight);
				dRLepCloseJet_2_Zinc1jet->Fill(mindRjmu, weight);
				if(jetsDR02[0].pt > 500){
					dRLepCloseJetCo_Zinc1jet->Fill(mindRjmu, weight);
					dRLepCloseJetCo_2_Zinc1jet->Fill(mindRjmu, weight);
				}
			}
			if (nJetsPt100DR04 >= 1 && jetsPt100DR04[0].pt > 300){
				dRptmin100LepCloseJetCo300dR04_Zinc1jet->Fill(mindRj04Pt100mu, weight);
				dRptmin100LepCloseJetCo300dR04_2_Zinc1jet->Fill(mindRj04Pt100mu, weight);
			}
			if (nJetsPt100DR04 >= 2 && jetsPt100DR04[0].pt > 300 && deltaPhi(leadJ, secondJ) > (PI - 0.3)){
				dRptmin100LepCloseDiJetCo300dR04_Zinc2jet->Fill(mindRdijet04Pt100mu, weight);
				dRptmin100LepCloseDiJetCo300dR04_2_Zinc2jet->Fill(mindRdijet04Pt100mu, weight);
			}
      
            if (nGoodJets >= 1){
                ZNGoodJets_Zinc->Fill(1., weight);
                ZNGoodJetsFull_Zinc->Fill(1., weight);
                ZNGoodJets_Zinc_NoWeight->Fill(1.);

                ZMass_Zinc1jet->Fill(Z.M(), weight);
                MET_Zinc1jet->Fill(METpt, weight);
                //MET_1_Zinc1jet->Fill(METpt, weight);
                METphi_Zinc1jet->Fill(METphi, weight);
                MT_Zinc1jet->Fill(MT, weight);
                //                ZPt_Zinc1jet->Fill(Z.Pt(), weight);
                ZRapidity_Zinc1jet->Fill(Z.Rapidity(), weight);
                ZEta_Zinc1jet->Fill(Z.Eta(), weight);

                lepPt_Zinc1jet->Fill(lepton1.pt, weight);
                lepEta_Zinc1jet->Fill(lepton1.eta, weight);
                lepPhi_Zinc1jet->Fill(lepton1.phi, weight);
                if (lepton1.charge > 0){
                    lepChargePlusPt_Zinc1jet->Fill(lepton1.pt, weight);
                    lepChargePlusEta_Zinc1jet->Fill(lepton1.eta, weight);
                    lepChargePlusPhi_Zinc1jet->Fill(lepton1.phi, weight);
                }
                if (lepton1.charge < 0){
                    lepChargeMinusPt_Zinc1jet->Fill(lepton1.pt, weight);
                    lepChargeMinusEta_Zinc1jet->Fill(lepton1.eta, weight);
                    lepChargeMinusPhi_Zinc1jet->Fill(lepton1.phi, weight);
                }
                if (doZ || doTT){
                    lepPt_Zinc1jet->Fill(lepton2.pt, weight);
                    lepEta_Zinc1jet->Fill(lepton2.eta, weight);
                    lepPhi_Zinc1jet->Fill(lepton2.phi, weight);
                }

                dPhiLeptons_Zinc1jet->Fill(deltaPhi(lep1, lep2), weight);
                dEtaLeptons_Zinc1jet->Fill(lepton1.eta - lepton2.eta, weight);
                dRLeptons_Zinc1jet->Fill(deltaR(lepton1.phi, lepton1.eta, lepton2.phi, lepton2.eta), weight);
                SpTLeptons_Zinc1jet->Fill(SpTsub(lep1, lep2), weight);
                
                FirstHighestJetPt_Zinc1jet->Fill(jets[0].pt, weight);
                FirstJetEta_Zinc1jet->Fill(fabs(jets[0].eta), weight);
                FirstJetEta_2_Zinc1jet->Fill(fabs(jets[0].eta), weight);
                FirstJetEtaFull_Zinc1jet->Fill(jets[0].eta, weight);
                FirstJetAbsRapidity_Zinc1jet->Fill(fabs(newLeadJ.Rapidity()), weight);
                FirstJetAbsRapidity_2_Zinc1jet->Fill(fabs(newLeadJ.Rapidity()), weight);
                FirstJetRapidityFull_Zinc1jet->Fill(newLeadJ.Rapidity(), weight);
                FirstJetmass_Zinc1jet->Fill(newLeadJ.M(), weight);
                FirstJetmass_1_Zinc1jet->Fill(newLeadJ.M(), weight);
                MeanNJetsHT_1D_Zinc1jet->Fill(jetsHT, weight*nGoodJets);
                MeanNJetsHT_Zinc1jet->Fill(jetsHT, nGoodJets, weight);
                FirstJetPtEta_Zinc1jet->Fill(jets[0].pt, fabs(jets[0].eta), weight);
                FirstJetPhi_Zinc1jet->Fill(jets[0].phi, weight);
                JetsHT_Zinc1jet->Fill(jetsHT, weight);
                JetsHT_20_30_Zinc1jet->Fill(jetsHT_20, weight);
                JetsHT_1_Zinc1jet->Fill(jetsHT, weight);
                JetsHT_2_Zinc1jet->Fill(jetsHT, weight);
                
				dPhiLepJet1_Zinc1jet->Fill(deltaPhi(lep1, newLeadJ), weight);
				dPhiLepJet1_2_Zinc1jet->Fill(deltaPhi(lep1, newLeadJ), weight);
                
                for (unsigned short j(0); j < nGoodJets; j++){
                    AllJetPt_Zinc1jet->Fill(jets[j].pt, weight);
                    AllJetEta_Zinc1jet->Fill(jets[j].eta, weight);
                    AllJetPhi_Zinc1jet->Fill(jets[j].phi, weight);
                }

                if ( doW ) dEtaBosonJet_Zinc1jet->Fill(fabs(jets[0].eta - lepton1.eta), weight);
                else dEtaBosonJet_Zinc1jet->Fill(fabs(jets[0].eta-Z.Eta()), weight);
                for ( int i =0 ; i < NbinsEta2D - 1 ; i++){
                    if ( fabs(jets[0].eta) >= j_Y_range[i] &&  fabs(jets[0].eta) < j_Y_range[i+1] ) FirstJetPt_Zinc1jet_Eta[i]->Fill(fabs(jets[0].pt), weight);
                }

                if (nGoodJets == 1){
                    // compute Delta pt between Z and jets
                    if (Z.Pt() > 0.8 * jetPtCutMin  && jets[0].pt/Z.Pt() < 1.2 && jets[0].pt/Z.Pt() > 0.8 && deltaPhi(leadJ, Z) > 2.7){
                        hPtEtaBackJet_Zexc1jet->Fill(leadJ.Pt(), leadJ.Eta(), weight);
                        if (JetAk04PuMva->at(jets[0].patIndex) > 0 ) hPtEtaBackJetMVA_Zexc1jet->Fill(leadJ.Pt(), leadJ.Eta(), weight);
                    }
           
                    nEventsExcl1Jets++;
                    ZNGoodJets_Zexc_NoWeight->Fill(1.);
                    ZMass_Zexc1jet->Fill(Z.M(), weight);
                    ZPt_Zexc1jet->Fill(Z.Pt(), weight);
                    ZRapidity_Zexc1jet->Fill(Z.Rapidity(), weight);
                    ZEta_Zexc1jet->Fill(Z.Eta(), weight);
                    lepPt_Zexc1jet->Fill(lepton1.pt, weight);
                    lepEta_Zexc1jet->Fill(lepton1.eta, weight);
                    if (doZ || doTT){
                        lepPt_Zexc1jet->Fill(lepton2.pt, weight);
                        lepEta_Zexc1jet->Fill(lepton2.eta, weight);
                    }
                    dPhiLeptons_Zexc1jet->Fill(deltaPhi(lep1, lep2), weight);
                    dEtaLeptons_Zexc1jet->Fill(lepton1.eta - lepton2.eta, weight);
                    SpTLeptons_Zexc1jet->Fill(SpTsub(lep1, lep2), weight);
                    FirstJetPt_Zexc1jet->Fill(jets[0].pt, weight);
                    FirstJetEta_Zexc1jet->Fill(jets[0].eta, weight);
                    FirstJetPhi_Zexc1jet->Fill(jets[0].phi, weight);
                    if ( doW ) dEtaBosonJet_Zexc1jet->Fill(fabs(jets[0].eta - lepton1.eta), weight);
                    else dEtaBosonJet_Zexc1jet->Fill(fabs(jets[0].eta-Z.Eta()), weight);
                    
                }
            }
            if (nGoodJetsAK8 >= 1){
                for (unsigned short j(0); j < nGoodJetsAK8; j++){
                    AllJetAK8Pt_Zinc1jet->Fill(jetsAK8[j].pt, weight);
                    AllJetAK8Eta_Zinc1jet->Fill(jetsAK8[j].eta, weight);
                    AllJetAK8Phi_Zinc1jet->Fill(jetsAK8[j].phi, weight);
                }
                LeadingJetAK8Pt_Zinc1jet->Fill(jetsAK8[0].pt, weight);
                LeadingJetAK8Pt_2_Zinc1jet->Fill(jetsAK8[0].pt, weight);
                LepPtPlusLeadingJetAK8Pt_Zinc1jet->Fill(lepton1.pt + jetsAK8[0].pt, weight);
                LepPtPlusLeadingJetAK8Pt_2_Zinc1jet->Fill(lepton1.pt + jetsAK8[0].pt, weight);
                if (nGoodJetsAK8 == 1){
                    LeadingJetAK8Pt_Zexc1jet->Fill(jetsAK8[0].pt, weight);
                    LeadingJetAK8Pt_2_Zexc1jet->Fill(jetsAK8[0].pt, weight);
                    LepPtPlusLeadingJetAK8Pt_Zexc1jet->Fill(lepton1.pt + jetsAK8[0].pt, weight);
                    LepPtPlusLeadingJetAK8Pt_2_Zexc1jet->Fill(lepton1.pt + jetsAK8[0].pt, weight);
                }
            }
            if (nGoodJets >= 2){
                TLorentzVector jet1Plus2PlusZ = jet1Plus2 + Z;
                ZNGoodJets_Zinc->Fill(2., weight);
                ZNGoodJetsFull_Zinc->Fill(2., weight);
                ZNGoodJets_Zinc_NoWeight->Fill(2.);
                ZMass_Zinc2jet->Fill(Z.M(), weight);
                MET_Zinc2jet->Fill(METpt, weight);
                METphi_Zinc2jet->Fill(METphi, weight);
                MT_Zinc2jet->Fill(MT, weight);
                TwoJetsPtDiff_Zinc2jet->Fill(jet1Minus2.Pt(), weight);
                // BestTwoJetsPtDiff_Zinc2jet->Fill(bestJet1Minus2.Pt(), weight);
                JetsMass_Zinc2jet->Fill(jet1Plus2.M(), weight);
                //                ZPt_Zinc2jet->Fill(Z.Pt(), weight);
                ZRapidity_Zinc2jet->Fill(Z.Rapidity(), weight);
                ZEta_Zinc2jet->Fill(Z.Eta(), weight);
                lepPt_Zinc2jet->Fill(lepton1.pt, weight);
                lepEta_Zinc2jet->Fill(lepton1.eta, weight);
                if (doZ || doTT){
                    lepPt_Zinc2jet->Fill(lepton2.pt, weight);
                    lepEta_Zinc2jet->Fill(lepton2.eta, weight);
                }
                dPhiLeptons_Zinc2jet->Fill(deltaPhi(lep1, lep2), weight);
                dEtaLeptons_Zinc2jet->Fill(lepton1.eta - lepton2.eta, weight);
                dRLeptons_Zinc2jet->Fill(deltaR(lepton1.phi, lepton1.eta, lepton2.phi, lepton2.eta), weight);
                SpTLeptons_Zinc2jet->Fill(SpTsub(lep1, lep2), weight);
                FirstHighestJetPt_Zinc2jet->Fill(jets[0].pt, weight);
                SecondHighestJetPt_Zinc2jet->Fill(jets[1].pt, weight);
                
                SecondJetEta_Zinc2jet->Fill(fabs(jets[1].eta), weight);
                SecondJetEta_2_Zinc2jet->Fill(fabs(jets[1].eta), weight);
                
                SecondJetEtaFull_Zinc2jet->Fill(jets[1].eta, weight);
                SecondJetAbsRapidity_Zinc2jet->Fill(fabs(newSecondJ.Rapidity()), weight);
                SecondJetAbsRapidity_2_Zinc2jet->Fill(fabs(newSecondJ.Rapidity()), weight);
                SecondJetRapidityFull_Zinc2jet->Fill(newSecondJ.Rapidity(), weight);
                SecondJetmass_Zinc2jet->Fill(newSecondJ.M(), weight);
                SecondJetmass_1_Zinc2jet->Fill(newSecondJ.M(), weight);
                dRapidityJets_Zinc2jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), weight);
                dRapidityJets_2_Zinc2jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), weight);
                
                dRapidityJetsFB_Zinc2jet->Fill(ForwardJetRapidity - BackwardJetRapidity, weight);
                dRapidityJetsFB_2_Zinc2jet->Fill(ForwardJetRapidity - BackwardJetRapidity, weight);
                
                diJetMass_Zinc2jet->Fill(jet1Plus2.M(), weight);
                diJetMass_2_Zinc2jet->Fill(jet1Plus2.M(), weight);
                
                dPhiJetsFB_Zinc2jet->Fill(deltaPhi(vJetYOrdered[0], vJetYOrdered[vJetYOrdered.size() - 1] ), weight);
                dPhiJetsFB_2_Zinc2jet->Fill(deltaPhi(vJetYOrdered[0], vJetYOrdered[vJetYOrdered.size() - 1] ), weight);
                
                dRJets_Zinc2jet->Fill(deltaRYPhi(newLeadJ, newSecondJ), weight);
                dRJets_2_Zinc2jet->Fill(deltaRYPhi(newLeadJ, newSecondJ), weight);
                
                diJetPt_Zinc2jet->Fill(jet1Plus2.Pt(), weight);
                diJetPt_2_Zinc2jet->Fill(jet1Plus2.Pt(), weight);
                
                dPhiLepJet2_Zinc2jet->Fill(deltaPhi(lep1, newSecondJ), weight);
                dPhiLepJet2_2_Zinc2jet->Fill(deltaPhi(lep1, newSecondJ), weight);
                
                MeanNJetsHT_1D_Zinc2jet->Fill(jetsHT, weight*nGoodJets);
                MeanNJetsdRapidity_1D_Zinc2jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), weight*nGoodJets);
                MeanNJetsdRapidityFB_1D_Zinc2jet->Fill(ForwardJetRapidity - BackwardJetRapidity, weight*nGoodJets);
                
                MeanNJetsHT_Zinc2jet->Fill(jetsHT, nGoodJets, weight);
                MeanNJetsdRapidity_Zinc2jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), nGoodJets, weight);
                MeanNJetsdRapidityFB_Zinc2jet->Fill(ForwardJetRapidity - BackwardJetRapidity, nGoodJets, weight);
                SecondJetPtEta_Zinc2jet->Fill(jets[1].pt, fabs(jets[1].eta), weight);
                RatioJetPt21_Zinc2jet->Fill(jets[1].pt/jets[0].pt, weight);
                SecondJetPhi_Zinc2jet->Fill(jets[1].phi, weight);
                JetsHT_Zinc2jet->Fill(jetsHT, weight);
                JetsHT_20_30_Zinc2jet->Fill(jetsHT_20, weight);
                JetsHT_1_Zinc2jet->Fill(jetsHT, weight);
                JetsHT_2_Zinc2jet->Fill(jetsHT, weight);
                
                ptBal_Zinc2jet->Fill(jet1Plus2PlusZ.Pt(), weight);
                
                dPhiJets_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                dPhiJets_2_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                
                // BestdPhiJets_Zinc2jet->Fill(deltaPhi(bestTwoJets.first, bestTwoJets.second), weight);
                dEtaJets_Zinc2jet->Fill(leadJ.Eta() - secondJ.Eta(), weight);
                dEtaFirstJetZ_Zinc2jet->Fill(leadJ.Eta() - Z.Eta(), weight);
                dEtaSecondJetZ_Zinc2jet->Fill(secondJ.Eta() - Z.Eta(), weight);
                dEtaJet1Plus2Z_Zinc2jet->Fill(jet1Plus2.Eta() - Z.Eta(), weight);
                PHI_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                // BestPHI_Zinc2jet->Fill(PHI(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                PHI_T_Zinc2jet->Fill(PHI_T(lep1, lep2, leadJ, secondJ), weight);
                // BestPHI_T_Zinc2jet->Fill(PHI_T(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                SpT_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                // BestSpT_Zinc2jet->Fill(SpT(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                SpTJets_Zinc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                // BestSpTJets_Zinc2jet->Fill(SpTsub(bestTwoJets.first, bestTwoJets.second), weight);
                SPhi_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                // BestSPhi_Zinc2jet->Fill(SPhi(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                
                for ( int i =0 ; i < NbinsEta2D - 1 ; i++){
                    if ( fabs(jets[1].eta) >= j_Y_range[i] &&  fabs(jets[1].eta) < j_Y_range[i+1]) SecondJetPt_Zinc2jet_Eta[i]->Fill(fabs(jets[0].pt), weight);
                }
                
                for (unsigned short j(0); j < nGoodJets; j++){
                    AllJetPt_Zinc2jet->Fill(jets[j].pt, weight);
                    AllJetEta_Zinc2jet->Fill(jets[j].eta, weight);
                    AllJetPhi_Zinc2jet->Fill(jets[j].phi, weight);
                }

                if (Z.Pt() < 25){
                    ptBal_LowPt_Zinc2jet->Fill(jet1Plus2PlusZ.Pt(), weight);
                    dPhiJets_LowPt_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                    // BestdPhiJets_LowPt_Zinc2jet->Fill(deltaPhi(bestTwoJets.first, bestTwoJets.second), weight);
                    dPhiLeptons_LowPt_Zinc2jet->Fill(deltaPhi(lep1, lep2), weight);
                    PHI_LowPt_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                    // BestPHI_LowPt_Zinc2jet->Fill(PHI(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                    PHI_T_LowPt_Zinc2jet->Fill(PHI_T(lep1, lep2, leadJ, secondJ), weight);
                    // BestPHI_T_LowPt_Zinc2jet->Fill(PHI_T(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                    SpT_LowPt_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    // BestSpT_LowPt_Zinc2jet->Fill(SpT(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                    SpTJets_LowPt_Zinc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                    // BestSpTJets_LowPt_Zinc2jet->Fill(SpTsub(bestTwoJets.first, bestTwoJets.second), weight);
                    SpTLeptons_LowPt_Zinc2jet->Fill(SpTsub(lep1, lep2), weight);
                    SPhi_LowPt_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    // BestSPhi_LowPt_Zinc2jet->Fill(SPhi(lep1, lep2, bestTwoJets.first, bestTwoJets.second), weight);
                    if (SpT(lep1, lep2, leadJ, secondJ) < 0.5){
                        PHI_LowSpT_LowPt_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        SPhi_LowSpT_LowPt_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    }
                    else {
                        PHI_HighSpT_LowPt_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        SPhi_HighSpT_LowPt_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    }
                    if (SPhi(lep1, lep2, leadJ, secondJ) < 0.5){
                        SpT_LowSPhi_LowPt_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    }
                    else {
                        SpT_HighSPhi_LowPt_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    }
                }
                else {
                    ptBal_HighPt_Zinc2jet->Fill(jet1Plus2PlusZ.Pt(),weight);
                    dPhiJets_HighPt_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                    dPhiLeptons_HighPt_Zinc2jet->Fill(deltaPhi(lep1, lep2), weight);
                    PHI_HighPt_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                    PHI_T_HighPt_Zinc2jet->Fill(PHI_T(lep1, lep2, leadJ, secondJ), weight);
                    SpT_HighPt_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    SpTJets_HighPt_Zinc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                    SpTLeptons_HighPt_Zinc2jet->Fill(SpTsub(lep1, lep2), weight);
                    SPhi_HighPt_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    if (SpT(lep1, lep2, leadJ, secondJ) < 0.5){
                        PHI_LowSpT_HighPt_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        SPhi_LowSpT_HighPt_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    }
                    else {
                        PHI_HighSpT_HighPt_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        SPhi_HighSpT_HighPt_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    }
                    if (SPhi(lep1, lep2, leadJ, secondJ) < 0.5){
                        SpT_LowSPhi_HighPt_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    }
                    else {
                        SpT_HighSPhi_HighPt_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    }
                }
                if (SPhi(lep1, lep2, leadJ, secondJ) < 0.5){
                    SpT_LowSPhi_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                }
                else {
                    SpT_HighSPhi_Zinc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                }
                if (SpT(lep1, lep2, leadJ, secondJ) < 0.5){
                    PHI_LowSpT_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                    SPhi_LowSpT_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                }
                else {
                    PHI_HighSpT_Zinc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                    SPhi_HighSpT_Zinc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                }
                if (nGoodJets == 2){
                    nEventsExcl2Jets++;
                    ZNGoodJets_Zexc_NoWeight->Fill(2.);
                    ZMass_Zexc2jet->Fill(Z.M(), weight);
                    ZPt_Zexc2jet->Fill(Z.Pt(), weight);
                    ZRapidity_Zexc2jet->Fill(Z.Rapidity(), weight);
                    ZEta_Zexc2jet->Fill(Z.Eta(), weight);
                    lepPt_Zexc2jet->Fill(lepton1.pt, weight);
                    lepEta_Zexc2jet->Fill(lepton1.eta, weight);
                    if (doZ || doTT){
                        lepPt_Zexc2jet->Fill(lepton2.pt, weight);
                        lepEta_Zexc2jet->Fill(lepton2.eta, weight);
                    }
                    dPhiLeptons_Zexc2jet->Fill(deltaPhi(lep1, lep2), weight);
                    dEtaLeptons_Zexc2jet->Fill(lepton1.eta - lepton2.eta, weight);
                    SpTLeptons_Zexc2jet->Fill(SpTsub(lep1, lep2), weight);
                    SecondJetPt_Zexc2jet->Fill(jets[1].pt, weight);
                    SecondJetEta_Zexc2jet->Fill(jets[1].eta, weight);
                    SecondJetPhi_Zexc2jet->Fill(jets[1].phi, weight);
                    
                    //-- DPS Histograms
                    TwoJetsPtDiff_Zexc2jet->Fill(jet1Minus2.Pt(), weight);
                    JetsMass_Zexc2jet->Fill(jet1Plus2.M(), weight);
                    ptBal_Zexc2jet->Fill(jet1Plus2PlusZ.Pt(), weight);
                    dPhiJets_Zexc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                    dEtaJets_Zexc2jet->Fill(leadJ.Eta() - secondJ.Eta(), weight);
                    dEtaFirstJetZ_Zexc2jet->Fill(leadJ.Eta() - Z.Eta(), weight);
                    dEtaSecondJetZ_Zexc2jet->Fill(secondJ.Eta() - Z.Eta(), weight);
                    dEtaJet1Plus2Z_Zexc2jet->Fill(jet1Plus2.Eta() - Z.Eta(), weight);
                    PHI_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                    PHI_T_Zexc2jet->Fill(PHI_T(lep1, lep2, leadJ, secondJ), weight);
                    SpT_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    SpTJets_Zexc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                    SPhi_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    if (Z.Pt() < 25){
                        ptBal_LowPt_Zexc2jet->Fill(jet1Plus2PlusZ.Pt(), weight);
                        dPhiJets_LowPt_Zexc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                        dPhiLeptons_LowPt_Zexc2jet->Fill(deltaPhi(lep1, lep2), weight);
                        PHI_LowPt_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        PHI_T_LowPt_Zexc2jet->Fill(PHI_T(lep1, lep2, leadJ, secondJ), weight);
                        SpT_LowPt_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                        SpTJets_LowPt_Zexc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                        SpTLeptons_LowPt_Zexc2jet->Fill(SpTsub(lep1, lep2), weight);
                        SPhi_LowPt_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                        if (SpT(lep1, lep2, leadJ, secondJ) < 0.5){
                            PHI_LowSpT_LowPt_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                            SPhi_LowSpT_LowPt_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                        }
                        else {
                            PHI_HighSpT_LowPt_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                            SPhi_HighSpT_LowPt_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                        }
                        if (SPhi(lep1, lep2, leadJ, secondJ) < 0.5){
                            SpT_LowSPhi_LowPt_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                        }
                        else {
                            SpT_HighSPhi_LowPt_Zexc2jet ->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                        }
                    }
                    else {
                        ptBal_HighPt_Zexc2jet->Fill(jet1Plus2PlusZ.Pt(),weight);
                        dPhiJets_HighPt_Zexc2jet->Fill(deltaPhi(leadJ, secondJ), weight);
                        dPhiLeptons_HighPt_Zexc2jet->Fill(deltaPhi(lep1, lep2), weight);
                        PHI_HighPt_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        PHI_T_HighPt_Zexc2jet->Fill(PHI_T(lep1, lep2, leadJ, secondJ), weight);
                        SpT_HighPt_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                        SpTJets_HighPt_Zexc2jet->Fill(SpTsub(leadJ, secondJ), weight);
                        SpTLeptons_HighPt_Zexc2jet->Fill(SpTsub(lep1, lep2), weight);
                        SPhi_HighPt_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                        if (SpT(lep1, lep2, leadJ, secondJ) < 0.5){
                            PHI_LowSpT_HighPt_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                            SPhi_LowSpT_HighPt_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                        }
                        else {
                            PHI_HighSpT_HighPt_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                            SPhi_HighSpT_HighPt_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                        }
                        if (SPhi(lep1, lep2, leadJ, secondJ) < 0.5){
                            SpT_LowSPhi_HighPt_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                        }
                        else {
                            SpT_HighSPhi_HighPt_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                        }
                    }
                    if (SPhi(lep1, lep2, leadJ, secondJ) < 0.5){
                        SpT_LowSPhi_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    }
                    else {
                        SpT_HighSPhi_Zexc2jet->Fill(SpT(lep1, lep2, leadJ, secondJ), weight);
                    }
                    if (SpT(lep1, lep2, leadJ, secondJ) < 0.5){
                        PHI_LowSpT_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        SPhi_LowSpT_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    }
                    else {
                        PHI_HighSpT_Zexc2jet->Fill(PHI(lep1, lep2, leadJ, secondJ), weight);
                        SPhi_HighSpT_Zexc2jet->Fill(SPhi(lep1, lep2, leadJ, secondJ), weight);
                    }
                }
            }
            if (nGoodJetsAK8 >= 2){
                LeadingJetAK8Pt_Zinc2jet->Fill(jetsAK8[0].pt, weight);
                LeadingJetAK8Pt_2_Zinc2jet->Fill(jetsAK8[0].pt, weight);
                LepPtPlusLeadingJetAK8Pt_Zinc2jet->Fill(lepton1.pt + jetsAK8[0].pt, weight);
                LepPtPlusLeadingJetAK8Pt_2_Zinc2jet->Fill(lepton1.pt + jetsAK8[0].pt, weight);
                if (nGoodJetsAK8 == 2){
                    LeadingJetAK8Pt_Zexc2jet->Fill(jetsAK8[0].pt, weight);
                    LeadingJetAK8Pt_2_Zexc2jet->Fill(jetsAK8[0].pt, weight);
                    LepPtPlusLeadingJetAK8Pt_Zexc2jet->Fill(lepton1.pt + jetsAK8[0].pt, weight);
                    LepPtPlusLeadingJetAK8Pt_2_Zexc2jet->Fill(lepton1.pt + jetsAK8[0].pt, weight);
                }
            }
            if (nGoodJets >= 3) {
                ZNGoodJets_Zinc->Fill(3., weight);
                ZNGoodJetsFull_Zinc->Fill(3., weight);
                ZNGoodJets_Zinc_NoWeight->Fill(3.);
                ZMass_Zinc3jet->Fill(Z.M(), weight);
                MET_Zinc3jet->Fill(METpt, weight);
                METphi_Zinc3jet->Fill(METphi, weight);
                MT_Zinc3jet->Fill(MT, weight);
                //                ZPt_Zinc3jet->Fill(Z.Pt(), weight);
                ZRapidity_Zinc3jet->Fill(Z.Rapidity(), weight);
                ZEta_Zinc3jet->Fill(Z.Eta(), weight);
                lepPt_Zinc3jet->Fill(lepton1.pt, weight);
                lepEta_Zinc3jet->Fill(lepton1.eta, weight);
                if ( doZ || doTT){
                    lepPt_Zinc3jet->Fill(lepton2.pt, weight);
                    lepEta_Zinc3jet->Fill(lepton2.eta, weight);
                }
                dPhiLeptons_Zinc3jet->Fill(deltaPhi(lep1, lep2), weight);
                dEtaLeptons_Zinc3jet->Fill(lepton1.eta - lepton2.eta, weight);
                dRLeptons_Zinc3jet->Fill(deltaR(lepton1.phi, lepton1.eta, lepton2.phi, lepton2.eta), weight);
                SpTLeptons_Zinc3jet->Fill(SpTsub(lep1, lep2), weight);
                
                FirstHighestJetPt_Zinc3jet->Fill(jets[0].pt, weight);
                SecondHighestJetPt_Zinc3jet->Fill(jets[1].pt, weight);
                ThirdHighestJetPt_Zinc3jet->Fill(jets[2].pt, weight);
                ThirdJetEta_Zinc3jet->Fill(fabs(jets[2].eta), weight);
                ThirdJetEta_2_Zinc3jet->Fill(fabs(jets[2].eta), weight);
                ThirdJetEtaFull_Zinc3jet->Fill(jets[2].eta, weight);

                ThirdJetAbsRapidity_Zinc3jet->Fill(fabs(newThirdJ.Rapidity()), weight);
                ThirdJetAbsRapidity_2_Zinc3jet->Fill(fabs(newThirdJ.Rapidity()), weight);
                ThirdJetRapidityFull_Zinc3jet->Fill(newThirdJ.Rapidity(), weight);
                ThirdJetmass_Zinc3jet->Fill(newThirdJ.M(), weight);
                ThirdJetmass_1_Zinc3jet->Fill(newThirdJ.M(), weight);
                dRapidityJets_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), weight);
                dRapidityJets_2_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), weight);
                
                dRapidityJets_First_Third_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newThirdJ.Rapidity()), weight);
                dRapidityJets_2_First_Third_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newThirdJ.Rapidity()), weight);
                
                dRapidityJets_Second_Third_Zinc3jet->Fill(fabs(newSecondJ.Rapidity() - newThirdJ.Rapidity()), weight);
                dRapidityJets_2_Second_Third_Zinc3jet->Fill(fabs(newSecondJ.Rapidity() - newThirdJ.Rapidity()), weight);
                
                dRapidityJetsFB_Zinc3jet->Fill(ForwardJetRapidity - BackwardJetRapidity, weight);
                dRapidityJetsFB_2_Zinc3jet->Fill(ForwardJetRapidity - BackwardJetRapidity, weight);
                
                diJetMass_Zinc3jet->Fill(jet1Plus2.M(), weight);
                diJetMass_2_Zinc3jet->Fill(jet1Plus2.M(), weight);
                
                diJetPt_Zinc3jet->Fill(jet1Plus2.Pt(), weight);
                diJetPt_2_Zinc3jet->Fill(jet1Plus2.Pt(), weight);

                dPhiLepJet3_Zinc3jet->Fill(deltaPhi(lep1, newThirdJ), weight);
                dPhiLepJet3_2_Zinc3jet->Fill(deltaPhi(lep1, newThirdJ), weight);

                ThirdJetPtEta_Zinc3jet->Fill(jets[2].pt, fabs(jets[2].eta), weight);
                RatioJetPt32_Zinc3jet->Fill(jets[2].pt/jets[1].pt, weight);
                ThirdJetPhi_Zinc3jet->Fill(jets[2].phi, weight);
                JetsHT_Zinc3jet->Fill(jetsHT, weight);
                JetsHT_20_30_Zinc3jet->Fill(jetsHT_20, weight);
                JetsHT_1_Zinc3jet->Fill(jetsHT, weight);
                JetsHT_2_Zinc3jet->Fill(jetsHT, weight);
                
                for (unsigned short j(0); j < nGoodJets; j++){
                    AllJetPt_Zinc3jet->Fill(jets[j].pt, weight);
                    AllJetEta_Zinc3jet->Fill(jets[j].eta, weight);
                    AllJetPhi_Zinc3jet->Fill(jets[j].phi, weight);
                }
                if (nGoodJets == 3){
                    nEventsExcl3Jets++;
                    ZNGoodJets_Zexc_NoWeight->Fill(3.);
                    ZMass_Zexc3jet->Fill(Z.M(), weight);
                    ZPt_Zexc3jet->Fill(Z.Pt(), weight);
                    ZRapidity_Zexc3jet->Fill(Z.Rapidity(), weight);
                    ZEta_Zexc3jet->Fill(Z.Eta(), weight);
                    lepPt_Zexc3jet->Fill(lepton1.pt, weight);
                    lepEta_Zexc3jet->Fill(lepton1.eta, weight);
                    if (doZ || doTT){
                        lepPt_Zexc3jet->Fill(lepton2.pt, weight);
                        lepEta_Zexc3jet->Fill(lepton2.eta, weight);
                    }
                    dPhiLeptons_Zexc3jet->Fill(deltaPhi(lep1, lep2), weight);
                    dEtaLeptons_Zexc3jet->Fill(lepton1.eta - lepton2.eta, weight);
                    SpTLeptons_Zexc3jet->Fill(SpTsub(lep1, lep2), weight);

                }
            }
            if (nGoodJetsAK8 >= 3){
                LeadingJetAK8Pt_Zinc3jet->Fill(jetsAK8[0].pt, weight);
                LeadingJetAK8Pt_2_Zinc3jet->Fill(jetsAK8[0].pt, weight);
                LepPtPlusLeadingJetAK8Pt_Zinc3jet->Fill(lepton1.pt + jetsAK8[0].pt, weight);
                LepPtPlusLeadingJetAK8Pt_2_Zinc3jet->Fill(lepton1.pt + jetsAK8[0].pt, weight);
                if (nGoodJetsAK8 == 3){
                    LeadingJetAK8Pt_Zexc3jet->Fill(jetsAK8[0].pt, weight);
                    LeadingJetAK8Pt_2_Zexc3jet->Fill(jetsAK8[0].pt, weight);
                    LepPtPlusLeadingJetAK8Pt_Zexc3jet->Fill(lepton1.pt + jetsAK8[0].pt, weight);
                    LepPtPlusLeadingJetAK8Pt_2_Zexc3jet->Fill(lepton1.pt + jetsAK8[0].pt, weight);
                }
            }
            if (nGoodJets >= 4){
                ZNGoodJets_Zinc->Fill(4., weight);
                ZNGoodJetsFull_Zinc->Fill(4., weight);
                ZNGoodJets_Zinc_NoWeight->Fill(4.);
                ZMass_Zinc4jet->Fill(Z.M(), weight);
                //                ZPt_Zinc4jet->Fill(Z.Pt(), weight);
                ZRapidity_Zinc4jet->Fill(Z.Rapidity(), weight);
                ZEta_Zinc4jet->Fill(Z.Eta(), weight);
                lepPt_Zinc4jet->Fill(lepton1.pt, weight);
                lepEta_Zinc4jet->Fill(lepton1.eta, weight);
                if (doZ || doTT){
                    lepEta_Zinc4jet->Fill(lepton2.eta, weight);
                    lepPt_Zinc4jet->Fill(lepton2.pt, weight);
                }
                dPhiLeptons_Zinc4jet->Fill(deltaPhi(lep1, lep2), weight);
                dEtaLeptons_Zinc4jet->Fill(lepton1.eta - lepton2.eta, weight);
                dRLeptons_Zinc4jet->Fill(deltaR(lepton1.phi, lepton1.eta, lepton2.phi, lepton2.eta), weight);
                SpTLeptons_Zinc4jet->Fill(SpTsub(lep1, lep2), weight);
                
                FirstHighestJetPt_Zinc4jet->Fill(jets[0].pt, weight);
                SecondHighestJetPt_Zinc4jet->Fill(jets[1].pt, weight);
                ThirdHighestJetPt_Zinc4jet->Fill(jets[2].pt, weight);
                FourthJetEta_Zinc4jet->Fill(fabs(jets[3].eta), weight);
                FourthJetEta_2_Zinc4jet->Fill(fabs(jets[3].eta), weight);
                FourthJetEtaFull_Zinc4jet->Fill(jets[3].eta, weight);
                FourthJetAbsRapidity_Zinc4jet->Fill(fabs(newFourthJ.Rapidity()), weight);
                FourthJetAbsRapidity_2_Zinc4jet->Fill(fabs(newFourthJ.Rapidity()), weight);
                FourthJetRapidityFull_Zinc4jet->Fill(newFourthJ.Rapidity(), weight);
                FourthJetmass_Zinc4jet->Fill(newFourthJ.M(), weight);
                FourthJetmass_1_Zinc4jet->Fill(newFourthJ.M(), weight);
                dRapidityJets_Zinc4jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), weight);
                dRapidityJets_2_Zinc4jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), weight);
                
                dRapidityJetsFB_Zinc4jet->Fill(ForwardJetRapidity - BackwardJetRapidity, weight);
                dRapidityJetsFB_2_Zinc4jet->Fill(ForwardJetRapidity - BackwardJetRapidity, weight);
                
                diJetMass_Zinc4jet->Fill(jet1Plus2.M(), weight);
                diJetMass_2_Zinc4jet->Fill(jet1Plus2.M(), weight);
                
                diJetPt_Zinc4jet->Fill(jet1Plus2.Pt(), weight);
                diJetPt_2_Zinc4jet->Fill(jet1Plus2.Pt(), weight);

                dPhiLepJet4_Zinc4jet->Fill(deltaPhi(lep1, newFourthJ), weight);
                dPhiLepJet4_2_Zinc4jet->Fill(deltaPhi(lep1, newFourthJ), weight);

                FourthJetPtEta_Zinc4jet->Fill(jets[3].pt, fabs(jets[3].eta), weight);
                FourthJetPhi_Zinc4jet->Fill(jets[3].phi, weight);
                JetsHT_Zinc4jet->Fill(jetsHT, weight);
                JetsHT_1_Zinc4jet->Fill(jetsHT, weight);
                JetsHT_2_Zinc4jet->Fill(jetsHT, weight);
                for (unsigned short j(0); j < nGoodJets; j++){
                    AllJetPt_Zinc4jet->Fill(jets[j].pt, weight);
                    AllJetEta_Zinc4jet->Fill(jets[j].eta, weight);
                    AllJetPhi_Zinc4jet->Fill(jets[j].phi, weight);
                }


                if (nGoodJets == 4){
                    ZNGoodJets_Zexc_NoWeight->Fill(4.);
                    ZMass_Zexc4jet->Fill(Z.M(), weight);
                    ZPt_Zexc4jet->Fill(Z.Pt(), weight);
                    ZRapidity_Zexc4jet->Fill(Z.Rapidity(), weight);
                    ZEta_Zexc4jet->Fill(Z.Eta(), weight);
                    lepPt_Zexc4jet->Fill(lepton1.pt, weight);
                    lepEta_Zexc4jet->Fill(lepton1.eta, weight);
                    if (doZ || doTT){
                        lepPt_Zexc4jet->Fill(lepton2.pt, weight);
                        lepEta_Zexc4jet->Fill(lepton2.eta, weight);
                    }
                    dPhiLeptons_Zexc4jet->Fill(deltaPhi(lep1, lep2), weight);
                    dEtaLeptons_Zexc4jet->Fill(lepton1.eta - lepton2.eta, weight);
                    SpTLeptons_Zexc4jet->Fill(SpTsub(lep1, lep2), weight);

                }
            }
            if (nGoodJets >= 5){
                ZNGoodJets_Zinc->Fill(5., weight);
                ZNGoodJetsFull_Zinc->Fill(5., weight);
                ZNGoodJets_Zinc_NoWeight->Fill(5.);
                ZMass_Zinc5jet->Fill(Z.M(), weight);
                ZPt_Zinc5jet->Fill(Z.Pt(), weight);
                ZRapidity_Zinc5jet->Fill(Z.Rapidity(), weight);
                ZEta_Zinc5jet->Fill(Z.Eta(), weight);
                lepPt_Zinc5jet->Fill(lepton1.pt, weight);
                lepEta_Zinc5jet->Fill(lepton1.eta, weight);
                if (doZ || doTT){
                    lepPt_Zinc5jet->Fill(lepton2.pt, weight);
                    lepEta_Zinc5jet->Fill(lepton2.eta, weight);
                }
                dPhiLeptons_Zinc5jet->Fill(deltaPhi(lep1, lep2), weight);
                dEtaLeptons_Zinc5jet->Fill(lepton1.eta - lepton2.eta, weight);
                dRLeptons_Zinc5jet->Fill(deltaR(lepton1.phi, lepton1.eta, lepton2.phi, lepton2.eta), weight);
                SpTLeptons_Zinc5jet->Fill(SpTsub(lep1, lep2), weight);
                
                FirstHighestJetPt_Zinc5jet->Fill(jets[0].pt, weight);
                SecondHighestJetPt_Zinc5jet->Fill(jets[1].pt, weight);
                ThirdHighestJetPt_Zinc5jet->Fill(jets[2].pt, weight);
                FifthJetEta_Zinc5jet->Fill(fabs(jets[4].eta), weight);
                FifthJetEta_2_Zinc5jet->Fill(fabs(jets[4].eta), weight);
                FifthJetEtaFull_Zinc5jet->Fill(jets[4].eta, weight);
                FifthJetPtEta_Zinc5jet->Fill(jets[4].pt, fabs(jets[4].eta), weight);
                FifthJetPhi_Zinc5jet->Fill(jets[4].phi, weight);
                JetsHT_Zinc5jet->Fill(jetsHT, weight);
                JetsHT_1_Zinc5jet->Fill(jetsHT, weight);
                JetsHT_2_Zinc5jet->Fill(jetsHT, weight);

				//dPhiLepJet5_Zinc5jet->Fill(deltaPhi(lep1, newFifthJ), weight);
				//dPhiLepJet5_2_Zinc5jet->Fill(deltaPhi(lep1, newFifthJ), weight);

                if (nGoodJets == 5){
                    ZNGoodJets_Zexc_NoWeight->Fill(5.);
                    ZMass_Zexc5jet->Fill(Z.M(), weight);
                    ZPt_Zexc5jet->Fill(Z.Pt(), weight);
                    ZRapidity_Zexc5jet->Fill(Z.Rapidity(), weight);
                    ZEta_Zexc5jet->Fill(Z.Eta(), weight);
                    lepPt_Zexc5jet->Fill(lepton1.pt, weight);
                    lepEta_Zexc5jet->Fill(lepton1.eta, weight);
                    if (doZ || doTT){
                        lepPt_Zexc5jet->Fill(lepton2.pt, weight);
                        lepEta_Zexc5jet->Fill(lepton2.eta, weight);
                    }
                    dPhiLeptons_Zexc5jet->Fill(deltaPhi(lep1, lep2), weight);
                    dEtaLeptons_Zexc5jet->Fill(lepton1.eta - lepton2.eta, weight);
                    SpTLeptons_Zexc5jet->Fill(SpTsub(lep1, lep2), weight);
                }
            }
            if (nGoodJets >= 6){
                ZNGoodJets_Zinc->Fill(6., weight);
                ZNGoodJetsFull_Zinc->Fill(6., weight);
                ZNGoodJets_Zinc_NoWeight->Fill(6.);
                ZMass_Zinc6jet->Fill(Z.M(), weight);
                ZPt_Zinc6jet->Fill(Z.Pt(), weight);
                ZRapidity_Zinc6jet->Fill(Z.Rapidity(), weight);
                ZEta_Zinc6jet->Fill(Z.Eta(), weight);
                
                FirstHighestJetPt_Zinc6jet->Fill(jets[0].pt, weight);
                SecondHighestJetPt_Zinc6jet->Fill(jets[1].pt, weight);
                ThirdHighestJetPt_Zinc6jet->Fill(jets[2].pt, weight);
                SixthJetEta_Zinc6jet->Fill(fabs(jets[5].eta), weight);
                SixthJetEtaFull_Zinc6jet->Fill(jets[5].eta, weight);
                SixthJetPtEta_Zinc6jet->Fill(jets[5].pt, fabs(jets[5].eta), weight);
                SixthJetPhi_Zinc6jet->Fill(jets[5].phi, weight);
                JetsHT_Zinc6jet->Fill(jetsHT, weight);
                JetsHT_1_Zinc6jet->Fill(jetsHT, weight);
                if (nGoodJets == 6){
                    ZNGoodJets_Zexc_NoWeight->Fill(6.);
                    ZMass_Zexc6jet->Fill(Z.M(), weight);
                    ZPt_Zexc6jet->Fill(Z.Pt(), weight);
                    ZRapidity_Zexc6jet->Fill(Z.Rapidity(), weight);
                    ZEta_Zexc6jet->Fill(Z.Eta(), weight);
                }
            }
            if (nGoodJets >= 7){
                ZNGoodJets_Zinc->Fill(7., weight);
                ZNGoodJetsFull_Zinc->Fill(7., weight);
            }
            //more bins in Single Lepton dataset than Double --> xsec bigger
            if (nGoodJets >= 8) ZNGoodJets_Zinc->Fill(8., weight);
            if (nGoodJets >= 9)  ZNGoodJets_Zinc->Fill(9., weight);
            if (nGoodJets >= 10) ZNGoodJets_Zinc->Fill(10., weight);
        }
        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        
        //=======================================================================================================//
        //             Unfolding              //
        //====================================//
        if (hasRecoInfo && hasGenInfo && PRINTEVENTINFO && jentry == eventOfInterest) {
            cout << __LINE__ << " PRINTEVENTINFO: Requirements for filling hresponse histos --- " << endl;
            cout << __LINE__ << " PRINTEVENTINFO: hasRecoInfo = " << hasRecoInfo << endl;
            cout << __LINE__ << " PRINTEVENTINFO: hasGenInfo = " << hasGenInfo << endl;
            cout << __LINE__ << " PRINTEVENTINFO: passesGenLeptonCut = " << passesGenLeptonCut << endl;
            cout << __LINE__ << " PRINTEVENTINFO: passesLeptonCut = " << passesLeptonCut << endl;
        }
        if (hasRecoInfo && hasGenInfo && passesGenLeptonCut && passesLeptonCut){
            //-- jet multiplicity exc
            hresponseZNGoodJets_Zexc->Fill(nGoodJets, nGoodGenJets, weight);
            hresponseZNGoodJetsFull_Zexc->Fill(nGoodJets, nGoodGenJets, weight);
            hresponseZNGoodJets_Zinc->Fill(0., 0., weight);
            hresponseZNGoodJetsFull_Zinc->Fill(0., 0., weight);
            
            if (nGoodGenJets_20 >= 1 && nGoodJets_20 >= 1){
                hresponseFirstJetPt_Zinc1jet  ->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                hresponseFirstJetPt_1_Zinc1jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                hresponseFirstJetPt_2_Zinc1jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                hresponseJetsHT_20_Zinc1jet->Fill(jetsHT_20, genJetsHT_20, weight);

                //andrew - alpha-s
                hresponseLeadingJetPt_Zinc1jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                // hresponseLeadingJetPt_1_Zinc1jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                hresponseLeadingJetPt_2_Zinc1jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                hresponseLepPtPlusLeadingJetPt_Zinc1jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                // hresponseLepPtPlusLeadingJetPt_1_Zinc1jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                hresponseLepPtPlusLeadingJetPt_2_Zinc1jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                hresponseZPt_Zinc1jet->Fill(Z.Pt(), genZ.Pt(), weight);
                // hresponseZPt_1_Zinc1jet->Fill(Z.Pt(), genZ.Pt(), weight);
                hresponseZPt_2_Zinc1jet->Fill(Z.Pt(), genZ.Pt(), weight);
                hresponseZPtPlusLeadingJetPt_Zinc1jet->Fill(Z.Pt()+jets_20[0].pt, genZ.Pt()+genJets_20[0].pt, weight);
                // hresponseZPtPlusLeadingJetPt_1_Zinc1jet->Fill(Z.Pt()+jets_20[0].pt, genZ.Pt()+genJets_20[0].pt, weight);
                hresponseZPtPlusLeadingJetPt_2_Zinc1jet->Fill(Z.Pt()+jets_20[0].pt, genZ.Pt()+genJets_20[0].pt, weight);

                // looking to study particle-level corrections of ratios
                LeadingJetPt_MIGRATIONS_Zinc1jet->Fill(jets_20[0].pt, weight);
                genLeadingJetPt_MIGRATIONS_Zinc1jet->Fill(genJets_20[0].pt, genWeight);
                LepPtPlusLeadingJetPt_MIGRATIONS_Zinc1jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                genLepPtPlusLeadingJetPt_MIGRATIONS_Zinc1jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);

                if (nGoodGenJets_20 == 1 && nGoodJets_20 == 1){
                    hresponseLeadingJetPt_Zexc1jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                    hresponseLeadingJetPt_2_Zexc1jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                    hresponseLepPtPlusLeadingJetPt_Zexc1jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                    hresponseLepPtPlusLeadingJetPt_2_Zexc1jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                }
            }
            if (nGoodGenJets_20 >= 2 && nGoodJets_20 >= 2){
                hresponseSecondJetPt_Zinc2jet  ->Fill(jets_20[1].pt, genJets_20[1].pt, weight);
                hresponseSecondJetPt_1_Zinc2jet->Fill(jets_20[1].pt, genJets_20[1].pt, weight);
                hresponseSecondJetPt_2_Zinc2jet->Fill(jets_20[1].pt, genJets_20[1].pt, weight);
                hresponseJetsHT_20_Zinc2jet->Fill(jetsHT_20, genJetsHT_20, weight);

                //andrew - alpha-s
                hresponseLeadingJetPt_Zinc2jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                // hresponseLeadingJetPt_1_Zinc2jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                hresponseLeadingJetPt_2_Zinc2jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                hresponseHTover2_Zinc2jet->Fill((jets_20[0].pt+jets_20[1].pt)/2., (genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                // hresponseHTover2_1_Zinc2jet->Fill((jets_20[0].pt+jets_20[1].pt)/2., (genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                hresponseHTover2_2_Zinc2jet->Fill((jets_20[0].pt+jets_20[1].pt)/2., (genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                hresponseLepPtPlusLeadingJetPt_Zinc2jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                // hresponseLepPtPlusLeadingJetPt_1_Zinc2jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                hresponseLepPtPlusLeadingJetPt_2_Zinc2jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                hresponseLepPtPlusHTover2_Zinc2jet->Fill(lepton1.pt + (jets_20[0].pt+jets_20[1].pt)/2., genLep1.Pt()+(genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                // hresponseLepPtPlusHTover2_1_Zinc2jet->Fill(lepton1.pt + (jets_20[0].pt+jets_20[1].pt)/2., genLep1.Pt()+(genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                hresponseLepPtPlusHTover2_2_Zinc2jet->Fill(lepton1.pt + (jets_20[0].pt+jets_20[1].pt)/2., genLep1.Pt()+(genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                hresponseZPt_Zinc2jet->Fill(Z.Pt(), genZ.Pt(), weight);
                // hresponseZPt_1_Zinc2jet->Fill(Z.Pt(), genZ.Pt(), weight);
                hresponseZPt_2_Zinc2jet->Fill(Z.Pt(), genZ.Pt(), weight);
                hresponseZPtPlusLeadingJetPt_Zinc2jet->Fill(Z.Pt()+jets_20[0].pt, genZ.Pt()+genJets_20[0].pt, weight);
                // hresponseZPtPlusLeadingJetPt_1_Zinc2jet->Fill(Z.Pt()+jets_20[0].pt, genZ.Pt()+genJets_20[0].pt, weight);
                hresponseZPtPlusLeadingJetPt_2_Zinc2jet->Fill(Z.Pt()+jets_20[0].pt, genZ.Pt()+genJets_20[0].pt, weight);
                hresponseZPtPlusHTover2_Zinc2jet->Fill(Z.Pt()+(jets_20[0].pt+jets_20[1].pt)/2., genZ.Pt()+(genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                // hresponseZPtPlusHTover2_1_Zinc2jet->Fill(Z.Pt()+(jets_20[0].pt+jets_20[1].pt)/2., genZ.Pt()+(genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                hresponseZPtPlusHTover2_2_Zinc2jet->Fill(Z.Pt()+(jets_20[0].pt+jets_20[1].pt)/2., genZ.Pt()+(genJets_20[0].pt+genJets_20[1].pt)/2., weight);

                // looking to study particle-level corrections of ratios
                LeadingJetPt_MIGRATIONS_Zinc2jet->Fill(jets_20[0].pt, weight);
                genLeadingJetPt_MIGRATIONS_Zinc2jet->Fill(genJets_20[0].pt, genWeight);
                LepPtPlusLeadingJetPt_MIGRATIONS_Zinc2jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                genLepPtPlusLeadingJetPt_MIGRATIONS_Zinc2jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);

                if (nGoodGenJets_20 == 2 && nGoodJets_20 == 2){
                    hresponseLeadingJetPt_Zexc2jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                    hresponseLeadingJetPt_2_Zexc2jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                    hresponseLepPtPlusLeadingJetPt_Zexc2jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                    hresponseLepPtPlusLeadingJetPt_2_Zexc2jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                }
            }
            if (nGoodGenJets_20 >= 3 && nGoodJets_20 >= 3){
                hresponseThirdJetPt_Zinc3jet  ->Fill(jets_20[2].pt, genJets_20[2].pt, weight);
                hresponseThirdJetPt_1_Zinc3jet->Fill(jets_20[2].pt, genJets_20[2].pt, weight);
                hresponseThirdJetPt_2_Zinc3jet->Fill(jets_20[2].pt, genJets_20[2].pt, weight);
                hresponseJetsHT_20_Zinc3jet->Fill(jetsHT_20, genJetsHT_20, weight);

                //andrew - alpha-s
                hresponseLeadingJetPt_Zinc3jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                // hresponseLeadingJetPt_1_Zinc3jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                hresponseLeadingJetPt_2_Zinc3jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                hresponseHTover2_Zinc3jet->Fill((jets_20[0].pt+jets_20[1].pt)/2., (genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                // hresponseHTover2_1_Zinc3jet->Fill((jets_20[0].pt+jets_20[1].pt)/2., (genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                hresponseHTover2_2_Zinc3jet->Fill((jets_20[0].pt+jets_20[1].pt)/2., (genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                hresponseLepPtPlusLeadingJetPt_Zinc3jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                // hresponseLepPtPlusLeadingJetPt_1_Zinc3jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                hresponseLepPtPlusLeadingJetPt_2_Zinc3jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                hresponseLepPtPlusHTover2_Zinc3jet->Fill(lepton1.pt + (jets_20[0].pt+jets_20[1].pt)/2., genLep1.Pt()+(genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                // hresponseLepPtPlusHTover2_1_Zinc3jet->Fill(lepton1.pt + (jets_20[0].pt+jets_20[1].pt)/2., genLep1.Pt()+(genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                hresponseLepPtPlusHTover2_2_Zinc3jet->Fill(lepton1.pt + (jets_20[0].pt+jets_20[1].pt)/2., genLep1.Pt()+(genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                hresponseZPt_Zinc3jet->Fill(Z.Pt(), genZ.Pt(), weight);
                // hresponseZPt_1_Zinc3jet->Fill(Z.Pt(), genZ.Pt(), weight);
                hresponseZPt_2_Zinc3jet->Fill(Z.Pt(), genZ.Pt(), weight);
                hresponseZPtPlusLeadingJetPt_Zinc3jet->Fill(Z.Pt()+jets_20[0].pt, genZ.Pt()+genJets_20[0].pt, weight);
                // hresponseZPtPlusLeadingJetPt_1_Zinc3jet->Fill(Z.Pt()+jets_20[0].pt, genZ.Pt()+genJets_20[0].pt, weight);
                hresponseZPtPlusLeadingJetPt_2_Zinc3jet->Fill(Z.Pt()+jets_20[0].pt, genZ.Pt()+genJets_20[0].pt, weight);
                hresponseZPtPlusHTover2_Zinc3jet->Fill(Z.Pt()+(jets_20[0].pt+jets_20[1].pt)/2., genZ.Pt()+(genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                // hresponseZPtPlusHTover2_1_Zinc3jet->Fill(Z.Pt()+(jets_20[0].pt+jets_20[1].pt)/2., genZ.Pt()+(genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                hresponseZPtPlusHTover2_2_Zinc3jet->Fill(Z.Pt()+(jets_20[0].pt+jets_20[1].pt)/2., genZ.Pt()+(genJets_20[0].pt+genJets_20[1].pt)/2., weight);

                // looking to study particle-level corrections of ratios
                LeadingJetPt_MIGRATIONS_Zinc3jet->Fill(jets_20[0].pt, weight);
                genLeadingJetPt_MIGRATIONS_Zinc3jet->Fill(genJets_20[0].pt, genWeight);
                LepPtPlusLeadingJetPt_MIGRATIONS_Zinc3jet->Fill(lepton1.pt + jets_20[0].pt, weight);
                genLepPtPlusLeadingJetPt_MIGRATIONS_Zinc3jet->Fill(genLep1.Pt()+genJets_20[0].pt, genWeight);

                if (nGoodGenJets_20 == 3 && nGoodJets_20 == 3){
                    hresponseLeadingJetPt_Zexc3jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                    hresponseLeadingJetPt_2_Zexc3jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                    hresponseLepPtPlusLeadingJetPt_Zexc3jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                    hresponseLepPtPlusLeadingJetPt_2_Zexc3jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                }
            }
            if (nGoodGenJets_20 >= 4 && nGoodJets_20 >= 4){
                hresponseFourthJetPt_Zinc4jet  ->Fill(jets_20[3].pt, genJets_20[3].pt, weight);
                hresponseFourthJetPt_1_Zinc4jet->Fill(jets_20[3].pt, genJets_20[3].pt, weight);
                hresponseFourthJetPt_2_Zinc4jet->Fill(jets_20[3].pt, genJets_20[3].pt, weight);

                //andrew - alpha-s
                hresponseLeadingJetPt_Zinc4jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                // hresponseLeadingJetPt_1_Zinc4jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                hresponseLeadingJetPt_2_Zinc4jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                hresponseHTover2_Zinc4jet->Fill((jets_20[0].pt+jets_20[1].pt)/2., (genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                // hresponseHTover2_1_Zinc4jet->Fill((jets_20[0].pt+jets_20[1].pt)/2., (genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                hresponseHTover2_2_Zinc4jet->Fill((jets_20[0].pt+jets_20[1].pt)/2., (genJets_20[0].pt+genJets_20[1].pt)/2., weight);
                hresponseLepPtPlusLeadingJetPt_Zinc4jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                // hresponseLepPtPlusLeadingJetPt_1_Zinc4jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                hresponseLepPtPlusLeadingJetPt_2_Zinc4jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                hresponseLepPtPlusHTover2_Zinc4jet->Fill(lepton1.pt + (jets_20[0].pt+jets_20[1].pt)/2., genLep1.Pt()+(genJets_20[0].pt + genJets_20[1].pt)/2., weight);
                // hresponseLepPtPlusHTover2_1_Zinc4jet->Fill(lepton1.pt + (jets_20[0].pt+jets_20[1].pt)/2., genLep1.Pt()+(genJets_20[0].pt + genJets_20[1].pt)/2., weight);
                hresponseLepPtPlusHTover2_2_Zinc4jet->Fill(lepton1.pt + (jets_20[0].pt+jets_20[1].pt)/2., genLep1.Pt()+(genJets_20[0].pt + genJets_20[1].pt)/2., weight);
                hresponseZPt_Zinc4jet->Fill(Z.Pt(), genZ.Pt(), weight);
                // hresponseZPt_1_Zinc4jet->Fill(Z.Pt(), genZ.Pt(), weight);
                hresponseZPt_2_Zinc4jet->Fill(Z.Pt(), genZ.Pt(), weight);
                hresponseZPtPlusLeadingJetPt_Zinc4jet->Fill(Z.Pt()+jets_20[0].pt, genZ.Pt()+genJets_20[0].pt, weight);
                // hresponseZPtPlusLeadingJetPt_1_Zinc4jet->Fill(Z.Pt()+jets_20[0].pt, genZ.Pt()+genJets_20[0].pt, weight);
                hresponseZPtPlusLeadingJetPt_2_Zinc4jet->Fill(Z.Pt()+jets_20[0].pt, genZ.Pt()+genJets_20[0].pt, weight);
                hresponseZPtPlusHTover2_Zinc4jet->Fill(Z.Pt()+(jets_20[0].pt+jets_20[1].pt)/2., genZ.Pt()+(genJets_20[0].pt + genJets_20[1].pt)/2., weight);
                // hresponseZPtPlusHTover2_1_Zinc4jet->Fill(Z.Pt()+(jets_20[0].pt+jets_20[1].pt)/2., genZ.Pt()+(genJets_20[0].pt + genJets_20[1].pt)/2., weight);
                hresponseZPtPlusHTover2_2_Zinc4jet->Fill(Z.Pt()+(jets_20[0].pt+jets_20[1].pt)/2., genZ.Pt()+(genJets_20[0].pt + genJets_20[1].pt)/2., weight);

                // if (nGoodGenJets_20 == 4 && nGoodJets_20 == 4){
                //     hresponseLeadingJetPt_Zexc4jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                //     hresponseLeadingJetPt_2_Zexc4jet->Fill(jets_20[0].pt, genJets_20[0].pt, weight);
                //     hresponseLepPtPlusLeadingJetPt_Zexc4jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                //     hresponseLepPtPlusLeadingJetPt_2_Zexc4jet->Fill(lepton1.pt + jets_20[0].pt, genLep1.Pt()+genJets_20[0].pt, weight);
                // }
            }
            if (nGoodGenJets_20 >= 5 && nGoodJets_20 >= 5){
                hresponseFifthJetPt_Zinc5jet  ->Fill(jets_20[4].pt, genJets_20[4].pt, weight);
                hresponseFifthJetPt_1_Zinc5jet->Fill(jets_20[4].pt, genJets_20[4].pt, weight);
                hresponseFifthJetPt_2_Zinc5jet->Fill(jets_20[4].pt, genJets_20[4].pt, weight);
            }

    	    if (nGenJetsDR02Cut >= 1 && nJetsDR02Cut >= 1){
                hresponsedRLepCloseJet_Zinc1jet->Fill(mindRjmu, genMindRjmu, weight);
                hresponsedRLepCloseJet_2_Zinc1jet->Fill(mindRjmu, genMindRjmu, weight);
                if(genJetsDR02[0].pt > 500 && jetsDR02[0].pt > 500){
                    hresponsedRLepCloseJetCo_Zinc1jet->Fill(mindRjmu, genMindRjmu, weight);
                    hresponsedRLepCloseJetCo_2_Zinc1jet->Fill(mindRjmu, genMindRjmu, weight);
                }
	        }
	    
            if (nGenJetsPt100DR04 >= 1 && nJetsPt100DR04 >= 1){
                if(genJetsPt100DR04[0].pt > 300 && jetsPt100DR04[0].pt > 300){
                    hresponsedRptmin100LepCloseJetCo300dR04_Zinc1jet->Fill(mindRj04Pt100mu, genMindRj04Pt100mu, weight);
                    hresponsedRptmin100LepCloseJetCo300dR04_2_Zinc1jet->Fill(mindRj04Pt100mu, genMindRj04Pt100mu, weight);
                }
            }

            if (nGenJetsPt100DR04 >= 2 && nJetsPt100DR04 >= 2){
                if(genJetsPt100DR04[0].pt > 300 && jetsPt100DR04[0].pt > 300){
                    if (deltaPhi(genLeadJ, genSecondJ) > (PI - 0.3) && deltaPhi(leadJ, secondJ) > (PI - 0.3)){
                    hresponsedRptmin100LepCloseDiJetCo300dR04_Zinc2jet->Fill(mindRdijet04Pt100mu, genMindRdijet04Pt100mu, weight);
                    hresponsedRptmin100LepCloseDiJetCo300dR04_2_Zinc2jet->Fill(mindRdijet04Pt100mu, genMindRdijet04Pt100mu, weight);
                    }
                }
            }
            
            //-- First Jet Pt
            if (nGoodGenJets >= 1 && nGoodJets >= 1){
                hresponseZNGoodJets_Zinc->Fill(1., 1., weight);
                hresponseZNGoodJetsFull_Zinc->Fill(1., 1., weight);
                hresponseFirstJetEta_Zinc1jet->Fill(fabs(jets[0].eta), fabs(genJets[0].eta), weight);
                hresponseJetsHT_Zinc1jet->Fill(jetsHT, genJetsHT, weight);
                hresponseJetsHT_20_30_Zinc1jet->Fill(jetsHT_20, genJetsHT_20, weight);

                hresponseFirstJetEta_2_Zinc1jet->Fill(fabs(jets[0].eta), fabs(genJets[0].eta), weight);
                hresponseJetsHT_1_Zinc1jet->Fill(jetsHT, genJetsHT, weight);
                hresponseJetsHT_2_Zinc1jet->Fill(jetsHT, genJetsHT, weight);
                
				hresponsedPhiLepJet1_Zinc1jet->Fill(deltaPhi(lep1, newLeadJ), deltaPhi(genLep1, genNewLeadJ), weight);
				hresponsedPhiLepJet1_2_Zinc1jet->Fill(deltaPhi(lep1, newLeadJ), deltaPhi(genLep1, genNewLeadJ), weight);
                
                hresponseFirstJetAbsRapidity_Zinc1jet->Fill(fabs(newLeadJ.Rapidity()), fabs(genNewLeadJ.Rapidity()), weight);
                hresponseFirstJetAbsRapidity_2_Zinc1jet->Fill(fabs(newLeadJ.Rapidity()), fabs(genNewLeadJ.Rapidity()), weight);

                hresponseLepPt_Zinc1jet->Fill(lepton1.pt, genLep1.Pt(), weight);
                hresponseMT_Zinc1jet->Fill(MT, genMT, weight);

            }
            if (nGoodGenJetsAK8 >= 1 && nGoodJetsAK8 >= 1){
                hresponseLeadingJetAK8Pt_Zinc1jet->Fill(jetsAK8[0].pt, genJetsAK8[0].pt, weight);
                hresponseLeadingJetAK8Pt_2_Zinc1jet->Fill(jetsAK8[0].pt, genJetsAK8[0].pt, weight);
                hresponseLepPtPlusLeadingJetAK8Pt_Zinc1jet->Fill(lepton1.pt + jetsAK8[0].pt, genLep1.Pt()+genJetsAK8[0].pt, weight);
                hresponseLepPtPlusLeadingJetAK8Pt_2_Zinc1jet->Fill(lepton1.pt + jetsAK8[0].pt, genLep1.Pt()+genJetsAK8[0].pt, weight);
                if (nGoodGenJetsAK8 == 1 && nGoodJetsAK8 == 1){
                    hresponseLeadingJetAK8Pt_Zexc1jet->Fill(jetsAK8[0].pt, genJetsAK8[0].pt, weight);
                    hresponseLeadingJetAK8Pt_2_Zexc1jet->Fill(jetsAK8[0].pt, genJetsAK8[0].pt, weight);
                    hresponseLepPtPlusLeadingJetAK8Pt_Zexc1jet->Fill(lepton1.pt + jetsAK8[0].pt, genLep1.Pt()+genJetsAK8[0].pt, weight);
                    hresponseLepPtPlusLeadingJetAK8Pt_2_Zexc1jet->Fill(lepton1.pt + jetsAK8[0].pt, genLep1.Pt()+genJetsAK8[0].pt, weight);
                }
            }
            //-- Second Jet Pt
            if (nGoodGenJets >= 2 && nGoodJets >= 2){
                hresponseZNGoodJets_Zinc->Fill(2., 2., weight);
                hresponseZNGoodJetsFull_Zinc->Fill(2., 2., weight);
                hresponseSecondJetEta_Zinc2jet->Fill(fabs(jets[1].eta), fabs(genJets[1].eta), weight);
                hresponseJetsHT_Zinc2jet->Fill(jetsHT, genJetsHT, weight);
                hresponseJetsHT_20_30_Zinc2jet->Fill(jetsHT_20, genJetsHT_20, weight);
                
                hresponseSecondJetEta_2_Zinc2jet->Fill(fabs(jets[1].eta), fabs(genJets[1].eta), weight);
                hresponseJetsHT_1_Zinc2jet->Fill(jetsHT, genJetsHT, weight);
                hresponseJetsHT_2_Zinc2jet->Fill(jetsHT, genJetsHT, weight);
                
                hresponsedRapidityJets_Zinc2jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), weight);
                hresponsedRapidityJets_2_Zinc2jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), weight);
                
                hresponsedRapidityJetsFB_Zinc2jet->Fill(ForwardJetRapidity - BackwardJetRapidity, genForwardJetRapidity - genBackwardJetRapidity, weight);
                hresponsedRapidityJetsFB_2_Zinc2jet->Fill(ForwardJetRapidity - BackwardJetRapidity, genForwardJetRapidity - genBackwardJetRapidity, weight);
                
                hresponsedPhiJets_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), deltaPhi(genLeadJ, genSecondJ), weight);
                hresponsedPhiJets_2_Zinc2jet->Fill(deltaPhi(leadJ, secondJ), deltaPhi(genLeadJ, genSecondJ), weight);
                
                hresponsedPhiJetsFB_Zinc2jet->Fill(deltaPhi(vJetYOrdered[0], vJetYOrdered[vJetYOrdered.size() - 1] ), deltaPhi(genvJetYOrdered[0], genvJetYOrdered[genvJetYOrdered.size() - 1] ), weight);
                hresponsedPhiJetsFB_2_Zinc2jet->Fill(deltaPhi(vJetYOrdered[0], vJetYOrdered[vJetYOrdered.size() - 1] ), deltaPhi(genvJetYOrdered[0], genvJetYOrdered[genvJetYOrdered.size() - 1] ), weight);
                
                hresponsedRJets_Zinc2jet->Fill(deltaRYPhi(newLeadJ, newSecondJ), deltaRYPhi(genNewLeadJ, genNewSecondJ), weight);
                hresponsedRJets_2_Zinc2jet->Fill(deltaRYPhi(newLeadJ, newSecondJ), deltaRYPhi(genNewLeadJ, genNewSecondJ), weight);
                
                hresponsediJetMass_Zinc2jet->Fill(jet1Plus2.M(), genJet1Plus2.M(), weight);
                hresponsediJetMass_2_Zinc2jet->Fill(jet1Plus2.M(), genJet1Plus2.M(), weight);
                
                hresponsediJetPt_Zinc2jet->Fill(jet1Plus2.Pt(), genJet1Plus2.Pt(), weight);
                hresponsediJetPt_2_Zinc2jet->Fill(jet1Plus2.Pt(), genJet1Plus2.Pt(), weight);
                
				hresponsedPhiLepJet2_Zinc2jet->Fill(deltaPhi(lep1, newSecondJ), deltaPhi(genLep1, genNewSecondJ), weight);
				hresponsedPhiLepJet2_2_Zinc2jet->Fill(deltaPhi(lep1, newSecondJ), deltaPhi(genLep1, genNewSecondJ), weight);
                
                hresponseSecondJetAbsRapidity_Zinc2jet->Fill(fabs(newSecondJ.Rapidity()), fabs(genNewSecondJ.Rapidity()), weight);
                hresponseSecondJetAbsRapidity_2_Zinc2jet->Fill(fabs(newSecondJ.Rapidity()), fabs(genNewSecondJ.Rapidity()), weight);


            }
            if (nGoodGenJetsAK8 >= 2 && nGoodJetsAK8 >= 2){
                hresponseLeadingJetAK8Pt_Zinc2jet->Fill(jetsAK8[0].pt, genJetsAK8[0].pt, weight);
                hresponseLeadingJetAK8Pt_2_Zinc2jet->Fill(jetsAK8[0].pt, genJetsAK8[0].pt, weight);
                hresponseLepPtPlusLeadingJetAK8Pt_Zinc2jet->Fill(lepton1.pt + jetsAK8[0].pt, genLep1.Pt()+genJetsAK8[0].pt, weight);
                hresponseLepPtPlusLeadingJetAK8Pt_2_Zinc2jet->Fill(lepton1.pt + jetsAK8[0].pt, genLep1.Pt()+genJetsAK8[0].pt, weight);
                if (nGoodGenJetsAK8 == 2 && nGoodJetsAK8 == 2){
                    hresponseLeadingJetAK8Pt_Zexc2jet->Fill(jetsAK8[0].pt, genJetsAK8[0].pt, weight);
                    hresponseLeadingJetAK8Pt_2_Zexc2jet->Fill(jetsAK8[0].pt, genJetsAK8[0].pt, weight);
                    hresponseLepPtPlusLeadingJetAK8Pt_Zexc2jet->Fill(lepton1.pt + jetsAK8[0].pt, genLep1.Pt()+genJetsAK8[0].pt, weight);
                    hresponseLepPtPlusLeadingJetAK8Pt_2_Zexc2jet->Fill(lepton1.pt + jetsAK8[0].pt, genLep1.Pt()+genJetsAK8[0].pt, weight);
                }
            }
            //-- Third Jet Pt
            if (nGoodGenJets >= 3 && nGoodJets >= 3){
                hresponseZNGoodJets_Zinc->Fill(3., 3., weight);
                hresponseZNGoodJetsFull_Zinc->Fill(3., 3., weight);
                hresponseThirdJetEta_Zinc3jet->Fill(fabs(jets[2].eta), fabs(genJets[2].eta), weight);
                hresponseJetsHT_Zinc3jet->Fill(jetsHT, genJetsHT, weight);
                hresponseJetsHT_20_30_Zinc3jet->Fill(jetsHT_20, genJetsHT_20, weight);
                
                hresponseThirdJetEta_2_Zinc3jet->Fill(fabs(jets[2].eta), fabs(genJets[2].eta), weight);
                hresponseJetsHT_1_Zinc3jet->Fill(jetsHT, genJetsHT, weight);
                hresponseJetsHT_2_Zinc3jet->Fill(jetsHT, genJetsHT, weight);
                
                hresponsedRapidityJets_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), weight);
                hresponsedRapidityJets_2_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), weight);
                
                hresponsedRapidityJetsFB_Zinc3jet->Fill(ForwardJetRapidity - BackwardJetRapidity, genForwardJetRapidity - genBackwardJetRapidity, weight);
                hresponsedRapidityJetsFB_2_Zinc3jet->Fill(ForwardJetRapidity - BackwardJetRapidity, genForwardJetRapidity - genBackwardJetRapidity, weight);
                
                hresponsedRapidityJets_First_Third_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newThirdJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewThirdJ.Rapidity()), weight);
                hresponsedRapidityJets_2_First_Third_Zinc3jet->Fill(fabs(newLeadJ.Rapidity() - newThirdJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewThirdJ.Rapidity()), weight);
                
                hresponsedRapidityJets_Second_Third_Zinc3jet->Fill(fabs(newSecondJ.Rapidity() - newThirdJ.Rapidity()), fabs(genNewSecondJ.Rapidity() - genNewThirdJ.Rapidity()), weight);
                hresponsedRapidityJets_2_Second_Third_Zinc3jet->Fill(fabs(newSecondJ.Rapidity() - newThirdJ.Rapidity()), fabs(genNewSecondJ.Rapidity() - genNewThirdJ.Rapidity()), weight);
                
                hresponsediJetMass_Zinc3jet->Fill(jet1Plus2.M(),genJet1Plus2.M(), weight);
                hresponsediJetMass_2_Zinc3jet->Fill(jet1Plus2.M(),genJet1Plus2.M(), weight);
                
                hresponsediJetPt_Zinc3jet->Fill(jet1Plus2.Pt(), genJet1Plus2.Pt(), weight);
                hresponsediJetPt_2_Zinc3jet->Fill(jet1Plus2.Pt(), genJet1Plus2.Pt(), weight);
                
				hresponsedPhiLepJet3_Zinc3jet->Fill(deltaPhi(lep1, newThirdJ), deltaPhi(genLep1, genNewThirdJ), weight);
				hresponsedPhiLepJet3_2_Zinc3jet->Fill(deltaPhi(lep1, newThirdJ), deltaPhi(genLep1, genNewThirdJ), weight);
                
                hresponseThirdJetAbsRapidity_Zinc3jet->Fill(fabs(newThirdJ.Rapidity()), fabs(genNewThirdJ.Rapidity()), weight);
                hresponseThirdJetAbsRapidity_2_Zinc3jet->Fill(fabs(newThirdJ.Rapidity()), fabs(genNewThirdJ.Rapidity()), weight);

            }
            if (nGoodGenJetsAK8 >= 3 && nGoodJetsAK8 >= 3){
                hresponseLeadingJetAK8Pt_Zinc3jet->Fill(jetsAK8[0].pt, genJetsAK8[0].pt, weight);
                hresponseLeadingJetAK8Pt_2_Zinc3jet->Fill(jetsAK8[0].pt, genJetsAK8[0].pt, weight);
                hresponseLepPtPlusLeadingJetAK8Pt_Zinc3jet->Fill(lepton1.pt + jetsAK8[0].pt, genLep1.Pt()+genJetsAK8[0].pt, weight);
                hresponseLepPtPlusLeadingJetAK8Pt_2_Zinc3jet->Fill(lepton1.pt + jetsAK8[0].pt, genLep1.Pt()+genJetsAK8[0].pt, weight);
                if (nGoodGenJetsAK8 == 3 && nGoodJetsAK8 == 3){
                    hresponseLeadingJetAK8Pt_Zexc3jet->Fill(jetsAK8[0].pt, genJetsAK8[0].pt, weight);
                    hresponseLeadingJetAK8Pt_2_Zexc3jet->Fill(jetsAK8[0].pt, genJetsAK8[0].pt, weight);
                    hresponseLepPtPlusLeadingJetAK8Pt_Zexc3jet->Fill(lepton1.pt + jetsAK8[0].pt, genLep1.Pt()+genJetsAK8[0].pt, weight);
                    hresponseLepPtPlusLeadingJetAK8Pt_2_Zexc3jet->Fill(lepton1.pt + jetsAK8[0].pt, genLep1.Pt()+genJetsAK8[0].pt, weight);
                }
            }
            //-- Fourth Jet Pt
            if (nGoodGenJets >= 4 && nGoodJets >= 4){
                hresponseZNGoodJets_Zinc->Fill(4., 4., weight);
                hresponseZNGoodJetsFull_Zinc->Fill(4., 4., weight);
                hresponseFourthJetEta_Zinc4jet->Fill(fabs(jets[3].eta), fabs(genJets[3].eta), weight);
                hresponseJetsHT_Zinc4jet->Fill(jetsHT, genJetsHT, weight);
                
                hresponseFourthJetEta_2_Zinc4jet->Fill(fabs(jets[3].eta), fabs(genJets[3].eta), weight);
                hresponseJetsHT_1_Zinc4jet->Fill(jetsHT, genJetsHT, weight);
                hresponseJetsHT_2_Zinc4jet->Fill(jetsHT, genJetsHT, weight);
                
                
                hresponsedRapidityJets_Zinc4jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), weight);
                hresponsedRapidityJets_2_Zinc4jet->Fill(fabs(newLeadJ.Rapidity() - newSecondJ.Rapidity()), fabs(genNewLeadJ.Rapidity() - genNewSecondJ.Rapidity()), weight);
                
                hresponsedRapidityJetsFB_Zinc4jet->Fill(ForwardJetRapidity - BackwardJetRapidity, genForwardJetRapidity - genBackwardJetRapidity, weight);
                hresponsedRapidityJetsFB_2_Zinc4jet->Fill(ForwardJetRapidity - BackwardJetRapidity, genForwardJetRapidity - genBackwardJetRapidity, weight);
                
                hresponsediJetMass_Zinc4jet->Fill(jet1Plus2.M(), genJet1Plus2.M(), weight);
                hresponsediJetMass_2_Zinc4jet->Fill(jet1Plus2.M(), genJet1Plus2.M(), weight);
                
                hresponsediJetPt_Zinc4jet->Fill(jet1Plus2.Pt(), genJet1Plus2.Pt(), weight);
                hresponsediJetPt_2_Zinc4jet->Fill(jet1Plus2.Pt(), genJet1Plus2.Pt(), weight);
                
                hresponsedPhiLepJet4_Zinc4jet->Fill(deltaPhi(lep1, newFourthJ), deltaPhi(genLep1, genNewFourthJ), weight);
                hresponsedPhiLepJet4_2_Zinc4jet->Fill(deltaPhi(lep1, newFourthJ), deltaPhi(genLep1, genNewFourthJ), weight);
                
                hresponseFourthJetAbsRapidity_Zinc4jet->Fill(fabs(newFourthJ.Rapidity()), fabs(genNewFourthJ.Rapidity()), weight);
                hresponseFourthJetAbsRapidity_2_Zinc4jet->Fill(fabs(newFourthJ.Rapidity()), fabs(genNewFourthJ.Rapidity()), weight);

            }
            //-- Fifth Jet Pt
            if (nGoodGenJets >= 5 && nGoodJets >= 5){
                hresponseZNGoodJets_Zinc->Fill(5., 5., weight);
                hresponseZNGoodJetsFull_Zinc->Fill(5., 5., weight);
                hresponseFifthJetEta_Zinc5jet->Fill(fabs(jets[4].eta), fabs(genJets[4].eta), weight);
                hresponseJetsHT_Zinc5jet->Fill(jetsHT, genJetsHT, weight);
                
                hresponseFifthJetEta_2_Zinc5jet->Fill(fabs(jets[4].eta), fabs(genJets[4].eta), weight);
                hresponseJetsHT_1_Zinc5jet->Fill(jetsHT, genJetsHT, weight);
                hresponseJetsHT_2_Zinc5jet->Fill(jetsHT, genJetsHT, weight);

				//hresponsedPhiLepJet5_Zinc5jet->Fill(deltaPhi(lep1, newFifthJ), deltaPhi(genLep1, genNewFifthJ), weight);
				//hresponsedPhiLepJet5_2_Zinc5jet->Fill(deltaPhi(lep1, newFifthJ), deltaPhi(genLep1, genNewFifthJ), weight);
            }
            //-- Sixth Jet Pt
            if (nGoodGenJets >= 6 && nGoodJets >= 6){
                hresponseZNGoodJets_Zinc->Fill(6., 6., weight);
                hresponseZNGoodJetsFull_Zinc->Fill(6., 6., weight);
            }
            //-- inc 7 jets case
            if (nGoodGenJets >= 7 && nGoodJets >= 7){
                hresponseZNGoodJets_Zinc->Fill(7., 7., weight);
                hresponseZNGoodJetsFull_Zinc->Fill(7., 7., weight);
            }
            //-- inc 8 jets case
            if (nGoodGenJets >= 8 && nGoodJets >= 8){
                hresponseZNGoodJets_Zinc->Fill(8., 8., weight);
            }
            //-- inc 9 jets case
            if (nGoodGenJets >= 9 && nGoodJets >= 9){
                hresponseZNGoodJets_Zinc->Fill(9., 9., weight);
            }
            //-- inc 10 jets case
            if (nGoodGenJets >= 10 && nGoodJets >= 10){
                hresponseZNGoodJets_Zinc->Fill(10., 10., weight);
            }
            
        }
        
        if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
        
        //=======================================================================================================//

    } //End of loop over all the events
    
    //==========================================================================================================//
	NEventsPassCuts->Fill(0., nEvents*1.0);
	NEventsPassCuts->Fill(1., nEventsPassMETFilter*1.0);
	NEventsPassCuts->Fill(2., countEventpassTrig*1.0);
	NEventsPassCuts->Fill(3., countEventpassLepReq*1.0);
	NEventsPassCuts->Fill(4., nEventsWithTwoGoodLeptons*1.0);
	NEventsPassCuts->Fill(5., countEventpassBveto*1.0);

    if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    //==========================================================================================================//
    //         Writing file             //
    //==================================//
    outputFile->cd();

    //--- Save all the histograms ---
    //listOfHistograms is defined in HistoSet.h
    unsigned short numbOfHistograms = listOfHistograms.size();

    cout << "\nskimAccep_[0] = " << skimAccep_[0] << "\n";
    cout << "sumEventW = " << sumEventW << "\n";
    if (hasRecoInfo && !isData){
        if (sumEventW > 0) {
            cout << "Scaling histograms by skimAccep_[0]/sumEventW:" << endl;
            cout << "skimAccep_[0]/sumEventW = " << skimAccep_[0]/sumEventW << "\n";
        }
    }
    
    for (unsigned short i(0); i < numbOfHistograms; i++){
        string hName = listOfHistograms[i]->GetName();
        if ( (!hasGenInfo && hName.find("gen") != string::npos ) || (!hasRecoInfo && hName.find("gen") == string::npos )) continue; 

        //global scaling of MC
        if (hasRecoInfo && !isData){
            if (sumEventW > 0) listOfHistograms[i]->Scale(skimAccep_[0]/sumEventW);
        }
        listOfHistograms[i]->Write();
        delete listOfHistograms[i];

        // need to make sure we don't skip over the deleting of the gen histos (for the reco files)
    }
    std::cout << "\nAll histograms written! Yay!" << std::endl;
    
    std::cout << "Writing the output file..." << std::endl;
    outputFile->Write();
    std::cout << "Closing the output file..." << std::endl;
    outputFile->Close();
    std::cout << "Output file " << outputFileName << " closed!" << std::endl;
    
    //==========================================================================================================//

    std::cout << "\nPrinting events' information:" << std::endl;
    cout << "Number of events                               : " << nEvents << endl;
    cout << "Number of events pass HasMET                   : " << countEvtpassHasMET << endl;
    cout << "Number of events pass MET Filter               : " << nEventsPassMETFilter << endl;
    cout << "Number of events pass trigger                  : " << countEventpassTrig << endl;
    cout << "Number of events pass Lepton requirements      : " << countEventpassLepReq << endl;
    cout << "Number of events pass MT cut                   : " << nEventsWithTwoGoodLeptons << endl;
    cout << "Number of events pass Btag veto                : " << countEventpassBveto << endl;
    cout << "Total GEN weight of all events                 : " << TotalGenWeight << endl;
    cout << "Total GEN pass: RECO weight of all events      : " << TotalGenWeightPassGENPU << endl;
    cout << "Total GEN pass: GEN weight of all events       : " << TotalGenWeightPassGEN << endl;
    cout << "Total RECO pass: RECO weight of all events     : " << TotalRecoWeightPassRECO << endl;
    cout << "Total RECO pass: GEN weight of all events      : " << TotalGenWeightPassRECO << endl;
    cout << "# Events: 0 jets inclusive                     : " << nEventsIncl0Jets << endl;
    cout << "# Events: 0 jets exclusive                     : " << nEventsExcl0Jets << endl;
    cout << "# Events: 1 jet exclusive                      : " << nEventsExcl1Jets << endl;
    cout << "# Events: 2 jets exclusive                     : " << nEventsExcl2Jets << endl;
    cout << "# Events: 3 jets exclusive                     : " << nEventsExcl3Jets << endl;
    cout << "# Events: 1 B-jet inclusive                    : " << nEventsIncBJets << endl;
    cout << "# GEN Events: 0 jets inclusive                 : " << GENnEventsIncl0Jets << endl;
    cout << "# GEN Events: 1 jets inclusive                 : " << GENnEventsIncl1Jets << endl;
    cout << "# GEN Events: 2 jets inclusive                 : " << GENnEventsIncl2Jets << endl;
    cout << "# GEN Events: 3 jets inclusive                 : " << GENnEventsIncl3Jets << endl;
    cout << "Total MC weight (sumEventW)                    : " << sumEventW << endl;
    std::cout << "\n=======================================================================================================" << std::endl;
}


ZJetsAndDPS::ZJetsAndDPS(string fileName_, int year_, float lumiScale_, float puScale_, bool useTriggerCorrection_, bool useEfficiencyCorrection_, 
        int systematics_, int direction_, float xsecfactor_, int jetPtCutMin_, int jetPtCutMax_, int ZPtCutMin_, int ZEtaCutMin_, int ZEtaCutMax_, 
        int METcut_, int jetEtaCutMin_, int jetEtaCutMax_): 

    HistoSet(fileName_.substr(0, fileName_.find("_"))), outputDirectory("HistoFiles/"),
    fileName(fileName_), year(year_), lumiScale(lumiScale_), puScale(puScale_), 
    useTriggerCorrection(useTriggerCorrection_), useEfficiencyCorrection(useEfficiencyCorrection_), 
    systematics(systematics_), direction(direction_), xsecfactor(xsecfactor_), 
    jetPtCutMin(jetPtCutMin_), jetPtCutMax(jetPtCutMax_), jetEtaCutMin(jetEtaCutMin_), jetEtaCutMax(jetEtaCutMax_), 
    ZPtCutMin(ZPtCutMin_), ZEtaCutMin(ZEtaCutMin_), ZEtaCutMax(ZEtaCutMax_), METcut(METcut_)

{
    std::cout << "\n >>>>>>>>>> ZJetsAndDPS::ZJetsAndDPS() >>>>>>>>>> " << std::endl;

    //This function serves to get the correct filenames/locations, and then read in the files to the TChain

    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.

    TChain *chain = new TChain("", "");
    TChain *BonzaiHeaderChain = new TChain("", "");
    
    isData = false;
    string fullFileName =  "../Data_Z_5311_New/" + fileName;
    
    string storageElement = "/eos/cms";
    string dirPath = "/tupel";
    string treeName = "/EventTree";

    std::string yearStr;
    std::stringstream yearSStr;
    yearSStr << year;
    yearStr = yearSStr.str();

    if (fileName.find("DMu_") == 0) leptonFlavor = "Muons";
    else if (fileName.find("DE_") == 0)  leptonFlavor = "Electrons"; 
    else if (fileName.find("SMu_") == 0) leptonFlavor = "SingleMuon";
    else if (fileName.find("SE_") == 0)  leptonFlavor = "SingleElectron";
    else if (fileName.find("SMuE_") == 0){
        leptonFlavor = "TTMuE";
        fullFileName =  "../DataTTbarEMu/" + fileName;
    }
    if (fileName.find("Data") != string::npos ) isData = true;
    if (fileName.find("SMu_") == 0 || fileName.find("SE_") == 0 ) fullFileName =  "DataW_txt_" + yearStr + "/" + fileName;
    if (fileName.find("Sherpa2") != string::npos) fullFileName =  "../DataSherpa2/" + fileName;
    if (fileName.find("List") == string::npos){
        if (fileName.find("Sherpa2") != string::npos){
            fullFileName += ".root";
            string treePath = fullFileName + "/tree";
            cout << "Loading file: " << fullFileName << endl;
            chain->Add(treePath.c_str());
        }
        else{
            fullFileName += ".root";
            string treePath = fullFileName + "/tupel/EventTree";
            cout << "Loading file: " << fullFileName << endl;
            chain->Add(treePath.c_str());
        }
    }
    else {
        //assuming that the filename locations are listed in a text file, which are read off from EOS
        fullFileName += ".txt";
        std::cout << "\nReading in .txt file: \n" << fullFileName << std::endl;
        ifstream infile(fullFileName.c_str());
        string line; 
        int countFiles(0);
        // std::cout << __LINE__ << std::endl;
        std::cout << std::endl;
        while (getline(infile, line)){
            
            // Loading the input files into TChains -----
            string treePath = storageElement + line + dirPath + treeName;
            string bonzaiHeaderPath = storageElement + line + dirPath + "/BonzaiHeader";
            cout << "Loading file #" << countFiles << ": " << line << endl;
            chain->Add(treePath.c_str());
            BonzaiHeaderChain->Add(bonzaiHeaderPath.c_str());
        
            countFiles++;
        }
    }
    fChain = chain;
    fBonzaiHeaderChain = BonzaiHeaderChain;
    
    if (!isData) getMcNorm();
    else skimAccep_ = std::vector<double>(1., 1.);
}

ZJetsAndDPS::~ZJetsAndDPS(){
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

string ZJetsAndDPS::CreateOutputFileName(bool useRoch, bool doFlat, int doPUStudy, bool doVarWidth, int doBJets, int doQCD, bool doSSign, bool doInvMassCut){
    ostringstream result;
    result << outputDirectory << fileName;
    result << "_EffiCorr_" << useEfficiencyCorrection;
    result << "_TrigCorr_" << useTriggerCorrection;
    result << "_Syst_" << systematics;
    if (direction == 1) result << "_Up";
    else if (direction == -1) result << "_Down";
    result << "_JetPtMin_" << jetPtCutMin;
    if (jetPtCutMax > jetPtCutMin) result << "_JetPtMax_" << jetPtCutMax;
    if (ZPtCutMin > 0) result << "_ZPtMin" << abs(ZPtCutMin);
    if (ZEtaCutMin > -999999 && ZEtaCutMin <  0) result << "_ZEtaMin_m" << abs(ZEtaCutMin);
    if (ZEtaCutMin > -999999 && ZEtaCutMin >= 0) result << "_ZEtaMin_"  << abs(ZEtaCutMin);
    if (ZEtaCutMax <  999999 && ZEtaCutMax >= 0) result << "_ZEtaMax_"  << abs(ZEtaCutMax);
    if (ZEtaCutMax <  999999 && ZEtaCutMax <  0) result << "_ZEtaMax_m" << abs(ZEtaCutMax);

    if (useRoch) result << "_rochester";
    if (!isData && doFlat) result << "_Flat";
    if (doPUStudy >= 0) result << "_Beta" << doPUStudy;
    if (doVarWidth) result << "_VarWidth";
    if (doInvMassCut) result << "_InvMass";
    if (doSSign) result << "_SS";
    if (doBJets > 0) result << "_BJets";
    if (doBJets < 0) result << "_BVeto";
    if (doQCD>0) result << "_QCD" << doQCD;
    if (METcut > 0) result << "_MET" << METcut;

    result << ".root";
    return result.str();
}

Int_t ZJetsAndDPS::GetEntry(Long64_t entry){
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

Long64_t ZJetsAndDPS::LoadTree(Long64_t entry){
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
        Notify();  
    }
    return centry;
}

void ZJetsAndDPS::Init(bool hasRecoInfo, bool hasGenInfo){
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).


    // Set object pointers -----

    EvtWeights = 0;
    mcSherpaWeights_ = 0 ;

    // TrigHltMu = 0;
    // TrigMET = 0;
    // TrigMETBit = 0;
    
    /////////////////////////////
    // GENERATOR VARIABLES

    GLepBarePt = 0;
    GLepBareEta = 0;
    GLepBarePhi = 0;
    GLepBareE = 0;
    GLepBareId = 0;
    GLepBareSt = 0;
    GLepBarePrompt = 0;

    GLepClosePhotPt = 0;
    GLepClosePhotEta = 0;
    GLepClosePhotPhi = 0;
    // GLepClosePhotE = 0;
    // GLepClosePhotM = 0;
    // GLepClosePhotId = 0;
    // GLepClosePhotMother0Id = 0;
    // GLepClosePhotMotherCnt = 0;
    GLepClosePhotSt = 0;

    // GMETPt = 0;
    // GMETPx = 0;
    // GMETPy = 0;
    // GMETE = 0;
    // GMETPhi = 0;

    GJetAk04Pt = 0;
    GJetAk04Eta = 0;
    GJetAk04Phi = 0;
    GJetAk04E = 0;

    GJetAk08Pt = 0;
    GJetAk08Eta = 0;
    GJetAk08Phi = 0;
    GJetAk08E = 0;

    /////////////////////////////
    // RECO VARIABLES

    MuPt = 0;
    MuEta = 0;
    MuPhi = 0;
    MuE = 0;
    // MuIdLoose = 0;
    // MuIdMedium = 0;
    MuIdTight = 0;
    MuCh = 0;
    // MuVtxZ = 0;
    // MuDxy = 0;
    MuPfIso = 0;
    // MuDz = 0;
    MuHltTrgPath1 = 0;
    MuHltTrgPath2 = 0;
    
    METPt = 0;
    METPx = 0;
    METPy = 0;
    METE = 0;
    METPhi = 0;
    METFilterPath1 = 0;
    METFilterPath2 = 0;
    METFilterPath3 = 0;
    METFilterPath4 = 0;
    METFilterPath5 = 0;
    METFilterPath6 = 0;
    METFilterPath7 = 0;
    
    JetAk04Pt = 0;   
    JetAk04Eta = 0;   
    JetAk04Phi = 0;   
    JetAk04E = 0;   
    JetAk04Id = 0;   
    JetAk04PuId = 0; 
    JetAk04PuIdLoose = 0; 
    JetAk04PuIdMedium = 0; 
    JetAk04PuIdTight = 0; 
    JetAk04PuMva = 0;   
    JetAk04BDiscCisvV2 = 0;
    JetAk04HadFlav = 0;
    JetAk04JecUncUp = 0;
    JetAk04JecUncDwn = 0;

    JetAk08Pt = 0;   
    JetAk08Eta = 0;   
    JetAk08Phi = 0;   
    JetAk08E = 0;   
    JetAk08Id = 0;     
    JetAk08BDiscCisvV2 = 0;
    JetAk08HadFlav = 0;
    JetAk08JecUncUp = 0;
    JetAk08JecUncDwn = 0;

    // Set branch addresses and branch pointers -----

    if (DEBUG) cout << "Stop after line " << __LINE__ << endl;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    if (hasRecoInfo){

        // General variables
        // used for pileup, pileup reweighting for MC, etc.
        fChain->SetBranchAddress("EvtVtxCnt", &EvtVtxCnt, &b_EvtVtxCnt);
        fChain->SetBranchAddress("EvtPuCntObs", &EvtPuCntObs, &b_EvtPuCntObs);
        fChain->SetBranchAddress("EvtPuCntTruth", &EvtPuCntTruth, &b_EvtPuCntTruth);
        // contains main weights for MC
        fChain->SetBranchAddress("EvtWeights", &EvtWeights, &b_EvtWeights);
        // old method for storing trigger bit information
        // fChain->SetBranchAddress("TrigHltMu", &TrigHltMu, &b_TrigHltMu);
        // fChain->SetBranchAddress("TrigMET", &TrigMET, &b_TrigMET);
        // fChain->SetBranchAddress("TrigMETBit", &TrigMETBit, &b_TrigMETBit);
        
        // Muons, Muon HLT Trigger Paths
        fChain->SetBranchAddress("MuPt", &MuPt, &b_MuPt);
        fChain->SetBranchAddress("MuEta", &MuEta, &b_MuEta);
        fChain->SetBranchAddress("MuPhi", &MuPhi, &b_MuPhi);
        fChain->SetBranchAddress("MuE", &MuE, &b_MuE);
        // fChain->SetBranchAddress("MuIdLoose", &MuIdLoose, &b_MuIdLoose);
        // fChain->SetBranchAddress("MuIdMedium", &MuIdMedium, &b_MuIdMedium);
        fChain->SetBranchAddress("MuIdTight", &MuIdTight, &b_MuIdTight);
        fChain->SetBranchAddress("MuCh", &MuCh, &b_MuCh);
        // fChain->SetBranchAddress("MuVtxZ", &MuVtxZ, &b_MuVtxZ);
        // fChain->SetBranchAddress("MuDxy", &MuDxy, &b_MuDxy);
        fChain->SetBranchAddress("MuPfIso", &MuPfIso, &b_MuPfIso);
        // fChain->SetBranchAddress("MuDz", &MuDz, &b_MuDz);
        fChain->SetBranchAddress("MuHltTrgPath1", &MuHltTrgPath1, &b_MuHltTrgPath1);
        fChain->SetBranchAddress("MuHltTrgPath2", &MuHltTrgPath2, &b_MuHltTrgPath2);
        
        //MET, MET filters
        fChain->SetBranchAddress("METPt", &METPt, &b_METPt);
        fChain->SetBranchAddress("METPx", &METPx, &b_METPx);
        fChain->SetBranchAddress("METPy", &METPy, &b_METPy);
        fChain->SetBranchAddress("METE", &METE, &b_METE);
        fChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
        fChain->SetBranchAddress("METFilterPath1", &METFilterPath1, &b_METFilterPath1);
        fChain->SetBranchAddress("METFilterPath2", &METFilterPath2, &b_METFilterPath2);
        fChain->SetBranchAddress("METFilterPath3", &METFilterPath3, &b_METFilterPath3);
        fChain->SetBranchAddress("METFilterPath4", &METFilterPath4, &b_METFilterPath4);
        fChain->SetBranchAddress("METFilterPath5", &METFilterPath5, &b_METFilterPath5);
        fChain->SetBranchAddress("METFilterPath6", &METFilterPath6, &b_METFilterPath6);
        fChain->SetBranchAddress("METFilterPath7", &METFilterPath7, &b_METFilterPath7);
        
        //PF Jets - AK4
        fChain->SetBranchAddress("JetAk04Pt", &JetAk04Pt, &b_JetAk04Pt);
        fChain->SetBranchAddress("JetAk04Eta", &JetAk04Eta, &b_JetAk04Eta);
        fChain->SetBranchAddress("JetAk04Phi", &JetAk04Phi, &b_JetAk04Phi);
        fChain->SetBranchAddress("JetAk04E", &JetAk04E, &b_JetAk04E);
        fChain->SetBranchAddress("JetAk04Id", &JetAk04Id, &b_JetAk04Id);
        fChain->SetBranchAddress("JetAk04PuId", &JetAk04PuId, &b_JetAk04PuId);
        fChain->SetBranchAddress("JetAk04PuIdLoose", &JetAk04PuIdLoose, &b_JetAk04PuIdLoose);
        fChain->SetBranchAddress("JetAk04PuIdMedium", &JetAk04PuIdMedium, &b_JetAk04PuIdMedium);
        fChain->SetBranchAddress("JetAk04PuIdTight", &JetAk04PuIdTight, &b_JetAk04PuIdTight);
        fChain->SetBranchAddress("JetAk04PuMva", &JetAk04PuMva, &b_JetAk04PuMva);
        fChain->SetBranchAddress("JetAk04BDiscCisvV2", &JetAk04BDiscCisvV2, &b_JetAk04BDiscCisvV2); 
	    fChain->SetBranchAddress("JetAk04HadFlav", &JetAk04HadFlav, &b_JetAk04HadFlav); 
        fChain->SetBranchAddress("JetAk04JecUncUp", &JetAk04JecUncUp, &b_JetAk04JecUncUp); 
        fChain->SetBranchAddress("JetAk04JecUncDwn", &JetAk04JecUncDwn, &b_JetAk04JecUncDwn); 

        //PF Jets - AK8
        fChain->SetBranchAddress("JetAk08Pt", &JetAk08Pt, &b_JetAk08Pt);
        fChain->SetBranchAddress("JetAk08Eta", &JetAk08Eta, &b_JetAk08Eta);
        fChain->SetBranchAddress("JetAk08Phi", &JetAk08Phi, &b_JetAk08Phi);
        fChain->SetBranchAddress("JetAk08E", &JetAk08E, &b_JetAk08E);
        fChain->SetBranchAddress("JetAk08Id", &JetAk08Id, &b_JetAk08Id);
        fChain->SetBranchAddress("JetAk08BDiscCisvV2", &JetAk08BDiscCisvV2, &b_JetAk08BDiscCisvV2); 
	    fChain->SetBranchAddress("JetAk08HadFlav", &JetAk08HadFlav, &b_JetAk08HadFlav); 
        fChain->SetBranchAddress("JetAk08JecUncUp", &JetAk08JecUncUp, &b_JetAk08JecUncUp); 
        fChain->SetBranchAddress("JetAk08JecUncDwn", &JetAk08JecUncDwn, &b_JetAk08JecUncDwn); 

    } //end hasRecoInfo

    if (hasGenInfo){

        //Generator level leptons, not-dressed
        fChain->SetBranchAddress("GLepBarePt", &GLepBarePt, &b_GLepBarePt);
        fChain->SetBranchAddress("GLepBareEta", &GLepBareEta, &b_GLepBareEta);
        fChain->SetBranchAddress("GLepBarePhi", &GLepBarePhi, &b_GLepBarePhi);
        fChain->SetBranchAddress("GLepBareE", &GLepBareE, &b_GLepBareE);
        fChain->SetBranchAddress("GLepBareId", &GLepBareId, &b_GLepBareId);
        fChain->SetBranchAddress("GLepBareSt", &GLepBareSt, &b_GLepBareSt);
        fChain->SetBranchAddress("GLepBarePrompt", &GLepBarePrompt, &b_GLepBarePrompt);

        //Photons in the vicinity of the leptons
        fChain->SetBranchAddress("GLepClosePhotPt", &GLepClosePhotPt, &b_GLepClosePhotPt);
        fChain->SetBranchAddress("GLepClosePhotEta", &GLepClosePhotEta, &b_GLepClosePhotEta);
        fChain->SetBranchAddress("GLepClosePhotPhi", &GLepClosePhotPhi, &b_GLepClosePhotPhi);
        // fChain->SetBranchAddress("GLepClosePhotE", &GLepClosePhotE, &b_GLepClosePhotE);
        // fChain->SetBranchAddress("GLepClosePhotM", &GLepClosePhotM, &b_GLepClosePhotM);
        // fChain->SetBranchAddress("GLepClosePhotId", &GLepClosePhotId, &b_GLepClosePhotId);
        // fChain->SetBranchAddress("GLepClosePhotMother0Id", &GLepClosePhotMother0Id, &b_GLepClosePhotMother0Id);
        // fChain->SetBranchAddress("GLepClosePhotMotherCnt", &GLepClosePhotMotherCnt, &b_GLepClosePhotMotherCnt);
        fChain->SetBranchAddress("GLepClosePhotSt", &GLepClosePhotSt, &b_GLepClosePhotSt);

        //GEN-level MET
        // fChain->SetBranchAddress("GMETPt", &GMETPt, &b_GMETPt);
        // fChain->SetBranchAddress("GMETPx", &GMETPx, &b_GMETPx);
        // fChain->SetBranchAddress("GMETPy", &GMETPy, &b_GMETPy);
        // fChain->SetBranchAddress("GMETE", &GMETE, &b_GMETE);
        // fChain->SetBranchAddress("GMETPhi", &GMETPhi, &b_GMETPhi);

        //Gen Jets, AK4
        fChain->SetBranchAddress("GJetAk04Pt", &GJetAk04Pt, &b_GJetAk04Pt);
        fChain->SetBranchAddress("GJetAk04Eta", &GJetAk04Eta, &b_GJetAk04Eta);
        fChain->SetBranchAddress("GJetAk04Phi", &GJetAk04Phi, &b_GJetAk04Phi);
        fChain->SetBranchAddress("GJetAk04E", &GJetAk04E, &b_GJetAk04E);

        //Gen Jets, AK8
        fChain->SetBranchAddress("GJetAk08Pt", &GJetAk08Pt, &b_GJetAk08Pt);
        fChain->SetBranchAddress("GJetAk08Eta", &GJetAk08Eta, &b_GJetAk08Eta);
        fChain->SetBranchAddress("GJetAk08Phi", &GJetAk08Phi, &b_GJetAk08Phi);
        fChain->SetBranchAddress("GJetAk08E", &GJetAk08E, &b_GJetAk08E);

    } //end hasGenInfo

    Notify();
    cout << "Branches are properly initialized." << endl;
}

Bool_t ZJetsAndDPS::Notify(){
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void ZJetsAndDPS::Show(Long64_t entry){
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}

Int_t ZJetsAndDPS::Cut(Long64_t entry){
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    printf("entry %lld", entry);
    return 1;
}

void ZJetsAndDPS::getMcNorm(){

    std::cout << "\n >>>>>>>>>> ZJetsAndDPS::getMcNorm() >>>>>>>>>> " << std::endl;
    Int_t InEvtCount = 0;
    std::vector<Double_t>* InEvtWeightSums  = 0;
    std::vector<Double_t>* EvtWeightSums = 0;
    fBonzaiHeaderChain->SetBranchAddress("InEvtWeightSums", &InEvtWeightSums);
    fBonzaiHeaderChain->SetBranchAddress("EvtWeightSums", &EvtWeightSums);

    int nheaders = fBonzaiHeaderChain->GetEntries(); //can be several in case files were merged with hadd
    cout << " >>>>> fBonzaiHeaderChain: Looping over " << nheaders << " TChain entries!" << endl;
    for(int ientry = 0; ientry < nheaders; ++ientry){
        cout << " -----> fBonzaiHeaderChain: getting entry #" << ientry << endl;
        fBonzaiHeaderChain->GetEntry(ientry);
        
        if(ientry == 0){
            InEvtWeightSums_ = std::vector<Double_t>(InEvtWeightSums->size(), 0);
            EvtWeightSums_ = std::vector<Double_t>(EvtWeightSums->size(), 0);
            if(InEvtWeightSums->size() != EvtWeightSums->size()){
                std::cerr << "InEvtWeightSums and EvtWeightSums branches "
                "of input BonzaiHeader tree have different size ("
                "resp. " << InEvtWeightSums->size() << " and "
                << EvtWeightSums->size() << ")\n";
                abort();
            }
        }
        if(InEvtWeightSums->size() != InEvtWeightSums_.size()){
            std::cerr << "Inconsistency in number of elements of InEvtWeightSums branch of input files!\n";
            abort();
        }
        if(EvtWeightSums->size() != EvtWeightSums_.size()){
            std::cerr << "Inconsistency in number of elements of EvtWeightSums branch of input files!\n";
            abort();
        }
        
        // Summing elements of InEvtWeightSums and EvtWeightSums vectors over all events
        for(size_t i = 0; i < InEvtWeightSums_.size(); ++i){
            InEvtWeightSums_[i] += (*InEvtWeightSums)[i];
        }
        for(size_t i = 0; i < EvtWeightSums_.size(); ++i){
            EvtWeightSums_[i] += (*EvtWeightSums)[i];
        }

        InEvtCount_ += InEvtCount;

        // Print out InEvtWeightSums_[0] and EvtWeightSums_[0] -----
        // we can see, per event, if there's an inf/NaN value here
        cout << "(*EvtWeightSums)[0]   = " << (*EvtWeightSums)[0] << endl;
        cout << "(*InEvtWeightSums)[0] = " << (*InEvtWeightSums)[0] << endl;
        // cout << "(*EvtWeightSums)[0]/(*InEvtWeightSums)[0] = " << ((*EvtWeightSums)[0])/((*InEvtWeightSums)[0]) << endl;
    }
    
    if(InEvtWeightSums_.size() > 0){
        cout << "\n >>>>> Sum totals over all Bonzai header events ----- " << endl;
        std::cout << "EvtWeightSums_[0] = " <<  EvtWeightSums_[0] << std::endl;
        std::cout << "InEvtWeightSums_[0] = " << InEvtWeightSums_[0] << std::endl;
    }
    
    if(EvtWeightSums_.size() == 0 || InEvtWeightSums_.size() == 0 || InEvtWeightSums_[0] == 0){
        EvtCount_ = fChain->GetEntries();
        if(InEvtCount_){
            skimAccep_ = std::vector<double>(1, EvtCount_/InEvtCount_);
        } 
        else{
            std::cout << "Warning: InEvtCount is equal to 0. Event yield normalization might be wrong!" << std::endl;
        }
    } 
    //andrew -- currently what our code uses for determining Bonzai acceptance -- 3 sept 2019
    else{
        skimAccep_ = std::vector<double>(InEvtWeightSums_.size());
        for(unsigned i = 0; i < InEvtWeightSums_.size(); ++i){
            skimAccep_[i] = EvtWeightSums_[i]/InEvtWeightSums_[i];
        }
    }

    cout << "\n >>>>> Sum totals over all Bonzai header events ----- " << endl;
    std::cerr << "skimAccep_[0] = " << skimAccep_[0] << "\n";
}
