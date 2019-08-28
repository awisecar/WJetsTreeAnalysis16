//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 22 08:07:46 2012 by ROOT version 5.32/01
// from TTree tree/tree
//////////////////////////////////////////////////////////

#ifndef _ZJetsAndDPS_H_
#define _ZJetsAndDPS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TDatime.h>
#include <TMath.h>
#include <TRandom3.h>

// Header file for the classes stored in the TTree if any.
#include "functions.h"
#include "getFilesAndHistograms.h"
#include "standalone_LumiReWeighting.h"
#include "HistoSet.h"

using namespace std;

// andrew - is this line necessary?
// ClassImp(HistoSet) 

class ZJetsAndDPS: public HistoSet {
  public:

    // Constructor and destructor --------------------
    ZJetsAndDPS(string fileName_, int year_ = 2017, float lumiScale_ = 1., float puScale_ = 1., bool useTriggerCorrection_ = 0, bool useEfficiencyCorrection_ = 0, int systematics_ = 0, int direction_ = 0, float xsecfactor_ = 1., int jetPtCutMin_ = 20, int jetPtCutMax_ = 0, int ZPtCutMin_ = 0 , int ZEtaCutMin_ = -999999, int ZEtaCutMax_ = 999999, int METcut_ = -30, int jetEtaCutMin_ = -24, int jetEtaCutMax_ = 24); 
    ~ZJetsAndDPS();

    //TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    TChain          *fChain;   //!pointer to the analyzed TTree or TChain
    TChain          *fBonzaiHeaderChain;
    Int_t           fCurrent; //!current Tree number in a TChain

    // Declaration of leaf types --------------------
    Int_t           EvtNum;
    Int_t           EvtRunNum;
    Int_t           EvtVtxCnt;
    Int_t           EvtPuCntObs;
    Int_t           EvtPuCntTruth;
    vector<double>  *EvtWeights;

    TBranch        *b_EvtNum;   //!
    TBranch        *b_EvtRunNum;   //!
    TBranch        *b_EvtVtxCnt;   //!
    TBranch        *b_EvtPuCntObs;   //!
    TBranch        *b_EvtPuCntTruth;   //!
    TBranch        *b_EvtWeights;   //!

    Int_t           GNup;
    Double_t        mcEveWeight_;
    Double_t        mcSherpaSumWeight3_ ;
    vector<double> *mcSherpaWeights_;

    TBranch        *b_GNup;   //!
    TBranch        *b_mcEveWeight_;   //!
    TBranch        *b_mcSherpaSumWeight3_;   //!
    TBranch        *b_mcSherpaWeights_;   //!

    ULong64_t      TrigHltMu;
    ULong64_t      TrigMET;
    ULong64_t      TrigMETBit;

    TBranch        *b_TrigHltMu;   //!
    TBranch        *b_TrigMET;   //!
    TBranch        *b_TrigMETBit;   //!
    
    // GENERATOR LEVEL (hasGenInfo) --------------------
    //Generator level leptons, not-dressed
    vector<float>  *GLepBarePt;
    vector<float>  *GLepBareEta;
    vector<float>  *GLepBarePhi;
    vector<float>  *GLepBareE;
    vector<int>    *GLepBareId;
    vector<int>    *GLepBareSt;
    vector<bool>   *GLepBarePrompt;

    TBranch        *b_GLepBarePt;   //!
    TBranch        *b_GLepBareEta;   //!
    TBranch        *b_GLepBarePhi;   //!
    TBranch        *b_GLepBareE;   //!
    TBranch        *b_GLepBareId;   //!
    TBranch        *b_GLepBareSt;   //!
    TBranch        *b_GLepBarePrompt;   //!

    //Photons in the vicinity of the leptons
    vector<float>  *GLepClosePhotPt;
    vector<float>  *GLepClosePhotEta;
    vector<float>  *GLepClosePhotPhi;
    vector<float>  *GLepClosePhotE;
    vector<float>  *GLepClosePhotM;
    vector<int>    *GLepClosePhotId;
    vector<int>    *GLepClosePhotMother0Id;
    vector<int>    *GLepClosePhotMotherCnt;
    vector<int>    *GLepClosePhotSt;

    TBranch        *b_GLepClosePhotPt;   //!
    TBranch        *b_GLepClosePhotEta;   //!
    TBranch        *b_GLepClosePhotPhi;   //!
    TBranch        *b_GLepClosePhotE;   //!
    TBranch        *b_GLepClosePhotM;   //!
    TBranch        *b_GLepClosePhotId;   //!
    TBranch        *b_GLepClosePhotMother0Id;   //!
    TBranch        *b_GLepClosePhotMotherCnt;   //!
    TBranch        *b_GLepClosePhotSt;   //!

    // GEN-level MET
    vector<float>  *GMETPt;
    vector<float>  *GMETPx;
    vector<float>  *GMETPy;
    vector<float>  *GMETE;
    vector<float>  *GMETPhi;

    TBranch        *b_GMETPt;
    TBranch        *b_GMETPx;
    TBranch        *b_GMETPy;
    TBranch        *b_GMETE;
    TBranch        *b_GMETPhi;   
   
    //Gen Jets, AK4
    vector<float>  *GJetAk04Pt;
    vector<float>  *GJetAk04Eta;
    vector<float>  *GJetAk04Phi;
    vector<float>  *GJetAk04E;

    TBranch        *b_GJetAk04Pt;   //!
    TBranch        *b_GJetAk04Eta;   //!
    TBranch        *b_GJetAk04Phi;   //!
    TBranch        *b_GJetAk04E;   //!

    //Gen Jets, AK8
    vector<float>  *GJetAk08Pt;
    vector<float>  *GJetAk08Eta;
    vector<float>  *GJetAk08Phi;
    vector<float>  *GJetAk08E;

    TBranch        *b_GJetAk08Pt;   //!
    TBranch        *b_GJetAk08Eta;   //!
    TBranch        *b_GJetAk08Phi;   //!
    TBranch        *b_GJetAk08E;   //!

    // RECONSTRUCTED LEVEL (hasRecoInfo) --------------------
    // Muons, Muon HLT Trigger Paths
    vector<float>        *MuPt;
    vector<float>        *MuEta;
    vector<float>        *MuPhi;
    vector<float>        *MuE;
    vector<bool>         *MuIdLoose;
    vector<bool>         *MuIdMedium;
    vector<bool>         *MuIdTight;
    vector<float>        *MuCh;
    vector<float>        *MuVtxZ;
    vector<float>        *MuDxy;
    vector<float>        *MuPfIso;
    vector<float>        *MuDz;
    vector<bool>         *MuHltTrgPath1;
    vector<bool>         *MuHltTrgPath2;

    TBranch        *b_MuPt;   //!
    TBranch        *b_MuEta;   //!
    TBranch        *b_MuPhi;   //!
    TBranch        *b_MuE;   //!
    TBranch        *b_MuIdLoose;   //!
    TBranch        *b_MuIdMedium;   //!
    TBranch        *b_MuIdTight;   //!
    TBranch        *b_MuCh;   //!
    TBranch        *b_MuVtxZ;   //!
    TBranch        *b_MuDxy;   //!
    TBranch        *b_MuPfIso;   //!
    TBranch        *b_MuDz;   //!
    TBranch        *b_MuHltTrgPath1;   //!
    TBranch        *b_MuHltTrgPath2;   //!
    
    //MET, MET filters
    vector<float>  *METPt;
    vector<float>  *METPx;
    vector<float>  *METPy;
    vector<float>  *METE;
    vector<float>  *METPhi;
    vector<bool>   *METFilterPath1;
    vector<bool>   *METFilterPath2;
    vector<bool>   *METFilterPath3;
    vector<bool>   *METFilterPath4;
    vector<bool>   *METFilterPath5;
    vector<bool>   *METFilterPath6;
    vector<bool>   *METFilterPath7;

    TBranch        *b_METPt;   //!
    TBranch        *b_METPx;   //!
    TBranch        *b_METPy;   //!
    TBranch        *b_METE;   //!
    TBranch        *b_METPhi;   //!
    TBranch        *b_METFilterPath1;   //!
    TBranch        *b_METFilterPath2;   //!
    TBranch        *b_METFilterPath3;   //!
    TBranch        *b_METFilterPath4;   //!
    TBranch        *b_METFilterPath5;   //!
    TBranch        *b_METFilterPath6;   //!
    TBranch        *b_METFilterPath7;   //!
    
    //PF Jets - AK4
    vector<float>  *JetAk04Pt;
    vector<float>  *JetAk04Eta;
    vector<float>  *JetAk04Phi;
    vector<float>  *JetAk04E;
    vector<float>  *JetAk04Id;
    vector<int>    *JetAk04PuId;
    vector<bool>   *JetAk04PuIdLoose;
    vector<bool>   *JetAk04PuIdMedium;
    vector<bool>   *JetAk04PuIdTight;
    vector<float>  *JetAk04PuMva;
    vector<float>  *JetAk04BDiscCisvV2;
    vector<float>  *JetAk04HadFlav;
    vector<float>  *JetAk04JecUncUp;
    vector<float>  *JetAk04JecUncDwn;

    TBranch        *b_JetAk04Pt;   
    TBranch        *b_JetAk04Eta;   
    TBranch        *b_JetAk04Phi;   
    TBranch        *b_JetAk04E;  
    TBranch        *b_JetAk04Id;   
    TBranch        *b_JetAk04PuId;
    TBranch        *b_JetAk04PuIdLoose;
    TBranch        *b_JetAk04PuIdMedium;
    TBranch        *b_JetAk04PuIdTight; 
    TBranch        *b_JetAk04PuMva;   
    TBranch        *b_JetAk04BDiscCisvV2;   
    TBranch        *b_JetAk04HadFlav;  
    TBranch        *b_JetAk04JecUncUp;
    TBranch        *b_JetAk04JecUncDwn;
    
    //PF Jets - AK8
    vector<float>  *JetAk08Pt;
    vector<float>  *JetAk08Eta;
    vector<float>  *JetAk08Phi;
    vector<float>  *JetAk08E;
    vector<float>  *JetAk08Id;
    vector<float>  *JetAk08BDiscCisvV2;
    vector<float>  *JetAk08HadFlav;
    vector<float>  *JetAk08JecUncUp;
    vector<float>  *JetAk08JecUncDwn;

    TBranch        *b_JetAk08Pt;   
    TBranch        *b_JetAk08Eta;   
    TBranch        *b_JetAk08Phi;   
    TBranch        *b_JetAk08E;  
    TBranch        *b_JetAk08Id;   
    TBranch        *b_JetAk08BDiscCisvV2;   
    TBranch        *b_JetAk08HadFlav;  
    TBranch        *b_JetAk08JecUncUp;
    TBranch        *b_JetAk08JecUncDwn;
    

    // Other functions
    string   CreateOutputFileName(bool, bool, int, bool, int, int, bool, bool);
    Int_t    Cut(Long64_t entry);
    Int_t    GetEntry(Long64_t entry);
    Long64_t LoadTree(Long64_t entry);
    void     Init(bool hasRecoInfo, bool hasGenInfo);
    void     Loop(bool hasRecoInfo = 1, bool hasGenInfo = 0, int year = 2017, int doQCD = 0, bool doSSign = 0, bool doInvMassCut = 1 ,   int doBJets = 0, int doPUStudy = -10,bool doFlat = 0, bool useRoch = 0, bool doVarWidth = 1);
    Bool_t   Notify();
    void     Show(Long64_t entry = -1);
    
    // Baobab->Bonzai acceptance
    void getMcNorm();
    std::vector<Double_t> InEvtWeightSums_;
    std::vector<Double_t> EvtWeightSums_;
    std::vector<Double_t> skimAccep_;
    Int_t InEvtCount_;
    Int_t EvtCount_;

    bool nEvents_10000;
    string outputDirectory;
    string fileName; 
    int year;
    float lumiScale;
    float puScale;
    bool useTriggerCorrection;
    bool useEfficiencyCorrection; 
    bool isData;
    int systematics;
    int direction;
    float xsecfactor;
    int jetPtCutMin;
    int jetPtCutMax;
    int jetEtaCutMin, jetEtaCutMax ;
    int ZPtCutMin;
    int ZEtaCutMin;
    int ZEtaCutMax;
    int METcut;
    string leptonFlavor;

    ClassDef(ZJetsAndDPS,2)
};

#endif


