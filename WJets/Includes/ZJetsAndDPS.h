//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 22 08:07:46 2012 by ROOT version 5.32/01
// from TTree tree/tree
//////////////////////////////////////////////////////////

#ifndef ZJetsAndDPS_h
#define ZJetsAndDPS_h

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
#include <RooUnfoldResponse.h>
#include <TDatime.h>
#include <TMath.h>
#include <TRandom3.h>

// Header file for the classes stored in the TTree if any.
#include "functions.h"
#include "getFilesAndHistograms.h"
#include "standalone_LumiReWeighting.h"
#include "HistoSet.h"

using namespace std;

ClassImp(HistoSet)

class ZJetsAndDPS: public HistoSet {

  public :

    //TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    TChain          *fChain;   //!pointer to the analyzed TTree or TChain
    TChain          *fBonzaiHeaderChain;
    //TTree          *tree;
    Int_t           fCurrent; //!current Tree number in a TChain


    // --- Declaration of leaf types ---
    Int_t           EvtPuCntTruth;
    Int_t           EvtPuCntObs;
    Int_t           GNup;

    // hasRecoInfo --- 
    vector<double> *EvtWeights;
    Int_t           EvtVtxCnt;
    Int_t           EvtRunNum;
    Int_t           EvtNum;
    
    vector<float>        *MuPt;
    vector<float>        *MuEta;
    vector<float>        *MuPhi;
    vector<float>        *MuVtxZ;
    vector<float>        *MuE;
    vector<float>        *MuCh;
    vector<float>        *MuDxy;
    vector<unsigned int> *MuIdTight;
    ULong64_t            TrigHltMu;
    vector<float>        *MuPfIso;

    vector<float>  *JetAk04Pt;
    vector<float>  *JetAk04Eta;
    vector<float>  *JetAk04Phi;
    vector<float>  *JetAk04E;
    vector<float>  *JetAk04Id;
    vector<int>    *JetAk04PuId;
    vector<float>  *JetAk04PuMva;
    // vector<float>  *JetAk04BTagCsv;
    vector<float>  *JetAk04BDiscCisvV2;
    vector<float>  *JetAk04HadFlav;
    // vector<float>  *JetAk04PartFlav;

    vector<float>  *METPt;
    vector<float>  *METPx;
    vector<float>  *METPy;
    vector<float>  *METPz;
    vector<float>  *METE;
    vector<float>  *METPhi;
    vector<float>  *METsig;
    ULong64_t      TrigMET;
    ULong64_t      TrigMETBit;
    vector<bool>  *EvtFilterbadChCandidate;
    vector<bool>  *EvtFilterbadPFMuon;

    // hasGenInfo ---
    Double_t        mcEveWeight_;
    Double_t        mcSherpaSumWeight3_ ;
    vector<double> *mcSherpaWeights_;

    vector<float>  *GLepBarePt;
    vector<float>  *GLepBareEta;
    vector<float>  *GLepBarePhi;
    vector<float>  *GLepBareE;
    vector<int>    *GLepBareId;
    vector<int>    *GLepBareSt;
    vector<bool>   *GLepBarePrompt;
    vector<bool>   *GLepBareTauProd;
    
    vector<float>  *GPhotPt;
    vector<float>  *GPhotEta;
    vector<float>  *GPhotPhi;
    vector<float>  *GPhotE;
    vector<int>    *GPhotMotherId;
    vector<int>    *GPhotSt;
    
    vector<float>   *GLepClosePhotPt;
    vector<float>   *GLepClosePhotEta;
    vector<float>   *GLepClosePhotPhi;
    vector<float>   *GLepClosePhotE;
    vector<int>     *GLepClosePhotId;
    vector<int>     *GLepClosePhotMother0Id;
    vector<int>     *GLepClosePhotMotherCnt;
    vector<int>     *GLepClosePhotSt;

    vector<float>  *GJetAk04Pt;
    vector<float>  *GJetAk04Eta;
    vector<float>  *GJetAk04Phi;
    vector<float>  *GJetAk04E;

    // hasPartonInfo ---
    vector<double>  *dpsParton_Pt;
    vector<double>  *dpsParton_Eta;
    vector<double>  *dpsParton_Phi;
    vector<double>  *dpsParton_E;
    vector<int>     *genMatchDPSpar;
    vector<double>  *dpsParton_dR;


    // --- List of branches ---
    
    TBranch        *b_EvtPuCntTruth;   //!
    TBranch        *b_EvtPuCntObs;   //!
    TBranch        *b_GNup;   //!

    // hasRecoInfo ---
    TBranch        *b_EvtWeights;   //!
    TBranch        *b_EvtVtxCnt;   //!
    TBranch        *b_EvtRunNum;   //!
    TBranch        *b_EvtNum;   //!

    TBranch        *b_MuPt;   //!
    TBranch        *b_MuEta;   //!
    TBranch        *b_MuPhi;   //!
    TBranch        *b_MuVtxZ;   //!
    TBranch        *b_MuE;   //!
    TBranch        *b_MuCh;   //!
    TBranch        *b_MuDxy;   //!
    TBranch        *b_MuIdTight;   //!
    TBranch        *b_TrigHltMu;   //!
    TBranch        *b_MuPfIso;   //!

    TBranch        *b_JetAk04Pt;   //!
    TBranch        *b_JetAk04Eta;   //!
    TBranch        *b_JetAk04Phi;   //!
    TBranch        *b_JetAk04E;   //!
    TBranch        *b_JetAk04Id;   //!
    TBranch        *b_JetAk04PuId;   //!
    TBranch        *b_JetAk04PuMva;   //!
    TBranch        *b_JetAk04BDiscCisvV2;   //!
    TBranch        *b_JetAk04HadFlav;   //!
  
    TBranch        *b_METPt;   //!
    TBranch        *b_METPx;   //!
    TBranch        *b_METPy;   //!
    TBranch        *b_METPz;   //!
    TBranch        *b_METE;   //!
    TBranch        *b_TrigMET;   //!
    TBranch        *b_TrigMETBit;   //!
    TBranch        *b_EvtFilterbadChCandidate;   //!
    TBranch        *b_EvtFilterbadPFMuon;   //!
    TBranch        *b_METPhi;   //!
    TBranch        *b_METsig;   //!

    // hasGenInfo ---
    TBranch        *b_mcEveWeight_;   //!
    TBranch        *b_mcSherpaSumWeight3_;   //!
    TBranch        *b_mcSherpaWeights_;   //!

    TBranch        *b_GLepBarePt;   //!
    TBranch        *b_GLepBareEta;   //!
    TBranch        *b_GLepBarePhi;   //!
    TBranch        *b_GLepBareE;   //!
    TBranch        *b_genLepQ_;   //!
    TBranch        *b_GLepBareId;   //!
    TBranch        *b_GLepBareSt;   //!
    TBranch        *b_GLepBarePrompt;   //!
    TBranch        *b_GLepBareTauProd;   //!
    
    TBranch        *b_GPhotPt;   //!
    TBranch        *b_GPhotEta;   //!
    TBranch        *b_GPhotPhi;   //!
    TBranch        *b_GPhotE;   //!
    TBranch        *b_GPhotMotherId;   //!
    TBranch        *b_GPhotSt;   //!
    
    TBranch        *b_GLepClosePhotPt;   //!
    TBranch        *b_GLepClosePhotEta;   //!
    TBranch        *b_GLepClosePhotPhi;   //!
    TBranch        *b_GLepClosePhotE;   //!
    TBranch        *b_GLepClosePhotId;   //!
    TBranch        *b_GLepClosePhotMother0Id;   //!
    TBranch        *b_GLepClosePhotMotherCnt;   //!
    TBranch        *b_GLepClosePhotSt;   //!
    
    TBranch        *b_GJetAk04Pt;   //!
    TBranch        *b_GJetAk04Eta;   //!
    TBranch        *b_GJetAk04Phi;   //!
    TBranch        *b_GJetAk04E;   //!

    // hasPartonInfo ---
    TBranch        *b_dpsParton_Pt;   //!
    TBranch        *b_dpsParton_Eta;   //!
    TBranch        *b_dpsParton_Phi;   //!
    TBranch        *b_dpsParton_E;   //!
    TBranch        *b_genMatchDPSpar;   //!
    TBranch        *b_dpsParton_dR;   //!

    // ...
    TBranch        *b_gsfElecPt_;   //!
    TBranch        *b_gsfElecEta_;   //!
    TBranch        *b_gsfElecPhi_;   //!
    TBranch        *b_gsfElecEnergy_;   //!
    
    // Constructor and destructor
    ZJetsAndDPS(string fileName_, int year_ = 2017, float lumiScale_ = 1., float puScale_ = 1., bool useTriggerCorrection_ = 0, bool useEfficiencyCorrection_ = 0, int systematics_ = 0, int direction_ = 0, float xsecfactor_ = 1., int jetPtCutMin_ = 20, int jetPtCutMax_ = 0, int ZPtCutMin_ = 0 , int ZEtaCutMin_ = -999999, int ZEtaCutMax_ = 999999, int METcut_ = -30, bool nEvents_10000_ = 0, int jetEtaCutMin_ = -24, int jetEtaCutMax_ = 24) ; 
    ~ZJetsAndDPS();

    // Other functions
    string   CreateOutputFileName(bool, bool, int, bool, int, int, bool, bool , string pdfSet = "", int aa = -1);
    Int_t    Cut(Long64_t entry);
    Int_t    GetEntry(Long64_t entry);
    Long64_t LoadTree(Long64_t entry);
    void     Init(bool hasRecoInfo, bool hasGenInfo, bool hasPartonInfo);
    void     Loop(bool hasRecoInfo = 1, bool hasGenInfo = 0, int year = 2017, int doQCD = 0, bool doSSign = 0, bool doInvMassCut = 1 ,   int doBJets = 0, int doPUStudy = -10,bool doFlat = 0, bool useRoch = 0, bool doVarWidth = 1, bool hasPartonInfo = 0, string pdfSet = "", int pdfMember = 0);
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

    TH1D *FlatNVtxWeight;

    ClassDef(ZJetsAndDPS,2)
};
#endif


