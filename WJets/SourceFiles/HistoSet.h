
#ifndef _HistoSet_H_
#define _HistoSet_H_

#include <iostream>
#include <TObject.h>
#include <TH1.h>
#include <TH2.h>
#include "TAxis.h" 
#include <TArray.h>
#include <vector>
#include <cstdarg>

using namespace std;

class HistoSet: public TObject{
public:

    HistoSet(string leptonFlavor = "DMu");
    ~HistoSet();

    vector<TH1*> listOfHistograms;

    vector<double> makeVector(int num, ...);
    void insertVector(vector<double>& veca, int num, ...);
    vector<double> buildVecFineBin(int nStdBin, double arrStdBin[], int factChop);
    TH1D* newTH1D(string, string, string, int, double*);
    TH1D* newTH1D(string, string, string, int, double, double);
    TH1D* newTH1D(string, string, string, vector<double>&);
    TH2D* newTH2D(string, string, int, double*, int, double*);
    TH2D* newTH2D(string, string, vector<double>&, vector<double>&);
    TH2D* newTH2D(string, string, int, double*, int, double, double);
    TH2D* newTH2D(string, string, int, double, double, int, double*);
    TH2D* newTH2D(string, string, int, double, double, int, double, double);
    int NbinsEta2Dtest; 
    static const int NbinsEta2D = 9 ;
    double j_pT_range[NbinsEta2D];
    double j_Y_range[NbinsEta2D];

    //***************************** Basic plots for Wjets *****************************//

    TH1D *NEventsPassCuts;
    
    //--- For calculating b-tagging efficiencies---
    TH2D *h_pt_eta_b;
    TH2D *h_pt_eta_b_tagged;
    TH2D *h_pt_eta_c;
    TH2D *h_pt_eta_c_tagged;
    TH2D *h_pt_eta_udsg;
    TH2D *h_pt_eta_udsg_tagged;
    
    TH1D *h_pt_b;
    TH1D *h_pt_b_tagged;
    TH1D *h_pt_c;
    TH1D *h_pt_c_tagged;
    TH1D *h_pt_udsg;
    TH1D *h_pt_udsg_tagged;

    TH1D *NumRecoVtx;
    TH1D *NumPUTruthVtx;
    TH1D *NumPUObsVtx;
    
    //--- Jet multiplicity
    TH1D *ZNGoodJets_Zexc;
    TH1D *genZNGoodJets_Zexc;
    TH2D *hresponseZNGoodJets_Zexc;
    
    TH1D *ZNGoodJets_Zinc;
    TH1D *genZNGoodJets_Zinc;
    TH2D *hresponseZNGoodJets_Zinc;
    
    TH1D *ZNGoodJetsFull_Zexc;
    TH1D *genZNGoodJetsFull_Zexc;
    TH2D *hresponseZNGoodJetsFull_Zexc;
    TH1D *ZNGoodJetsFull_Zinc;
    TH1D *genZNGoodJetsFull_Zinc;
    TH2D *hresponseZNGoodJetsFull_Zinc;

    //--- Jet pt
    TH1D *FirstJetPt_Zinc1jet;
    TH1D *FirstJetPt_1_Zinc1jet;
    TH1D *FirstJetPt_2_Zinc1jet;
    TH1D *SecondJetPt_Zinc2jet;
    TH1D *SecondJetPt_1_Zinc2jet;
    TH1D *SecondJetPt_2_Zinc2jet;
    TH1D *ThirdJetPt_Zinc3jet;
    TH1D *ThirdJetPt_1_Zinc3jet;
    TH1D *ThirdJetPt_2_Zinc3jet;
    TH1D *FourthJetPt_Zinc4jet;
    TH1D *FourthJetPt_1_Zinc4jet;
    TH1D *FourthJetPt_2_Zinc4jet;
    TH1D *FifthJetPt_Zinc5jet;
    TH1D *FifthJetPt_1_Zinc5jet;
    TH1D *FifthJetPt_2_Zinc5jet;
    TH1D *SixthJetPt_Zinc6jet;
    TH1D *SixthJetPt_1_Zinc6jet;
    TH1D *SixthJetPt_2_Zinc6jet;

    TH1D *genFirstJetPt_Zinc1jet;
    TH1D *genFirstJetPt_1_Zinc1jet;
    TH1D *genFirstJetPt_2_Zinc1jet;
    TH1D *genSecondJetPt_Zinc2jet;
    TH1D *genSecondJetPt_1_Zinc2jet;
    TH1D *genSecondJetPt_2_Zinc2jet;
    TH1D *genThirdJetPt_Zinc3jet;
    TH1D *genThirdJetPt_1_Zinc3jet;
    TH1D *genThirdJetPt_2_Zinc3jet;
    TH1D *genFourthJetPt_Zinc4jet;
    TH1D *genFourthJetPt_1_Zinc4jet;
    TH1D *genFourthJetPt_2_Zinc4jet;
    TH1D *genFifthJetPt_Zinc5jet;
    TH1D *genFifthJetPt_1_Zinc5jet;
    TH1D *genFifthJetPt_2_Zinc5jet;
    TH1D *genSixthJetPt_Zinc6jet;
    TH1D *genSixthJetPt_1_Zinc6jet;
    TH1D *genSixthJetPt_2_Zinc6jet;

    TH2D *hresponseFirstJetPt_Zinc1jet;
    TH2D *hresponseFirstJetPt_1_Zinc1jet;
    TH2D *hresponseFirstJetPt_2_Zinc1jet;
    TH2D *hresponseSecondJetPt_Zinc2jet;
    TH2D *hresponseSecondJetPt_1_Zinc2jet;
    TH2D *hresponseSecondJetPt_2_Zinc2jet;
    TH2D *hresponseThirdJetPt_Zinc3jet;
    TH2D *hresponseThirdJetPt_1_Zinc3jet;
    TH2D *hresponseThirdJetPt_2_Zinc3jet;
    TH2D *hresponseFourthJetPt_Zinc4jet;
    TH2D *hresponseFourthJetPt_1_Zinc4jet;
    TH2D *hresponseFourthJetPt_2_Zinc4jet;
    TH2D *hresponseFifthJetPt_Zinc5jet;
    TH2D *hresponseFifthJetPt_1_Zinc5jet;
    TH2D *hresponseFifthJetPt_2_Zinc5jet;
    TH2D *hresponseSixthJetPt_Zinc6jet;
    TH2D *hresponseSixthJetPt_1_Zinc6jet;
    TH2D *hresponseSixthJetPt_2_Zinc6jet;

    ///////////////////////////////////////////////////////////////////////////

    ///////// AK4 Jet Distributions

    //leading jet pt -- alpha-s
    TH1D *LeadingJetPt_Zinc1jet;
    TH1D *LeadingJetPt_Zinc2jet;
    TH1D *LeadingJetPt_Zinc3jet;
    TH1D *LeadingJetPt_Zinc4jet;
    // TH1D *LeadingJetPt_1_Zinc1jet;
    // TH1D *LeadingJetPt_1_Zinc2jet;
    // TH1D *LeadingJetPt_1_Zinc3jet;
    // TH1D *LeadingJetPt_1_Zinc4jet;
    TH1D *LeadingJetPt_2_Zinc1jet;
    TH1D *LeadingJetPt_2_Zinc2jet;
    TH1D *LeadingJetPt_2_Zinc3jet;
    TH1D *LeadingJetPt_2_Zinc4jet;

    TH1D *genLeadingJetPt_Zinc1jet;
    TH1D *genLeadingJetPt_Zinc2jet;
    TH1D *genLeadingJetPt_Zinc3jet;
    TH1D *genLeadingJetPt_Zinc4jet;
    // TH1D *genLeadingJetPt_1_Zinc1jet;
    // TH1D *genLeadingJetPt_1_Zinc2jet;
    // TH1D *genLeadingJetPt_1_Zinc3jet;
    // TH1D *genLeadingJetPt_1_Zinc4jet;
    TH1D *genLeadingJetPt_2_Zinc1jet;
    TH1D *genLeadingJetPt_2_Zinc2jet;
    TH1D *genLeadingJetPt_2_Zinc3jet;
    TH1D *genLeadingJetPt_2_Zinc4jet;

    TH2D *hresponseLeadingJetPt_Zinc1jet;
    TH2D *hresponseLeadingJetPt_Zinc2jet;
    TH2D *hresponseLeadingJetPt_Zinc3jet;
    TH2D *hresponseLeadingJetPt_Zinc4jet;
    // TH2D *hresponseLeadingJetPt_1_Zinc1jet;
    // TH2D *hresponseLeadingJetPt_1_Zinc2jet;
    // TH2D *hresponseLeadingJetPt_1_Zinc3jet;
    // TH2D *hresponseLeadingJetPt_1_Zinc4jet;
    TH2D *hresponseLeadingJetPt_2_Zinc1jet;
    TH2D *hresponseLeadingJetPt_2_Zinc2jet;
    TH2D *hresponseLeadingJetPt_2_Zinc3jet;
    TH2D *hresponseLeadingJetPt_2_Zinc4jet;

    // exclusive jet requirement
    TH1D *LeadingJetPt_Zexc1jet;
    TH1D *LeadingJetPt_Zexc2jet;
    TH1D *LeadingJetPt_Zexc3jet;
    TH1D *LeadingJetPt_2_Zexc1jet;
    TH1D *LeadingJetPt_2_Zexc2jet;
    TH1D *LeadingJetPt_2_Zexc3jet;

    TH1D *genLeadingJetPt_Zexc1jet;
    TH1D *genLeadingJetPt_Zexc2jet;
    TH1D *genLeadingJetPt_Zexc3jet;
    TH1D *genLeadingJetPt_2_Zexc1jet;
    TH1D *genLeadingJetPt_2_Zexc2jet;
    TH1D *genLeadingJetPt_2_Zexc3jet;

    TH2D *hresponseLeadingJetPt_Zexc1jet;
    TH2D *hresponseLeadingJetPt_Zexc2jet;
    TH2D *hresponseLeadingJetPt_Zexc3jet;
    TH2D *hresponseLeadingJetPt_2_Zexc1jet;
    TH2D *hresponseLeadingJetPt_2_Zexc2jet;
    TH2D *hresponseLeadingJetPt_2_Zexc3jet;

    // looking to study particle-level corrections of ratios
    TH1D *LeadingJetPt_MIGRATIONS_Zinc1jet;
    TH1D *LeadingJetPt_MIGRATIONS_Zinc2jet;
    TH1D *LeadingJetPt_MIGRATIONS_Zinc3jet;
    TH1D *genLeadingJetPt_MIGRATIONS_Zinc1jet;
    TH1D *genLeadingJetPt_MIGRATIONS_Zinc2jet;
    TH1D *genLeadingJetPt_MIGRATIONS_Zinc3jet;

    //--- Jet HT/2 -- alpha-s
    TH1D *HTover2_Zinc2jet;
    TH1D *HTover2_Zinc3jet;
    TH1D *HTover2_Zinc4jet;
    // TH1D *HTover2_1_Zinc2jet;
    // TH1D *HTover2_1_Zinc3jet;
    // TH1D *HTover2_1_Zinc4jet;
    TH1D *HTover2_2_Zinc2jet;
    TH1D *HTover2_2_Zinc3jet;
    TH1D *HTover2_2_Zinc4jet;

    TH1D *genHTover2_Zinc2jet;
    TH1D *genHTover2_Zinc3jet;
    TH1D *genHTover2_Zinc4jet;
    // TH1D *genHTover2_1_Zinc2jet;
    // TH1D *genHTover2_1_Zinc3jet;
    // TH1D *genHTover2_1_Zinc4jet;
    TH1D *genHTover2_2_Zinc2jet;
    TH1D *genHTover2_2_Zinc3jet;
    TH1D *genHTover2_2_Zinc4jet;
    
    TH2D *hresponseHTover2_Zinc2jet;
    TH2D *hresponseHTover2_Zinc3jet;
    TH2D *hresponseHTover2_Zinc4jet;
    // TH2D *hresponseHTover2_1_Zinc2jet;
    // TH2D *hresponseHTover2_1_Zinc3jet;
    // TH2D *hresponseHTover2_1_Zinc4jet;
    TH2D *hresponseHTover2_2_Zinc2jet;
    TH2D *hresponseHTover2_2_Zinc3jet;
    TH2D *hresponseHTover2_2_Zinc4jet;

    //Lepton Pt + LJ Pt -- alpha-s
    TH1D *LepPtPlusLeadingJetPt_Zinc1jet;
    TH1D *LepPtPlusLeadingJetPt_Zinc2jet;
    TH1D *LepPtPlusLeadingJetPt_Zinc3jet;
    TH1D *LepPtPlusLeadingJetPt_Zinc4jet;
    // TH1D *LepPtPlusLeadingJetPt_1_Zinc1jet;
    // TH1D *LepPtPlusLeadingJetPt_1_Zinc2jet;
    // TH1D *LepPtPlusLeadingJetPt_1_Zinc3jet;
    // TH1D *LepPtPlusLeadingJetPt_1_Zinc4jet;
    TH1D *LepPtPlusLeadingJetPt_2_Zinc1jet;
    TH1D *LepPtPlusLeadingJetPt_2_Zinc2jet;
    TH1D *LepPtPlusLeadingJetPt_2_Zinc3jet;
    TH1D *LepPtPlusLeadingJetPt_2_Zinc4jet;

    TH1D *genLepPtPlusLeadingJetPt_Zinc1jet;
    TH1D *genLepPtPlusLeadingJetPt_Zinc2jet;
    TH1D *genLepPtPlusLeadingJetPt_Zinc3jet;
    TH1D *genLepPtPlusLeadingJetPt_Zinc4jet;
    // TH1D *genLepPtPlusLeadingJetPt_1_Zinc1jet;
    // TH1D *genLepPtPlusLeadingJetPt_1_Zinc2jet;
    // TH1D *genLepPtPlusLeadingJetPt_1_Zinc3jet;
    // TH1D *genLepPtPlusLeadingJetPt_1_Zinc4jet;
    TH1D *genLepPtPlusLeadingJetPt_2_Zinc1jet;
    TH1D *genLepPtPlusLeadingJetPt_2_Zinc2jet;
    TH1D *genLepPtPlusLeadingJetPt_2_Zinc3jet;
    TH1D *genLepPtPlusLeadingJetPt_2_Zinc4jet;

    TH2D *hresponseLepPtPlusLeadingJetPt_Zinc1jet;
    TH2D *hresponseLepPtPlusLeadingJetPt_Zinc2jet;
    TH2D *hresponseLepPtPlusLeadingJetPt_Zinc3jet;
    TH2D *hresponseLepPtPlusLeadingJetPt_Zinc4jet;
    // TH2D *hresponseLepPtPlusLeadingJetPt_1_Zinc1jet;
    // TH2D *hresponseLepPtPlusLeadingJetPt_1_Zinc2jet;
    // TH2D *hresponseLepPtPlusLeadingJetPt_1_Zinc3jet;
    // TH2D *hresponseLepPtPlusLeadingJetPt_1_Zinc4jet;
    TH2D *hresponseLepPtPlusLeadingJetPt_2_Zinc1jet;
    TH2D *hresponseLepPtPlusLeadingJetPt_2_Zinc2jet;
    TH2D *hresponseLepPtPlusLeadingJetPt_2_Zinc3jet;
    TH2D *hresponseLepPtPlusLeadingJetPt_2_Zinc4jet;

    // looking to study particle-level corrections of ratios
    TH1D *LepPtPlusLeadingJetPt_MIGRATIONS_Zinc1jet;
    TH1D *LepPtPlusLeadingJetPt_MIGRATIONS_Zinc2jet;
    TH1D *LepPtPlusLeadingJetPt_MIGRATIONS_Zinc3jet;
    TH1D *genLepPtPlusLeadingJetPt_MIGRATIONS_Zinc1jet;
    TH1D *genLepPtPlusLeadingJetPt_MIGRATIONS_Zinc2jet;
    TH1D *genLepPtPlusLeadingJetPt_MIGRATIONS_Zinc3jet;

    //exclusive jet requirement
    TH1D *LepPtPlusLeadingJetPt_Zexc1jet;
    TH1D *LepPtPlusLeadingJetPt_Zexc2jet;
    TH1D *LepPtPlusLeadingJetPt_Zexc3jet;
    TH1D *LepPtPlusLeadingJetPt_2_Zexc1jet;
    TH1D *LepPtPlusLeadingJetPt_2_Zexc2jet;
    TH1D *LepPtPlusLeadingJetPt_2_Zexc3jet;

    TH1D *genLepPtPlusLeadingJetPt_Zexc1jet;
    TH1D *genLepPtPlusLeadingJetPt_Zexc2jet;
    TH1D *genLepPtPlusLeadingJetPt_Zexc3jet;
    TH1D *genLepPtPlusLeadingJetPt_2_Zexc1jet;
    TH1D *genLepPtPlusLeadingJetPt_2_Zexc2jet;
    TH1D *genLepPtPlusLeadingJetPt_2_Zexc3jet;

    TH2D *hresponseLepPtPlusLeadingJetPt_Zexc1jet;
    TH2D *hresponseLepPtPlusLeadingJetPt_Zexc2jet;
    TH2D *hresponseLepPtPlusLeadingJetPt_Zexc3jet;
    TH2D *hresponseLepPtPlusLeadingJetPt_2_Zexc1jet;
    TH2D *hresponseLepPtPlusLeadingJetPt_2_Zexc2jet;
    TH2D *hresponseLepPtPlusLeadingJetPt_2_Zexc3jet;

    //Lepton Pt + HT/2 -- alpha-s
    TH1D *LepPtPlusHTover2_Zinc2jet;
    TH1D *LepPtPlusHTover2_Zinc3jet;
    TH1D *LepPtPlusHTover2_Zinc4jet;
    // TH1D *LepPtPlusHTover2_1_Zinc2jet;
    // TH1D *LepPtPlusHTover2_1_Zinc3jet;
    // TH1D *LepPtPlusHTover2_1_Zinc4jet;
    TH1D *LepPtPlusHTover2_2_Zinc2jet;
    TH1D *LepPtPlusHTover2_2_Zinc3jet;
    TH1D *LepPtPlusHTover2_2_Zinc4jet;

    TH1D *genLepPtPlusHTover2_Zinc2jet;
    TH1D *genLepPtPlusHTover2_Zinc3jet;
    TH1D *genLepPtPlusHTover2_Zinc4jet;
    // TH1D *genLepPtPlusHTover2_1_Zinc2jet;
    // TH1D *genLepPtPlusHTover2_1_Zinc3jet;
    // TH1D *genLepPtPlusHTover2_1_Zinc4jet;
    TH1D *genLepPtPlusHTover2_2_Zinc2jet;
    TH1D *genLepPtPlusHTover2_2_Zinc3jet;
    TH1D *genLepPtPlusHTover2_2_Zinc4jet;

    TH2D *hresponseLepPtPlusHTover2_Zinc2jet;
    TH2D *hresponseLepPtPlusHTover2_Zinc3jet;
    TH2D *hresponseLepPtPlusHTover2_Zinc4jet;
    // TH2D *hresponseLepPtPlusHTover2_1_Zinc2jet;
    // TH2D *hresponseLepPtPlusHTover2_1_Zinc3jet;
    // TH2D *hresponseLepPtPlusHTover2_1_Zinc4jet;
    TH2D *hresponseLepPtPlusHTover2_2_Zinc2jet;
    TH2D *hresponseLepPtPlusHTover2_2_Zinc3jet;
    TH2D *hresponseLepPtPlusHTover2_2_Zinc4jet;

    // WpT -- alpha-s
    TH1D *ZPt_Zinc1jet;
    TH1D *ZPt_Zinc2jet;
    TH1D *ZPt_Zinc3jet;
    TH1D *ZPt_Zinc4jet;
    // TH1D *ZPt_1_Zinc1jet;
    // TH1D *ZPt_1_Zinc2jet;
    // TH1D *ZPt_1_Zinc3jet;
    // TH1D *ZPt_1_Zinc4jet;
    TH1D *ZPt_2_Zinc1jet;
    TH1D *ZPt_2_Zinc2jet;
    TH1D *ZPt_2_Zinc3jet;
    TH1D *ZPt_2_Zinc4jet;

    TH1D *genZPt_Zinc1jet;
    TH1D *genZPt_Zinc2jet;
    TH1D *genZPt_Zinc3jet;
    TH1D *genZPt_Zinc4jet;
    // TH1D *genZPt_1_Zinc1jet;
    // TH1D *genZPt_1_Zinc2jet;
    // TH1D *genZPt_1_Zinc3jet;
    // TH1D *genZPt_1_Zinc4jet;
    TH1D *genZPt_2_Zinc1jet;
    TH1D *genZPt_2_Zinc2jet;
    TH1D *genZPt_2_Zinc3jet;
    TH1D *genZPt_2_Zinc4jet;

    TH2D *hresponseZPt_Zinc1jet;
    TH2D *hresponseZPt_Zinc2jet;
    TH2D *hresponseZPt_Zinc3jet;
    TH2D *hresponseZPt_Zinc4jet;
    // TH2D *hresponseZPt_1_Zinc1jet;
    // TH2D *hresponseZPt_1_Zinc2jet;
    // TH2D *hresponseZPt_1_Zinc3jet;
    // TH2D *hresponseZPt_1_Zinc4jet;
    TH2D *hresponseZPt_2_Zinc1jet;
    TH2D *hresponseZPt_2_Zinc2jet;
    TH2D *hresponseZPt_2_Zinc3jet;
    TH2D *hresponseZPt_2_Zinc4jet;

    //WpT + LJ pT -- alpha-s
    TH1D *ZPtPlusLeadingJetPt_Zinc1jet;
    TH1D *ZPtPlusLeadingJetPt_Zinc2jet;
    TH1D *ZPtPlusLeadingJetPt_Zinc3jet;
    TH1D *ZPtPlusLeadingJetPt_Zinc4jet;
    // TH1D *ZPtPlusLeadingJetPt_1_Zinc1jet;
    // TH1D *ZPtPlusLeadingJetPt_1_Zinc2jet;
    // TH1D *ZPtPlusLeadingJetPt_1_Zinc3jet;
    // TH1D *ZPtPlusLeadingJetPt_1_Zinc4jet;
    TH1D *ZPtPlusLeadingJetPt_2_Zinc1jet;
    TH1D *ZPtPlusLeadingJetPt_2_Zinc2jet;
    TH1D *ZPtPlusLeadingJetPt_2_Zinc3jet;
    TH1D *ZPtPlusLeadingJetPt_2_Zinc4jet;

    TH1D *genZPtPlusLeadingJetPt_Zinc1jet;
    TH1D *genZPtPlusLeadingJetPt_Zinc2jet;
    TH1D *genZPtPlusLeadingJetPt_Zinc3jet;
    TH1D *genZPtPlusLeadingJetPt_Zinc4jet;
    // TH1D *genZPtPlusLeadingJetPt_1_Zinc1jet;
    // TH1D *genZPtPlusLeadingJetPt_1_Zinc2jet;
    // TH1D *genZPtPlusLeadingJetPt_1_Zinc3jet;
    // TH1D *genZPtPlusLeadingJetPt_1_Zinc4jet;
    TH1D *genZPtPlusLeadingJetPt_2_Zinc1jet;
    TH1D *genZPtPlusLeadingJetPt_2_Zinc2jet;
    TH1D *genZPtPlusLeadingJetPt_2_Zinc3jet;
    TH1D *genZPtPlusLeadingJetPt_2_Zinc4jet;

    TH2D *hresponseZPtPlusLeadingJetPt_Zinc1jet;
    TH2D *hresponseZPtPlusLeadingJetPt_Zinc2jet;
    TH2D *hresponseZPtPlusLeadingJetPt_Zinc3jet;
    TH2D *hresponseZPtPlusLeadingJetPt_Zinc4jet;
    // TH2D *hresponseZPtPlusLeadingJetPt_1_Zinc1jet;
    // TH2D *hresponseZPtPlusLeadingJetPt_1_Zinc2jet;
    // TH2D *hresponseZPtPlusLeadingJetPt_1_Zinc3jet;
    // TH2D *hresponseZPtPlusLeadingJetPt_1_Zinc4jet;
    TH2D *hresponseZPtPlusLeadingJetPt_2_Zinc1jet;
    TH2D *hresponseZPtPlusLeadingJetPt_2_Zinc2jet;
    TH2D *hresponseZPtPlusLeadingJetPt_2_Zinc3jet;
    TH2D *hresponseZPtPlusLeadingJetPt_2_Zinc4jet;

    //WpT + HT,2/2 -- alpha-s
    TH1D *ZPtPlusHTover2_Zinc2jet;
    TH1D *ZPtPlusHTover2_Zinc3jet;
    TH1D *ZPtPlusHTover2_Zinc4jet;
    // TH1D *ZPtPlusHTover2_1_Zinc2jet;
    // TH1D *ZPtPlusHTover2_1_Zinc3jet;
    // TH1D *ZPtPlusHTover2_1_Zinc4jet;
    TH1D *ZPtPlusHTover2_2_Zinc2jet;
    TH1D *ZPtPlusHTover2_2_Zinc3jet;
    TH1D *ZPtPlusHTover2_2_Zinc4jet;

    TH1D *genZPtPlusHTover2_Zinc2jet;
    TH1D *genZPtPlusHTover2_Zinc3jet;
    TH1D *genZPtPlusHTover2_Zinc4jet;
    // TH1D *genZPtPlusHTover2_1_Zinc2jet;
    // TH1D *genZPtPlusHTover2_1_Zinc3jet;
    // TH1D *genZPtPlusHTover2_1_Zinc4jet;
    TH1D *genZPtPlusHTover2_2_Zinc2jet;
    TH1D *genZPtPlusHTover2_2_Zinc3jet;
    TH1D *genZPtPlusHTover2_2_Zinc4jet;

    TH2D *hresponseZPtPlusHTover2_Zinc2jet;
    TH2D *hresponseZPtPlusHTover2_Zinc3jet;
    TH2D *hresponseZPtPlusHTover2_Zinc4jet;
    // TH2D *hresponseZPtPlusHTover2_1_Zinc2jet;
    // TH2D *hresponseZPtPlusHTover2_1_Zinc3jet;
    // TH2D *hresponseZPtPlusHTover2_1_Zinc4jet;
    TH2D *hresponseZPtPlusHTover2_2_Zinc2jet;
    TH2D *hresponseZPtPlusHTover2_2_Zinc3jet;
    TH2D *hresponseZPtPlusHTover2_2_Zinc4jet;


    ///////// AK8 Jet Distributions

    // Leading AK8 Jet Pt -- alpha-s
    TH1D *LeadingJetAK8Pt_Zinc1jet;
    TH1D *LeadingJetAK8Pt_Zinc2jet;
    TH1D *LeadingJetAK8Pt_Zinc3jet;
    TH1D *LeadingJetAK8Pt_2_Zinc1jet;
    TH1D *LeadingJetAK8Pt_2_Zinc2jet;
    TH1D *LeadingJetAK8Pt_2_Zinc3jet;

    TH1D *genLeadingJetAK8Pt_Zinc1jet;
    TH1D *genLeadingJetAK8Pt_Zinc2jet;
    TH1D *genLeadingJetAK8Pt_Zinc3jet;
    TH1D *genLeadingJetAK8Pt_2_Zinc1jet;
    TH1D *genLeadingJetAK8Pt_2_Zinc2jet;
    TH1D *genLeadingJetAK8Pt_2_Zinc3jet;

    TH2D *hresponseLeadingJetAK8Pt_Zinc1jet;
    TH2D *hresponseLeadingJetAK8Pt_Zinc2jet;
    TH2D *hresponseLeadingJetAK8Pt_Zinc3jet;
    TH2D *hresponseLeadingJetAK8Pt_2_Zinc1jet;
    TH2D *hresponseLeadingJetAK8Pt_2_Zinc2jet;
    TH2D *hresponseLeadingJetAK8Pt_2_Zinc3jet;

    // exclusive jet requirement
    TH1D *LeadingJetAK8Pt_Zexc1jet;
    TH1D *LeadingJetAK8Pt_Zexc2jet;
    TH1D *LeadingJetAK8Pt_Zexc3jet;
    TH1D *LeadingJetAK8Pt_2_Zexc1jet;
    TH1D *LeadingJetAK8Pt_2_Zexc2jet;
    TH1D *LeadingJetAK8Pt_2_Zexc3jet;

    TH1D *genLeadingJetAK8Pt_Zexc1jet;
    TH1D *genLeadingJetAK8Pt_Zexc2jet;
    TH1D *genLeadingJetAK8Pt_Zexc3jet;
    TH1D *genLeadingJetAK8Pt_2_Zexc1jet;
    TH1D *genLeadingJetAK8Pt_2_Zexc2jet;
    TH1D *genLeadingJetAK8Pt_2_Zexc3jet;

    TH2D *hresponseLeadingJetAK8Pt_Zexc1jet;
    TH2D *hresponseLeadingJetAK8Pt_Zexc2jet;
    TH2D *hresponseLeadingJetAK8Pt_Zexc3jet;
    TH2D *hresponseLeadingJetAK8Pt_2_Zexc1jet;
    TH2D *hresponseLeadingJetAK8Pt_2_Zexc2jet;
    TH2D *hresponseLeadingJetAK8Pt_2_Zexc3jet;

    // Lepton Pt + Leading AK8 Jet Pt -- alpha-s
    TH1D *LepPtPlusLeadingJetAK8Pt_Zinc1jet;
    TH1D *LepPtPlusLeadingJetAK8Pt_Zinc2jet;
    TH1D *LepPtPlusLeadingJetAK8Pt_Zinc3jet;
    TH1D *LepPtPlusLeadingJetAK8Pt_2_Zinc1jet;
    TH1D *LepPtPlusLeadingJetAK8Pt_2_Zinc2jet;
    TH1D *LepPtPlusLeadingJetAK8Pt_2_Zinc3jet;

    TH1D *genLepPtPlusLeadingJetAK8Pt_Zinc1jet;
    TH1D *genLepPtPlusLeadingJetAK8Pt_Zinc2jet;
    TH1D *genLepPtPlusLeadingJetAK8Pt_Zinc3jet;
    TH1D *genLepPtPlusLeadingJetAK8Pt_2_Zinc1jet;
    TH1D *genLepPtPlusLeadingJetAK8Pt_2_Zinc2jet;
    TH1D *genLepPtPlusLeadingJetAK8Pt_2_Zinc3jet;

    TH2D *hresponseLepPtPlusLeadingJetAK8Pt_Zinc1jet;
    TH2D *hresponseLepPtPlusLeadingJetAK8Pt_Zinc2jet;
    TH2D *hresponseLepPtPlusLeadingJetAK8Pt_Zinc3jet;
    TH2D *hresponseLepPtPlusLeadingJetAK8Pt_2_Zinc1jet;
    TH2D *hresponseLepPtPlusLeadingJetAK8Pt_2_Zinc2jet;
    TH2D *hresponseLepPtPlusLeadingJetAK8Pt_2_Zinc3jet;

    // exclusive jet requirement
    TH1D *LepPtPlusLeadingJetAK8Pt_Zexc1jet;
    TH1D *LepPtPlusLeadingJetAK8Pt_Zexc2jet;
    TH1D *LepPtPlusLeadingJetAK8Pt_Zexc3jet;
    TH1D *LepPtPlusLeadingJetAK8Pt_2_Zexc1jet;
    TH1D *LepPtPlusLeadingJetAK8Pt_2_Zexc2jet;
    TH1D *LepPtPlusLeadingJetAK8Pt_2_Zexc3jet;

    TH1D *genLepPtPlusLeadingJetAK8Pt_Zexc1jet;
    TH1D *genLepPtPlusLeadingJetAK8Pt_Zexc2jet;
    TH1D *genLepPtPlusLeadingJetAK8Pt_Zexc3jet;
    TH1D *genLepPtPlusLeadingJetAK8Pt_2_Zexc1jet;
    TH1D *genLepPtPlusLeadingJetAK8Pt_2_Zexc2jet;
    TH1D *genLepPtPlusLeadingJetAK8Pt_2_Zexc3jet;

    TH2D *hresponseLepPtPlusLeadingJetAK8Pt_Zexc1jet;
    TH2D *hresponseLepPtPlusLeadingJetAK8Pt_Zexc2jet;
    TH2D *hresponseLepPtPlusLeadingJetAK8Pt_Zexc3jet;
    TH2D *hresponseLepPtPlusLeadingJetAK8Pt_2_Zexc1jet;
    TH2D *hresponseLepPtPlusLeadingJetAK8Pt_2_Zexc2jet;
    TH2D *hresponseLepPtPlusLeadingJetAK8Pt_2_Zexc3jet;

    ///////////////////////////////////////////////////////////////////////////

    //--- Jet Ht
    TH1D *JetsHT_Zinc1jet;
    TH1D *JetsHT_Zinc2jet;
    TH1D *JetsHT_Zinc3jet;
    TH1D *JetsHT_20_Zinc1jet;
    TH1D *JetsHT_20_Zinc2jet;
    TH1D *JetsHT_20_Zinc3jet; 
    TH1D *JetsHT_20_30_Zinc1jet;
    TH1D *JetsHT_20_30_Zinc2jet;
    TH1D *JetsHT_20_30_Zinc3jet; 
    TH1D *JetsHT_Zinc4jet;
    TH1D *JetsHT_Zinc5jet;
    TH1D *JetsHT_Zinc6jet;
    TH1D *JetsHT_1_Zinc1jet;
    TH1D *JetsHT_1_Zinc2jet;
    TH1D *JetsHT_1_Zinc3jet;
    TH1D *JetsHT_1_Zinc4jet;
    TH1D *JetsHT_1_Zinc5jet;
    TH1D *JetsHT_1_Zinc6jet;
    TH1D *JetsHT_2_Zinc1jet;
    TH1D *JetsHT_2_Zinc2jet;
    TH1D *JetsHT_2_Zinc3jet;
    TH1D *JetsHT_2_Zinc4jet;
    TH1D *JetsHT_2_Zinc5jet;
    TH1D *JetsHT_2_Zinc6jet;

    TH1D *genJetsHT_Zinc1jet;
    TH1D *genJetsHT_Zinc2jet;
    TH1D *genJetsHT_Zinc3jet;
    TH1D *genJetsHT_20_Zinc1jet;
    TH1D *genJetsHT_20_Zinc2jet;
    TH1D *genJetsHT_20_Zinc3jet;
    TH1D *genJetsHT_20_30_Zinc1jet;
    TH1D *genJetsHT_20_30_Zinc2jet;
    TH1D *genJetsHT_20_30_Zinc3jet;
    TH1D *genJetsHT_Zinc4jet;
    TH1D *genJetsHT_Zinc5jet;
    TH1D *genJetsHT_Zinc6jet;
    TH1D *genJetsHT_1_Zinc1jet;
    TH1D *genJetsHT_1_Zinc2jet;
    TH1D *genJetsHT_1_Zinc3jet;
    TH1D *genJetsHT_1_Zinc4jet;
    TH1D *genJetsHT_1_Zinc5jet;
    TH1D *genJetsHT_1_Zinc6jet;
    TH1D *genJetsHT_2_Zinc1jet;
    TH1D *genJetsHT_2_Zinc2jet;
    TH1D *genJetsHT_2_Zinc3jet;
    TH1D *genJetsHT_2_Zinc4jet;
    TH1D *genJetsHT_2_Zinc5jet;
    TH1D *genJetsHT_2_Zinc6jet;

    TH2D *hresponseJetsHT_Zinc1jet;
    TH2D *hresponseJetsHT_Zinc2jet;
    TH2D *hresponseJetsHT_Zinc3jet;
    TH2D *hresponseJetsHT_20_Zinc1jet;
    TH2D *hresponseJetsHT_20_Zinc2jet;
    TH2D *hresponseJetsHT_20_Zinc3jet;
    TH2D *hresponseJetsHT_20_30_Zinc1jet;
    TH2D *hresponseJetsHT_20_30_Zinc2jet;
    TH2D *hresponseJetsHT_20_30_Zinc3jet;
    TH2D *hresponseJetsHT_Zinc4jet;
    TH2D *hresponseJetsHT_Zinc5jet;
    TH2D *hresponseJetsHT_Zinc6jet;
    TH2D *hresponseJetsHT_1_Zinc1jet;
    TH2D *hresponseJetsHT_1_Zinc2jet;
    TH2D *hresponseJetsHT_1_Zinc3jet;
    TH2D *hresponseJetsHT_1_Zinc4jet;
    TH2D *hresponseJetsHT_1_Zinc5jet;
    TH2D *hresponseJetsHT_1_Zinc6jet;
    TH2D *hresponseJetsHT_2_Zinc1jet;
    TH2D *hresponseJetsHT_2_Zinc2jet;
    TH2D *hresponseJetsHT_2_Zinc3jet;
    TH2D *hresponseJetsHT_2_Zinc4jet;
    TH2D *hresponseJetsHT_2_Zinc5jet;
    TH2D *hresponseJetsHT_2_Zinc6jet;

    //--- Jet eta
    TH1D *FirstJetEta_Zinc1jet;
    TH1D *FirstJetEta_2_Zinc1jet;
    TH1D *SecondJetEta_Zinc2jet;
    TH1D *SecondJetEta_2_Zinc2jet;
    TH1D *ThirdJetEta_Zinc3jet;
    TH1D *ThirdJetEta_2_Zinc3jet;
    TH1D *FourthJetEta_Zinc4jet;
    TH1D *FourthJetEta_2_Zinc4jet;
    TH1D *FifthJetEta_Zinc5jet;
    TH1D *FifthJetEta_2_Zinc5jet;
    TH1D *SixthJetEta_Zinc6jet;
    TH1D *SixthJetEta_2_Zinc6jet;
    
    TH1D *genFirstJetEta_Zinc1jet;
    TH1D *genFirstJetEta_2_Zinc1jet;
    TH1D *genSecondJetEta_Zinc2jet;
    TH1D *genSecondJetEta_2_Zinc2jet;
    TH1D *genThirdJetEta_Zinc3jet;
    TH1D *genThirdJetEta_2_Zinc3jet;
    TH1D *genFourthJetEta_Zinc4jet;
    TH1D *genFourthJetEta_2_Zinc4jet;
    TH1D *genFifthJetEta_Zinc5jet;
    TH1D *genFifthJetEta_2_Zinc5jet;
    TH1D *genSixthJetEta_Zinc6jet;
    TH1D *genSixthJetEta_2_Zinc6jet;

    TH2D *hresponseFirstJetEta_Zinc1jet;
    TH2D *hresponseFirstJetEta_2_Zinc1jet;
    TH2D *hresponseSecondJetEta_Zinc2jet;
    TH2D *hresponseSecondJetEta_2_Zinc2jet;
    TH2D *hresponseThirdJetEta_Zinc3jet;
    TH2D *hresponseThirdJetEta_2_Zinc3jet;
    TH2D *hresponseFourthJetEta_Zinc4jet;
    TH2D *hresponseFourthJetEta_2_Zinc4jet;
    TH2D *hresponseFifthJetEta_Zinc5jet;
    TH2D *hresponseFifthJetEta_2_Zinc5jet;
    TH2D *hresponseSixthJetEta_Zinc6jet;
    TH2D *hresponseSixthJetEta_2_Zinc6jet;
    //***************************** End Basic plots for Wjets *****************************//


    //***************************** Additional plots *****************************//
    //--- dijet Rapidity
    TH1D *dRapidityJets_Zinc2jet;
    TH1D *dRapidityJets_Zinc3jet;
    TH1D *dRapidityJets_Zinc4jet;
    TH1D *dRapidityJets_2_Zinc2jet;
    TH1D *dRapidityJets_2_Zinc3jet;
    TH1D *dRapidityJets_2_Zinc4jet;
    
    TH1D *gendRapidityJets_Zinc2jet;
    TH1D *gendRapidityJets_Zinc3jet;
    TH1D *gendRapidityJets_Zinc4jet;
    TH1D *gendRapidityJets_2_Zinc2jet;
    TH1D *gendRapidityJets_2_Zinc3jet;
    TH1D *gendRapidityJets_2_Zinc4jet;
    
    TH2D *hresponsedRapidityJets_Zinc2jet;
    TH2D *hresponsedRapidityJets_Zinc3jet;
    TH2D *hresponsedRapidityJets_Zinc4jet;
    TH2D *hresponsedRapidityJets_2_Zinc2jet;
    TH2D *hresponsedRapidityJets_2_Zinc3jet;
    TH2D *hresponsedRapidityJets_2_Zinc4jet;
    
    //--- dijet Rapidity FB
    TH1D *dRapidityJetsFB_Zinc2jet;
    TH1D *dRapidityJetsFB_Zinc3jet;
    TH1D *dRapidityJetsFB_Zinc4jet;
    TH1D *dRapidityJetsFB_2_Zinc2jet;
    TH1D *dRapidityJetsFB_2_Zinc3jet;
    TH1D *dRapidityJetsFB_2_Zinc4jet;
    
    TH1D *gendRapidityJetsFB_Zinc2jet;
    TH1D *gendRapidityJetsFB_Zinc3jet;
    TH1D *gendRapidityJetsFB_Zinc4jet;
    TH1D *gendRapidityJetsFB_2_Zinc2jet;
    TH1D *gendRapidityJetsFB_2_Zinc3jet;
    TH1D *gendRapidityJetsFB_2_Zinc4jet;
    
    TH2D *hresponsedRapidityJetsFB_Zinc2jet;
    TH2D *hresponsedRapidityJetsFB_Zinc3jet;
    TH2D *hresponsedRapidityJetsFB_Zinc4jet;
    TH2D *hresponsedRapidityJetsFB_2_Zinc2jet;
    TH2D *hresponsedRapidityJetsFB_2_Zinc3jet;
    TH2D *hresponsedRapidityJetsFB_2_Zinc4jet;

    //--- dijet Rapidity 1,3
    TH1D *dRapidityJets_First_Third_Zinc3jet;
    TH1D *dRapidityJets_First_Third_Zinc4jet;
    TH1D *dRapidityJets_2_First_Third_Zinc3jet;
    TH1D *dRapidityJets_2_First_Third_Zinc4jet;

    TH1D *gendRapidityJets_First_Third_Zinc3jet;
    TH1D *gendRapidityJets_First_Third_Zinc4jet;
    TH1D *gendRapidityJets_2_First_Third_Zinc3jet;
    TH1D *gendRapidityJets_2_First_Third_Zinc4jet;
    
    TH2D *hresponsedRapidityJets_First_Third_Zinc3jet;
    TH2D *hresponsedRapidityJets_First_Third_Zinc4jet;
    TH2D *hresponsedRapidityJets_2_First_Third_Zinc3jet;
    TH2D *hresponsedRapidityJets_2_First_Third_Zinc4jet;
    
    //--- dijet Rapidity 2,3
    TH1D *dRapidityJets_Second_Third_Zinc3jet;
    TH1D *dRapidityJets_Second_Third_Zinc4jet;
    TH1D *dRapidityJets_2_Second_Third_Zinc3jet;
    TH1D *dRapidityJets_2_Second_Third_Zinc4jet;
    
    TH1D *gendRapidityJets_Second_Third_Zinc3jet;
    TH1D *gendRapidityJets_Second_Third_Zinc4jet;
    TH1D *gendRapidityJets_2_Second_Third_Zinc3jet;
    TH1D *gendRapidityJets_2_Second_Third_Zinc4jet;

    TH2D *hresponsedRapidityJets_Second_Third_Zinc3jet;
    TH2D *hresponsedRapidityJets_Second_Third_Zinc4jet;
    TH2D *hresponsedRapidityJets_2_Second_Third_Zinc3jet;
    TH2D *hresponsedRapidityJets_2_Second_Third_Zinc4jet;

    //--- dPhi dijet
    TH1D *dPhiJets_Zinc2jet;
    TH1D *dPhiJets_Zinc3jet;
    TH1D *dPhiJets_Zinc4jet;
    TH1D *dPhiJets_2_Zinc2jet;
    TH1D *dPhiJets_2_Zinc3jet;
    TH1D *dPhiJets_2_Zinc4jet;
    
    TH1D *gendPhiJets_Zinc2jet;
    TH1D *gendPhiJets_Zinc3jet;
    TH1D *gendPhiJets_Zinc4jet;
    TH1D *gendPhiJets_2_Zinc2jet;
    TH1D *gendPhiJets_2_Zinc3jet;
    TH1D *gendPhiJets_2_Zinc4jet;
    
    TH2D *hresponsedPhiJets_Zinc2jet;
    TH2D *hresponsedPhiJets_Zinc3jet;
    TH2D *hresponsedPhiJets_Zinc4jet;
    TH2D *hresponsedPhiJets_2_Zinc2jet;
    TH2D *hresponsedPhiJets_2_Zinc3jet;
    TH2D *hresponsedPhiJets_2_Zinc4jet;
    
    //--- dPhi dijet FB
    TH1D *dPhiJetsFB_Zinc2jet;
    TH1D *dPhiJetsFB_Zinc3jet;
    TH1D *dPhiJetsFB_Zinc4jet;
    TH1D *dPhiJetsFB_2_Zinc2jet;
    TH1D *dPhiJetsFB_2_Zinc3jet;
    TH1D *dPhiJetsFB_2_Zinc4jet;
    
    TH1D *gendPhiJetsFB_Zinc2jet;
    TH1D *gendPhiJetsFB_Zinc3jet;
    TH1D *gendPhiJetsFB_Zinc4jet;
    TH1D *gendPhiJetsFB_2_Zinc2jet;
    TH1D *gendPhiJetsFB_2_Zinc3jet;
    TH1D *gendPhiJetsFB_2_Zinc4jet;

    TH2D *hresponsedPhiJetsFB_Zinc2jet;
    TH2D *hresponsedPhiJetsFB_Zinc3jet;
    TH2D *hresponsedPhiJetsFB_Zinc4jet;
    TH2D *hresponsedPhiJetsFB_2_Zinc2jet;
    TH2D *hresponsedPhiJetsFB_2_Zinc3jet;
    TH2D *hresponsedPhiJetsFB_2_Zinc4jet;
    
    //--- dPhi (muon, nth_jet)
    TH1D *dPhiLepJet1_Zinc1jet;
    TH1D *dPhiLepJet2_Zinc2jet;
    TH1D *dPhiLepJet3_Zinc3jet;
    TH1D *dPhiLepJet4_Zinc4jet;
    TH1D *dPhiLepJet5_Zinc5jet;
    TH1D *dPhiLepJet1_2_Zinc1jet;
    TH1D *dPhiLepJet2_2_Zinc2jet;
    TH1D *dPhiLepJet3_2_Zinc3jet;
    TH1D *dPhiLepJet4_2_Zinc4jet;
    TH1D *dPhiLepJet5_2_Zinc5jet;
    
    TH1D *gendPhiLepJet1_Zinc1jet;
    TH1D *gendPhiLepJet2_Zinc2jet;
    TH1D *gendPhiLepJet3_Zinc3jet;
    TH1D *gendPhiLepJet4_Zinc4jet;
    TH1D *gendPhiLepJet5_Zinc5jet;
    TH1D *gendPhiLepJet1_2_Zinc1jet;
    TH1D *gendPhiLepJet2_2_Zinc2jet;
    TH1D *gendPhiLepJet3_2_Zinc3jet;
    TH1D *gendPhiLepJet4_2_Zinc4jet;
    TH1D *gendPhiLepJet5_2_Zinc5jet;
    
    TH2D *hresponsedPhiLepJet1_Zinc1jet;
    TH2D *hresponsedPhiLepJet2_Zinc2jet;
    TH2D *hresponsedPhiLepJet3_Zinc3jet;
    TH2D *hresponsedPhiLepJet4_Zinc4jet;
    TH2D *hresponsedPhiLepJet5_Zinc5jet;
    TH2D *hresponsedPhiLepJet1_2_Zinc1jet;
    TH2D *hresponsedPhiLepJet2_2_Zinc2jet;
    TH2D *hresponsedPhiLepJet3_2_Zinc3jet;
    TH2D *hresponsedPhiLepJet4_2_Zinc4jet;
    TH2D *hresponsedPhiLepJet5_2_Zinc5jet;
    	
	//--- dR (muon, closest jet)
	TH1D *dRLepCloseJet_Zinc1jet;
	TH1D *dRLepCloseJet_2_Zinc1jet;
	
	TH1D *gendRLepCloseJet_Zinc1jet;
	TH1D *gendRLepCloseJet_2_Zinc1jet;
	
	TH2D *hresponsedRLepCloseJet_Zinc1jet;
	TH2D *hresponsedRLepCloseJet_2_Zinc1jet;
	
	TH1D *dRLepCloseJetCo_Zinc1jet;
	TH1D *dRLepCloseJetCo_2_Zinc1jet;
	
	TH1D *gendRLepCloseJetCo_Zinc1jet;
	TH1D *gendRLepCloseJetCo_2_Zinc1jet;
	
	TH2D *hresponsedRLepCloseJetCo_Zinc1jet;
	TH2D *hresponsedRLepCloseJetCo_2_Zinc1jet;
	
	TH1D *dRptmin100LepCloseJetCo300dR04_Zinc1jet;
	TH1D *dRptmin100LepCloseJetCo300dR04_2_Zinc1jet;
	
	TH1D *gendRptmin100LepCloseJetCo300dR04_Zinc1jet;
	TH1D *gendRptmin100LepCloseJetCo300dR04_2_Zinc1jet;
	
	TH2D *hresponsedRptmin100LepCloseJetCo300dR04_Zinc1jet;
	TH2D *hresponsedRptmin100LepCloseJetCo300dR04_2_Zinc1jet;

	//--- dR (muon, closest jet in dijet)
	TH1D *dRptmin100LepCloseDiJetCo300dR04_Zinc2jet;
	TH1D *dRptmin100LepCloseDiJetCo300dR04_2_Zinc2jet;

	TH1D *gendRptmin100LepCloseDiJetCo300dR04_Zinc2jet;
	TH1D *gendRptmin100LepCloseDiJetCo300dR04_2_Zinc2jet;

	TH2D *hresponsedRptmin100LepCloseDiJetCo300dR04_Zinc2jet;
	TH2D *hresponsedRptmin100LepCloseDiJetCo300dR04_2_Zinc2jet;
	
    //--- dR dijet
    TH1D *dRJets_Zinc2jet;
    TH1D *dRJets_2_Zinc2jet;

    TH1D *gendRJets_Zinc2jet;
    TH1D *gendRJets_2_Zinc2jet;

    TH2D *hresponsedRJets_Zinc2jet;
    TH2D *hresponsedRJets_2_Zinc2jet;

    //--- dijet mass
    TH1D *diJetMass_Zinc2jet;
    TH1D *diJetMass_Zinc3jet;
    TH1D *diJetMass_Zinc4jet;
    TH1D *diJetMass_2_Zinc2jet;
    TH1D *diJetMass_2_Zinc3jet;
    TH1D *diJetMass_2_Zinc4jet;

    TH1D *gendiJetMass_Zinc2jet;
    TH1D *gendiJetMass_Zinc3jet;
    TH1D *gendiJetMass_Zinc4jet;
    TH1D *gendiJetMass_2_Zinc2jet;
    TH1D *gendiJetMass_2_Zinc3jet;
    TH1D *gendiJetMass_2_Zinc4jet;

    TH2D *hresponsediJetMass_Zinc2jet;
    TH2D *hresponsediJetMass_Zinc3jet;
    TH2D *hresponsediJetMass_Zinc4jet;
    TH2D *hresponsediJetMass_2_Zinc2jet;
    TH2D *hresponsediJetMass_2_Zinc3jet;
    TH2D *hresponsediJetMass_2_Zinc4jet;

    //--- dijet pt
    TH1D *diJetPt_Zinc2jet;
    TH1D *diJetPt_Zinc3jet;
    TH1D *diJetPt_Zinc4jet;
    TH1D *diJetPt_2_Zinc2jet;
    TH1D *diJetPt_2_Zinc3jet;
    TH1D *diJetPt_2_Zinc4jet;

    TH1D *gendiJetPt_Zinc2jet;
    TH1D *gendiJetPt_Zinc3jet;
    TH1D *gendiJetPt_Zinc4jet;
    TH1D *gendiJetPt_2_Zinc2jet;
    TH1D *gendiJetPt_2_Zinc3jet;
    TH1D *gendiJetPt_2_Zinc4jet;

    TH2D *hresponsediJetPt_Zinc2jet;
    TH2D *hresponsediJetPt_Zinc3jet;
    TH2D *hresponsediJetPt_Zinc4jet;
    TH2D *hresponsediJetPt_2_Zinc2jet;
    TH2D *hresponsediJetPt_2_Zinc3jet;
    TH2D *hresponsediJetPt_2_Zinc4jet;

    //--- mean number of jets
    TH2D *MeanNJetsHT_Zinc1jet;
    TH2D *MeanNJetsHT_Zinc2jet;
    TH2D *MeanNJetsdRapidity_Zinc2jet;
    TH2D *MeanNJetsdRapidityFB_Zinc2jet;
    TH2D *genMeanNJetsHT_Zinc1jet;
    TH2D *genMeanNJetsHT_Zinc2jet;
    TH2D *genMeanNJetsdRapidity_Zinc2jet;
    TH2D *genMeanNJetsdRapidityFB_Zinc2jet;

    //--- Jets Rapidity ---
    TH1D *FirstJetAbsRapidity_Zinc1jet;
    TH1D *SecondJetAbsRapidity_Zinc2jet;
    TH1D *ThirdJetAbsRapidity_Zinc3jet;
    TH1D *FourthJetAbsRapidity_Zinc4jet;
    TH1D *FifthJetAbsRapidity_Zinc5jet;
    TH1D *genFirstJetAbsRapidity_Zinc1jet;
    TH1D *genSecondJetAbsRapidity_Zinc2jet;
    TH1D *genThirdJetAbsRapidity_Zinc3jet;
    TH1D *genFourthJetAbsRapidity_Zinc4jet;
    TH1D *genFifthJetAbsRapidity_Zinc5jet;
    TH2D *hresponseFirstJetAbsRapidity_Zinc1jet;
    TH2D *hresponseSecondJetAbsRapidity_Zinc2jet;
    TH2D *hresponseThirdJetAbsRapidity_Zinc3jet;
    TH2D *hresponseFourthJetAbsRapidity_Zinc4jet;
    TH2D *hresponseFifthJetAbsRapidity_Zinc5jet;
    
    TH1D *FirstJetAbsRapidity_2_Zinc1jet;
    TH1D *SecondJetAbsRapidity_2_Zinc2jet;
    TH1D *ThirdJetAbsRapidity_2_Zinc3jet;
    TH1D *FourthJetAbsRapidity_2_Zinc4jet;
    TH1D *FifthJetAbsRapidity_2_Zinc5jet;
    TH1D *genFirstJetAbsRapidity_2_Zinc1jet;
    TH1D *genSecondJetAbsRapidity_2_Zinc2jet;
    TH1D *genThirdJetAbsRapidity_2_Zinc3jet;
    TH1D *genFourthJetAbsRapidity_2_Zinc4jet;
    TH1D *genFifthJetAbsRapidity_2_Zinc5jet;
    TH2D *hresponseFirstJetAbsRapidity_2_Zinc1jet;
    TH2D *hresponseSecondJetAbsRapidity_2_Zinc2jet;
    TH2D *hresponseThirdJetAbsRapidity_2_Zinc3jet;
    TH2D *hresponseFourthJetAbsRapidity_2_Zinc4jet;
    TH2D *hresponseFifthJetAbsRapidity_2_Zinc5jet;

    
    // Jets Rapidity -- not used
    TH1D *FirstJetRapidityFull_Zinc1jet;
    TH1D *SecondJetRapidityFull_Zinc2jet;
    TH1D *ThirdJetRapidityFull_Zinc3jet;
    TH1D *FourthJetRapidityFull_Zinc4jet;
    TH1D *genFirstJetRapidityFull_Zinc1jet;
    TH1D *genSecondJetRapidityFull_Zinc2jet;
    TH1D *genThirdJetRapidityFull_Zinc3jet;
    TH1D *genFourthJetRapidityFull_Zinc4jet;

    // Jets Mass -- not used
    TH1D *FirstJetmass_Zinc1jet;
    TH1D *SecondJetmass_Zinc2jet;
    TH1D *ThirdJetmass_Zinc3jet;
    TH1D *FourthJetmass_Zinc4jet;
    TH1D *FirstJetmass_1_Zinc1jet;
    TH1D *SecondJetmass_1_Zinc2jet;
    TH1D *ThirdJetmass_1_Zinc3jet;
    TH1D *FourthJetmass_1_Zinc4jet;

    // mean number of jets -- not used
    TH1D *MeanNJetsHT_1D_Zinc1jet;
    TH1D *MeanNJetsHT_1D_Zinc2jet;
    TH1D *MeanNJetsdRapidity_1D_Zinc2jet;
    TH1D *MeanNJetsdRapidityFB_1D_Zinc2jet;
    TH1D *genMeanNJetsHT_1D_Zinc1jet;
    TH1D *genMeanNJetsHT_1D_Zinc2jet;
    TH1D *genMeanNJetsdRapidity_1D_Zinc2jet;
    TH1D *genMeanNJetsdRapidityFB_1D_Zinc2jet;
    //***************************** end additional plots *****************************//


    TH1D *NumberPFcandidates;
    TH1D *ZMass_lowDeltaR;

    TH1D *ZMassAllPassLep;
    TH1D *AllPassLepID;
    TH1D *AllPassWithMassCutLepID;
    TH1D *AllPassWithMassCutLepIDCharge;

    TH1D *ZMass_Zinc0jet;
    TH1D *ZMass_Zinc1jet;
    TH1D *ZMass_Zinc2jet;
    TH1D *ZMass_Zinc3jet;
    TH1D *ZMass_Zinc4jet;
    TH1D *ZMass_Zinc5jet;
    TH1D *ZMass_Zinc6jet;
    TH1D *genZMass_Zinc0jet;
    TH1D *genZMass_Zinc1jet;
    TH1D *genZMass_Zinc2jet;
    TH1D *genZMass_Zinc3jet;
    TH1D *genZMass_Zinc4jet;
    TH1D *genZMass_Zinc5jet;
    TH1D *genZMass_Zinc6jet;
    TH1D *ZMass_Zexc0jet;
    TH1D *ZMass_Zexc1jet;
    TH1D *ZMass_Zexc2jet;
    TH1D *ZMass_Zexc3jet;
    TH1D *ZMass_Zexc4jet;
    TH1D *ZMass_Zexc5jet;
    TH1D *ZMass_Zexc6jet;

    TH1D *ZPt_Zinc0jet;
    // TH1D *ZPt_Zinc1jet;
    // TH1D *ZPt_Zinc2jet;
    // TH1D *ZPt_Zinc3jet;
    // TH1D *ZPt_Zinc4jet;
    TH1D *ZPt_Zinc5jet;
    TH1D *ZPt_Zinc6jet;

    TH1D *genZPt_Zinc0jet;
    // TH1D *genZPt_Zinc1jet;
    // TH1D *genZPt_Zinc2jet;
    // TH1D *genZPt_Zinc3jet;
    // TH1D *genZPt_Zinc4jet;
    TH1D *genZPt_Zinc5jet;
    TH1D *genZPt_Zinc6jet;


    /////


    TH1D *ZPt_Zexc0jet;
    TH1D *ZPt_Zexc1jet;
    TH1D *ZPt_Zexc2jet;
    TH1D *ZPt_Zexc3jet;
    TH1D *ZPt_Zexc4jet;
    TH1D *ZPt_Zexc5jet;
    TH1D *ZPt_Zexc6jet;
    TH1D *ZRapidity_Zinc0jet;
    TH1D *ZRapidity_Zinc1jet;
    TH1D *ZRapidity_Zinc2jet;
    TH1D *ZRapidity_Zinc3jet;
    TH1D *ZRapidity_Zinc4jet;
    TH1D *ZRapidity_Zinc5jet;
    TH1D *ZRapidity_Zinc6jet;
    TH1D *genZRapidity_Zinc0jet;
    TH1D *genZRapidity_Zinc1jet;
    TH1D *genZRapidity_Zinc2jet;
    TH1D *genZRapidity_Zinc3jet;
    TH1D *genZRapidity_Zinc4jet;
    TH1D *genZRapidity_Zinc5jet;
    TH1D *genZRapidity_Zinc6jet;
    TH1D *ZRapidity_Zexc0jet;
    TH1D *ZRapidity_Zexc1jet;
    TH1D *ZRapidity_Zexc2jet;
    TH1D *ZRapidity_Zexc3jet;
    TH1D *ZRapidity_Zexc4jet;
    TH1D *ZRapidity_Zexc5jet;
    TH1D *ZRapidity_Zexc6jet;
    TH1D *ZEta_Zinc0jet;
    TH1D *ZEta_Zinc1jet;
    TH1D *ZEta_Zinc2jet;
    TH1D *ZEta_Zinc3jet;
    TH1D *ZEta_Zinc4jet;
    TH1D *ZEta_Zinc5jet;
    TH1D *ZEta_Zinc6jet;
    TH1D *genZEta_Zinc0jet;
    TH1D *genZEta_Zinc1jet;
    TH1D *genZEta_Zinc2jet;
    TH1D *genZEta_Zinc3jet;
    TH1D *genZEta_Zinc4jet;
    TH1D *genZEta_Zinc5jet;
    TH1D *genZEta_Zinc6jet;
    TH1D *ZEta_Zexc0jet;
    TH1D *ZEta_Zexc1jet;
    TH1D *ZEta_Zexc2jet;
    TH1D *ZEta_Zexc3jet;
    TH1D *ZEta_Zexc4jet;
    TH1D *ZEta_Zexc5jet;
    TH1D *ZEta_Zexc6jet;

    TH1D *lepEta_Zinc0jet;
    TH1D *lepEta_Zinc1jet;
    TH1D *lepEta_Zinc2jet;
    TH1D *lepEta_Zinc3jet;
    TH1D *lepEta_Zinc4jet;
    TH1D *lepEta_Zinc5jet;
    TH1D *lepPhi_Zexc0jet;
    TH1D *lepPhi_Zexc1jet;
    TH1D *lepPhi_Zexc2jet;
    TH1D *lepPhi_Zexc3jet;
    TH1D *lepPhi_Zexc4jet;
    TH1D *lepPhi_Zexc5jet;
    TH1D *lepPhi_Zinc0jet;
    TH1D *lepPhi_Zinc1jet;
    TH1D *GLepBareEtaZinc0jet;
    TH1D *GLepBareEtaZinc1jet;
    TH1D *GLepBareEtaZinc2jet;
    TH1D *lepEta_Zexc0jet;
    TH1D *lepEta_Zexc1jet;
    TH1D *lepEta_Zexc2jet;
    TH1D *lepEta_Zexc3jet;
    TH1D *lepEta_Zexc4jet;
    TH1D *lepEta_Zexc5jet;
    TH2D *lepEtaEta_Zinc0jet;

    TH1D *lepChargePlusEta_Zinc1jet;
    TH1D *lepChargeMinusEta_Zinc1jet;
    TH1D *lepChargePlusPhi_Zinc1jet;
    TH1D *lepChargeMinusPhi_Zinc1jet;

    TH1D *FirstJetEtaFull_Zinc1jet;
    TH1D *SecondJetEtaFull_Zinc2jet;
    TH1D *ThirdJetEtaFull_Zinc3jet;
    TH1D *FourthJetEtaFull_Zinc4jet;
    TH1D *FifthJetEtaFull_Zinc5jet;
    TH1D *SixthJetEtaFull_Zinc6jet;

    TH1D *FirstJetEta_Zexc1jet;
    TH1D *SecondJetEta_Zexc2jet;

    TH1D *AllJetEta_Zinc1jet;
    TH1D *AllJetEta_Zinc2jet;
    TH1D *AllJetEta_Zinc3jet;
    TH1D *AllJetEta_Zinc4jet;

    TH1D *AllJetAK8Eta_Zinc1jet;

    TH1D *FirstJetPhi_Zinc1jet;
    TH1D *SecondJetPhi_Zinc2jet;
    TH1D *ThirdJetPhi_Zinc3jet;
    TH1D *FourthJetPhi_Zinc4jet;
    TH1D *FifthJetPhi_Zinc5jet;
    TH1D *SixthJetPhi_Zinc6jet;
    TH1D *FirstJetPhi_Zexc1jet;
    TH1D *SecondJetPhi_Zexc2jet;

    TH1D *AllJetPhi_Zinc1jet;
    TH1D *AllJetPhi_Zinc2jet;
    TH1D *AllJetPhi_Zinc3jet;
    TH1D *AllJetPhi_Zinc4jet;

    TH1D *AllJetAK8Phi_Zinc1jet;

    TH1D *lepPt_Zinc0jet;
    TH1D *lepPt_Zinc1jet;
    TH1D *lepPt_Zinc2jet;
    TH1D *lepPt_Zinc3jet;
    TH1D *lepPt_Zinc4jet;
    TH1D *lepPt_Zinc5jet;
    TH1D *GLepBarePtZinc0jet;
    TH1D *GLepBarePtZinc1jet;
    TH1D *GLepBarePtZinc2jet;
    //andrew
    TH1D *genMT_Zinc1jet;
    TH2D *hresponseLepPt_Zinc1jet;
    //    TH2D *hresponseMET_Zinc1jet;
    TH2D *hresponseMT_Zinc1jet;
    
    TH1D *lepPt_Zexc0jet;
    TH1D *lepPt_Zexc1jet;
    TH1D *lepPt_Zexc2jet;
    TH1D *lepPt_Zexc3jet;
    TH1D *lepPt_Zexc4jet;
    TH1D *lepPt_Zexc5jet;

    TH1D *lepChargePlusPt_Zinc1jet;
    TH1D *lepChargeMinusPt_Zinc1jet;

    TH1D *dPhiLeptons_Zexc0jet;
    TH1D *dPhiLeptons_Zexc1jet;
    TH1D *dPhiLeptons_Zexc2jet;
    TH1D *dPhiLeptons_Zexc3jet;
    TH1D *dPhiLeptons_Zexc4jet;
    TH1D *dPhiLeptons_Zexc5jet;
    TH1D *dPhiLeptons_Zinc0jet;
    TH1D *dPhiLeptons_Zinc1jet;
    TH1D *dPhiLeptons_Zinc2jet;
    TH1D *dPhiLeptons_Zinc3jet;
    TH1D *dPhiLeptons_Zinc4jet;
    TH1D *dPhiLeptons_Zinc5jet;
    TH1D *dEtaLeptons_Zexc0jet;
    TH1D *dEtaLeptons_Zexc1jet;
    TH1D *dEtaLeptons_Zexc2jet;
    TH1D *dEtaLeptons_Zexc3jet;
    TH1D *dEtaLeptons_Zexc4jet;
    TH1D *dEtaLeptons_Zexc5jet;
    TH1D *dEtaLeptons_Zinc0jet;
    TH1D *dEtaLeptons_Zinc1jet;
    TH1D *dEtaLeptons_Zinc2jet;
    TH1D *dEtaLeptons_Zinc3jet;
    TH1D *dEtaLeptons_Zinc4jet;
    TH1D *dEtaLeptons_Zinc5jet;
    TH1D *dRLeptons_Zinc0jet;
    TH1D *dRLeptons_Zinc1jet;
    TH1D *dRLeptons_Zinc2jet;
    TH1D *dRLeptons_Zinc3jet;
    TH1D *dRLeptons_Zinc4jet;
    TH1D *dRLeptons_Zinc5jet;
    TH1D *SpTLeptons_Zexc0jet;
    TH1D *SpTLeptons_Zexc1jet;
    TH1D *SpTLeptons_Zexc2jet;
    TH1D *SpTLeptons_Zexc3jet;
    TH1D *SpTLeptons_Zexc4jet;
    TH1D *SpTLeptons_Zexc5jet;
    TH1D *genSpTLeptons_Zexc2jet;
    TH1D *SpTLeptons_Zinc0jet;
    TH1D *SpTLeptons_Zinc1jet;
    TH1D *SpTLeptons_Zinc2jet;
    TH1D *SpTLeptons_Zinc3jet;
    TH1D *SpTLeptons_Zinc4jet;
    TH1D *SpTLeptons_Zinc5jet;
    TH1D *genSpTLeptons_Zinc2jet;

    ///

    TH1D *RatioJetPt21_Zinc2jet;
    TH1D *RatioJetPt32_Zinc3jet;
    TH1D *genRatioJetPt21_Zinc2jet;
    TH1D *genRatioJetPt32_Zinc3jet;

    TH1D *FirstJetPt_Zexc1jet;
    TH1D *SecondJetPt_Zexc2jet;
    TH1D *genFirstJetPt_Zexc1jet;
    TH1D *genSecondJetPt_Zexc2jet;
    TH1D *FirstHighestJetPt_Zinc1jet;
    TH1D *FirstHighestJetPt_Zinc2jet;
    TH1D *FirstHighestJetPt_Zinc3jet;
    TH1D *FirstHighestJetPt_Zinc4jet;
    TH1D *FirstHighestJetPt_Zinc5jet;
    TH1D *FirstHighestJetPt_Zinc6jet;
    TH1D *genFirstHighestJetPt_Zinc1jet;
    TH1D *genFirstHighestJetPt_Zinc2jet;
    TH1D *genFirstHighestJetPt_Zinc3jet;
    TH1D *genFirstHighestJetPt_Zinc4jet;
    TH1D *genFirstHighestJetPt_Zinc5jet;
    TH1D *genFirstHighestJetPt_Zinc6jet;

    TH1D *SecondHighestJetPt_Zinc2jet;
    TH1D *SecondHighestJetPt_Zinc3jet;
    TH1D *SecondHighestJetPt_Zinc4jet;
    TH1D *SecondHighestJetPt_Zinc5jet;
    TH1D *SecondHighestJetPt_Zinc6jet;
    TH1D *genSecondHighestJetPt_Zinc2jet;
    TH1D *genSecondHighestJetPt_Zinc3jet;
    TH1D *genSecondHighestJetPt_Zinc4jet;
    TH1D *genSecondHighestJetPt_Zinc5jet;
    TH1D *genSecondHighestJetPt_Zinc6jet;

    TH1D *ThirdHighestJetPt_Zinc3jet;
    TH1D *ThirdHighestJetPt_Zinc4jet;
    TH1D *ThirdHighestJetPt_Zinc5jet;
    TH1D *ThirdHighestJetPt_Zinc6jet;
    TH1D *genThirdHighestJetPt_Zinc3jet;
    TH1D *genThirdHighestJetPt_Zinc4jet;
    TH1D *genThirdHighestJetPt_Zinc5jet;
    TH1D *genThirdHighestJetPt_Zinc6jet;

    TH1D *AllJetPt_Zinc1jet;
    TH1D *AllJetPt_Zinc2jet;
    TH1D *AllJetPt_Zinc3jet;
    TH1D *AllJetPt_Zinc4jet;

    TH1D *AllJetAK8Pt_Zinc1jet;

    TH2D *hPtEtaBackJet_Zexc1jet;
    TH2D *hPtEtaBackJetMVA_Zexc1jet;
    TH2D *FirstJetPtEta_Zinc1jet;
    TH2D *SecondJetPtEta_Zinc2jet;
    TH2D *ThirdJetPtEta_Zinc3jet;
    TH2D *FourthJetPtEta_Zinc4jet;
    TH2D *FifthJetPtEta_Zinc5jet;
    TH2D *SixthJetPtEta_Zinc6jet;
    TH2D *genFirstJetPtEta_Zinc1jet;
    TH2D *genSecondJetPtEta_Zinc2jet;
    TH2D *genThirdJetPtEta_Zinc3jet;
    TH2D *genFourthJetPtEta_Zinc4jet;
    TH2D *genFifthJetPtEta_Zinc5jet;
    TH2D *genSixthJetPtEta_Zinc6jet;

    TH2D *ZNGoodJetsNVtx_Zexc;
    TH1D *ZNGoodJets_Zexc_check;
    TH1D *ZNGoodJets_Zexc_NoWeight;
    TH1D *ZNGoodJets_Zinc_NoWeight;

    TH1D *TwoJetsPtDiff_Zexc2jet;
    TH1D *genTwoJetsPtDiff_Zexc2jet;
    TH1D *JetsMass_Zexc2jet;
    TH1D *genJetsMass_Zexc2jet;

    TH1D *ptBal_Zexc2jet;
    TH1D *genptBal_Zexc2jet;
    TH1D *dPhiJets_Zexc2jet;
    TH1D *gendPhiJets_Zexc2jet;
    TH1D *dEtaJets_Zexc2jet;
    TH1D *gendEtaJets_Zexc2jet;
    TH1D *dEtaFirstJetZ_Zexc2jet;
    TH1D *gendEtaFirstJetZ_Zexc2jet;
    TH1D *dEtaSecondJetZ_Zexc2jet;
    TH1D *gendEtaSecondJetZ_Zexc2jet;
    TH1D *dEtaJet1Plus2Z_Zexc2jet;
    TH1D *gendEtaJet1Plus2Z_Zexc2jet;
    TH1D *PHI_Zexc2jet;
    TH1D *genPHI_Zexc2jet;
    TH1D *PHI_T_Zexc2jet;
    TH1D *genPHI_T_Zexc2jet;
    TH1D *SpT_Zexc2jet;
    TH1D *genSpT_Zexc2jet;
    TH1D *SpTJets_Zexc2jet;
    TH1D *genSpTJets_Zexc2jet;
    TH1D *SPhi_Zexc2jet;
    TH1D *genSPhi_Zexc2jet;

    TH1D *TwoJetsPtDiff_Zinc2jet;
    TH1D *genTwoJetsPtDiff_Zinc2jet;
    TH1D *BestTwoJetsPtDiff_Zinc2jet;
    TH1D *genBestTwoJetsPtDiff_Zinc2jet;
    TH1D *JetsMass_Zinc2jet;
    TH1D *genJetsMass_Zinc2jet;
    TH1D *BestJetsMass_Zinc2jet;
    TH1D *genBestJetsMass_Zinc2jet;

    TH1D *ptBal_Zinc2jet;
    TH1D *genptBal_Zinc2jet;

    TH1D *BestdPhiJets_Zinc2jet;
    TH1D *genBestdPhiJets_Zinc2jet;
    TH1D *dEtaJets_Zinc2jet;
    TH1D *gendEtaJets_Zinc2jet;
    TH1D *dEtaFirstJetZ_Zinc2jet;
    TH1D *gendEtaFirstJetZ_Zinc2jet;
    TH1D *dEtaSecondJetZ_Zinc2jet;
    TH1D *gendEtaSecondJetZ_Zinc2jet;
    TH1D *dEtaJet1Plus2Z_Zinc2jet;
    TH1D *gendEtaJet1Plus2Z_Zinc2jet;
    TH1D *PHI_Zinc2jet;
    TH1D *genPHI_Zinc2jet;
    TH1D *BestPHI_Zinc2jet;
    TH1D *genBestPHI_Zinc2jet;
    TH1D *PHI_T_Zinc2jet;
    TH1D *genPHI_T_Zinc2jet;
    TH1D *BestPHI_T_Zinc2jet;
    TH1D *genBestPHI_T_Zinc2jet;
    TH1D *SpT_Zinc2jet;
    TH1D *genSpT_Zinc2jet;
    TH1D *BestSpT_Zinc2jet;
    TH1D *genBestSpT_Zinc2jet;
    TH1D *SpTJets_Zinc2jet;
    TH1D *genSpTJets_Zinc2jet;
    TH1D *BestSpTJets_Zinc2jet;
    TH1D *genBestSpTJets_Zinc2jet;
    TH1D *SPhi_Zinc2jet;
    TH1D *genSPhi_Zinc2jet;
    TH1D *BestSPhi_Zinc2jet;
    TH1D *genBestSPhi_Zinc2jet;

    //-- low Z pT;
    TH1D *TwoJetsPtDiff_LowPt_Zexc2jet;
    TH1D *genTwoJetsPtDiff_LowPt_Zexc2jet;
    TH1D *JetsMass_LowPt_Zexc2jet;
    TH1D *genJetsMass_LowPt_Zexc2jet;
    TH1D *ptBal_LowPt_Zexc2jet;
    TH1D *genptBal_LowPt_Zexc2jet;
    TH1D *dPhiJets_LowPt_Zexc2jet;
    TH1D *gendPhiJets_LowPt_Zexc2jet;
    TH1D *dPhiLeptons_LowPt_Zexc2jet;
    TH1D *gendPhiLeptons_LowPt_Zexc2jet;
    TH1D *PHI_LowPt_Zexc2jet;
    TH1D *genPHI_LowPt_Zexc2jet;
    TH1D *PHI_T_LowPt_Zexc2jet;
    TH1D *genPHI_T_LowPt_Zexc2jet;
    TH1D *SpT_LowPt_Zexc2jet;
    TH1D *genSpT_LowPt_Zexc2jet;
    TH1D *SpTJets_LowPt_Zexc2jet;
    TH1D *genSpTJets_LowPt_Zexc2jet;
    TH1D *SpTLeptons_LowPt_Zexc2jet;
    TH1D *genSpTLeptons_LowPt_Zexc2jet;
    TH1D *SPhi_LowPt_Zexc2jet;
    TH1D *genSPhi_LowPt_Zexc2jet;

    TH1D *TwoJetsPtDiff_LowPt_Zinc2jet;
    TH1D *genTwoJetsPtDiff_LowPt_Zinc2jet;
    TH1D *BestTwoJetsPtDiff_LowPt_Zinc2jet;
    TH1D *genBestTwoJetsPtDiff_LowPt_Zinc2jet;

    TH1D *JetsMass_LowPt_Zinc2jet;
    TH1D *genJetsMass_LowPt_Zinc2jet;
    TH1D *BestJetsMass_LowPt_Zinc2jet;
    TH1D *genBestJetsMass_LowPt_Zinc2jet;
    TH1D *ptBal_LowPt_Zinc2jet;
    TH1D *genptBal_LowPt_Zinc2jet;
    TH1D *dPhiJets_LowPt_Zinc2jet;
    TH1D *gendPhiJets_LowPt_Zinc2jet;
    TH1D *BestdPhiJets_LowPt_Zinc2jet;
    TH1D *genBestdPhiJets_LowPt_Zinc2jet;
    TH1D *dPhiLeptons_LowPt_Zinc2jet;
    TH1D *gendPhiLeptons_LowPt_Zinc2jet;
    TH1D *PHI_LowPt_Zinc2jet;
    TH1D *genPHI_LowPt_Zinc2jet;
    TH1D *BestPHI_LowPt_Zinc2jet;
    TH1D *genBestPHI_LowPt_Zinc2jet;
    TH1D *PHI_T_LowPt_Zinc2jet;
    TH1D *genPHI_T_LowPt_Zinc2jet;
    TH1D *BestPHI_T_LowPt_Zinc2jet;
    TH1D *genBestPHI_T_LowPt_Zinc2jet;
    TH1D *SpT_LowPt_Zinc2jet;
    TH1D *genSpT_LowPt_Zinc2jet;
    TH1D *BestSpT_LowPt_Zinc2jet;
    TH1D *genBestSpT_LowPt_Zinc2jet;
    TH1D *SpTJets_LowPt_Zinc2jet;
    TH1D *genSpTJets_LowPt_Zinc2jet;
    TH1D *BestSpTJets_LowPt_Zinc2jet;
    TH1D *genBestSpTJets_LowPt_Zinc2jet;
    TH1D *SpTLeptons_LowPt_Zinc2jet;
    TH1D *genSpTLeptons_LowPt_Zinc2jet;
    TH1D *SPhi_LowPt_Zinc2jet;
    TH1D *genSPhi_LowPt_Zinc2jet;
    TH1D *BestSPhi_LowPt_Zinc2jet;
    TH1D *genBestSPhi_LowPt_Zinc2jet;

    //-- low Z pT and low SpT;
    TH1D *PHI_LowSpT_LowPt_Zexc2jet;
    TH1D *genPHI_LowSpT_LowPt_Zexc2jet;
    TH1D *SPhi_LowSpT_LowPt_Zexc2jet;
    TH1D *genSPhi_LowSpT_LowPt_Zexc2jet;

    TH1D *PHI_LowSpT_LowPt_Zinc2jet;
    TH1D *genPHI_LowSpT_LowPt_Zinc2jet;
    TH1D *SPhi_LowSpT_LowPt_Zinc2jet;
    TH1D *genSPhi_LowSpT_LowPt_Zinc2jet;

    //-- low Z pT and high SpT;
    TH1D *PHI_HighSpT_LowPt_Zexc2jet;
    TH1D *genPHI_HighSpT_LowPt_Zexc2jet;
    TH1D *SPhi_HighSpT_LowPt_Zexc2jet;
    TH1D *genSPhi_HighSpT_LowPt_Zexc2jet;

    TH1D *PHI_HighSpT_LowPt_Zinc2jet;
    TH1D *genPHI_HighSpT_LowPt_Zinc2jet;
    TH1D *SPhi_HighSpT_LowPt_Zinc2jet;
    TH1D *genSPhi_HighSpT_LowPt_Zinc2jet;

    //-- low Z pT and low SPhi;
    TH1D *SpT_LowSPhi_LowPt_Zexc2jet;
    TH1D *genSpT_LowSPhi_LowPt_Zexc2jet;

    TH1D *SpT_LowSPhi_LowPt_Zinc2jet;
    TH1D *genSpT_LowSPhi_LowPt_Zinc2jet;

    //-- low Z pT and high SPhi;
    TH1D *SpT_HighSPhi_LowPt_Zexc2jet;
    TH1D *genSpT_HighSPhi_LowPt_Zexc2jet;
    ;
    TH1D *SpT_HighSPhi_LowPt_Zinc2jet;
    TH1D *genSpT_HighSPhi_LowPt_Zinc2jet;

    //-- high Z pT;
    TH1D *ptBal_HighPt_Zexc2jet;
    TH1D *genptBal_HighPt_Zexc2jet;
    TH1D *dPhiJets_HighPt_Zexc2jet;
    TH1D *gendPhiJets_HighPt_Zexc2jet;
    TH1D *dPhiLeptons_HighPt_Zexc2jet;
    TH1D *gendPhiLeptons_HighPt_Zexc2jet;
    TH1D *PHI_HighPt_Zexc2jet;
    TH1D *genPHI_HighPt_Zexc2jet;
    TH1D *PHI_T_HighPt_Zexc2jet;
    TH1D *genPHI_T_HighPt_Zexc2jet;
    TH1D *SpT_HighPt_Zexc2jet;
    TH1D *genSpT_HighPt_Zexc2jet;
    TH1D *SpTJets_HighPt_Zexc2jet;
    TH1D *genSpTJets_HighPt_Zexc2jet;
    TH1D *SpTLeptons_HighPt_Zexc2jet;
    TH1D *genSpTLeptons_HighPt_Zexc2jet;
    TH1D *SPhi_HighPt_Zexc2jet;
    TH1D *genSPhi_HighPt_Zexc2jet;

    TH1D *ptBal_HighPt_Zinc2jet;
    TH1D *genptBal_HighPt_Zinc2jet;
    TH1D *dPhiJets_HighPt_Zinc2jet;
    TH1D *gendPhiJets_HighPt_Zinc2jet;
    TH1D *dPhiLeptons_HighPt_Zinc2jet;
    TH1D *gendPhiLeptons_HighPt_Zinc2jet;
    TH1D *PHI_HighPt_Zinc2jet;
    TH1D *genPHI_HighPt_Zinc2jet;
    TH1D *PHI_T_HighPt_Zinc2jet;
    TH1D *genPHI_T_HighPt_Zinc2jet;
    TH1D *SpT_HighPt_Zinc2jet;
    TH1D *genSpT_HighPt_Zinc2jet;
    TH1D *SpTJets_HighPt_Zinc2jet;
    TH1D *genSpTJets_HighPt_Zinc2jet;
    TH1D *SpTLeptons_HighPt_Zinc2jet;
    TH1D *genSpTLeptons_HighPt_Zinc2jet;
    TH1D *SPhi_HighPt_Zinc2jet;
    TH1D *genSPhi_HighPt_Zinc2jet;

    //-- high Z pT and low SpT
    TH1D *PHI_LowSpT_HighPt_Zexc2jet;
    TH1D *SPhi_LowSpT_HighPt_Zexc2jet;

    TH1D *PHI_LowSpT_HighPt_Zinc2jet;
    TH1D *SPhi_LowSpT_HighPt_Zinc2jet;

    //-- high Z pT and high SpT
    TH1D *PHI_HighSpT_HighPt_Zexc2jet;
    TH1D *SPhi_HighSpT_HighPt_Zexc2jet;

    TH1D *PHI_HighSpT_HighPt_Zinc2jet;
    TH1D *SPhi_HighSpT_HighPt_Zinc2jet;

    //-- high Z pT and low SPhi
    TH1D *SpT_LowSPhi_HighPt_Zexc2jet;
    TH1D *SpT_LowSPhi_HighPt_Zinc2jet;

    //-- high Z pT and high SPhi
    TH1D *SpT_HighSPhi_HighPt_Zexc2jet;
    TH1D *SpT_HighSPhi_HighPt_Zinc2jet;

    //-- low SPhi
    TH1D *SpT_LowSPhi_Zexc2jet;
    TH1D *SpT_LowSPhi_Zinc2jet;

    //-- high SPhi
    TH1D *SpT_HighSPhi_Zexc2jet;
    TH1D *SpT_HighSPhi_Zinc2jet;

    //-- low SpT
    TH1D *PHI_LowSpT_Zexc2jet;
    TH1D *SPhi_LowSpT_Zexc2jet;

    TH1D *PHI_LowSpT_Zinc2jet;
    TH1D *SPhi_LowSpT_Zinc2jet;

    //-- high SpT
    TH1D *PHI_HighSpT_Zexc2jet; 
    TH1D *SPhi_HighSpT_Zexc2jet;

    TH1D *PHI_HighSpT_Zinc2jet; 
    TH1D *SPhi_HighSpT_Zinc2jet;

    //-- gen stuff
    TH1D *gendPhiJetsDeltaR_Zexc2jet;
    TH1D *resdPhiJetsDeltaR_Zexc2jet;
    TH1D *genPHI_TDeltaR_Zexc2jet;
    TH1D *resPHI_TDeltaR_Zexc2jet;
    TH1D *genSpTJetsDeltaR_Zexc2jet;
    TH1D *resSpTJetsDeltaR_Zexc2jet;
    TH1D *genSpTDeltaR_Zexc2jet;
    TH1D *resSpTDeltaR_Zexc2jet;

    TH1D *gendPhiJetsDPS_Zexc2jet;
    TH1D *gendPhiJetsDPSDeltaR_Zexc2jet;
    TH1D *genPHI_TDPS_Zexc2jet;
    TH1D *genPHI_TDPSDeltaR_Zexc2jet;
    TH1D *genSpTJetsDPS_Zexc2jet;
    TH1D *genSpTJetsDPSDeltaR_Zexc2jet;
    TH1D *genSpTDPS_Zexc2jet;
    TH1D *genSpTDPSDeltaR_Zexc2jet;
    TH1D *genSpTDPSPartons_Zexc2jet;


    TH2D *gendPhiJetsDPSDeltaR_ZpT_Zexc2jet;
    TH2D *partonX2D;

    TH1D *gendeltaRjetMu;

    /// additional information
    // Muoisolation

    TH1D *MuDetIsoRhoCorr;
    TH1D *MuPFIsoDBetaCorr;

    TH1D *MuPFIso_Zinc0jet;
    TH1D *MuPFIso_2ndZinc0jet;
    TH1D *MuPFIso_3rdZinc0jet;

    TH1D *deltaRjetMu;
    TH1D *deltaPtjetMu;

    TH1D *Beta;
    TH1D *BetaStar;
    TH2D *ZNGoodJetsBeta_Zexc;
    TH1D *puBeta_JetsMatchGenJets;
    TH1D *puBetaStar_JetsMatchGenJets;
    TH1D *puBeta_JetsNoMatchGenJets;
    TH1D *puBetaStar_JetsNoMatchGenJets;
    TH1D *puMVA;
    TH1D *puMVA_JetsMatchGenJets;
    TH1D *puMVA_JetsNoMatchGenJets;
    TH1D *jetsEta_JetsMatchGenJets;
    TH1D *jetsEta_JetsNoMatchGenJets;
    TH1D *FirstJetdEtaGenReco_Zinc1;
    TH1D *FourthJetdEtaGenReco_Zinc4; 
    TH2D *puMVAvsBeta;
    TH1D *PUWeight;
    TH1D *fullMET;
    TH1D *fullMET_pfMETPFlow;
    TH1D *fullMET_pfMet;
    TH1D *fullMET_pfType1CorrectedMet;
    TH1D *fullMET_pfType1p2CorrectedMet;
    TH1D *fullMT;
    TH2D *METvslepIso;
    TH2D *MTvslepIso;
    TH1D *PUWeight0;
    TH1D *PUWeight1;
    TH1D *MuPlusPt;
    TH1D *MuMinusPt;
    TH1D *MuPlusEta;
    TH1D *MuMinusEta;

    TH1D *MET_Zinc0jet;
    TH1D *MET_Zinc1jet;
    TH1D *MET_Zinc2jet;
    TH1D *MET_Zinc3jet;

    TH1D *METphi_Zinc0jet;
    TH1D *METphi_Zinc1jet;
    TH1D *METphi_Zinc2jet;
    TH1D *METphi_Zinc3jet;

    TH1D *MT_Zinc0jet;
    TH1D *MT_Zinc1jet;
    TH1D *MT_Zinc2jet;
    TH1D *MT_Zinc3jet;

    //--- number of addtion genrated parton
    TH1D *partonsN;         
    TH1D *partonsNWeighted;  
    TH1D *partonsNAfterGenCut;        
    TH1D *partonsNAfterGenCutWeighted; 

    // -- Vector boson jet properties
    TH1D * dEtaBosonJet_Zexc1jet;
    TH1D * dEtaBosonJet_Zinc1jet;
    TH1D * gendEtaBosonJet_Zexc1jet;
    TH1D * gendEtaBosonJet_Zinc1jet;

    // for 2D unfolding --- create 1D histograms for different rapidity ranges
    //int NbinsEta = 9 ;
    //double j_pT_range[9]={30.,50.,70.,100.,130.,170.,220.,350.,1000.};
    //double j_Y_range[9]={0.,0.5,1.0,1.5,2.0,2.5,3.0,3.2,4.7};
    TH1D *FirstJetPt_Zinc1jet_Eta[10];
    TH1D *genFirstJetPt_Zinc1jet_Eta[10];
    TH1D *SecondJetPt_Zinc2jet_Eta[10];
    TH1D *genSecondJetPt_Zinc2jet_Eta[10];

    ClassDef(HistoSet,0);
};

#endif