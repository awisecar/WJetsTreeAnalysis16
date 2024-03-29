#ifndef _functions_H_
#define _functions_H_

#include <iostream>
#include <cstdarg>
#include <vector>
#include <TLorentzVector.h>

using namespace std;

// void barre_de_progression(int);

struct leptonStruct{
    double pt, eta, phi, energy, charge, iso, scEta ;
};

struct jetStruct{
    double pt, eta, phi, energy;
    int patIndex;
    bool isBJet;
	bool passDR04;
	bool passDR02;
    float btagDiscScore;
    bool hasGoodSVIVF;
    bool hasGoodSVSSV;
};

bool LepDescendingOrder(leptonStruct, leptonStruct);
bool JetDescendingOrder(jetStruct, jetStruct);
bool JetYDescendingOrder(TLorentzVector, TLorentzVector);
double deltaRYPhi(TLorentzVector, TLorentzVector);

vector<double> makeVector(int num, ...);
void insertVector(vector<double>& veca, int num, ...);
TH1D* newTH1D(string, string, string, int, double*);
TH1D* newTH1D(string, string, string, int, double, double);
TH1D* newTH1D(string, string, string, vector<double>&);
TH2D* newTH2D(string, string, int, double*, int, double*);
TH2D* newTH2D(string, string, int, double*, int, double, double);
TH2D* newTH2D(string, string, int, double, double, int, double*);
TH2D* newTH2D(string, string, int, double, double, int, double, double);

double phi0to2pi(double);

double deltaPhi(TLorentzVector, TLorentzVector);

double deltaPhi(double, double);
double deltaR(TLorentzVector, TLorentzVector);
double deltaR(double, double, double, double);
double PHI(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);
double PHI_T(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);
double SpTsub(TLorentzVector, TLorentzVector);
double SpT(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);
double SPhi(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);

class record{
    public:
        double ptLow, ptHi, etaLow, etaHi, effi, effiErrorLow, effiErrorHigh;
        record();
        record(double, double, double, double, double, double, double);
        bool belongTo(double, double);
};

class table{
    public:
        table();
        table(string);
        double getEfficiency(double, double, int);
        double getEfficiency(double, double);  
        double getEfficiencyLow(double, double);  
        double getEfficiencyHigh(double, double);  

    private:
        std::vector<record> recd;
};

double SmearJetPt(double recoPt, double genPt, double eta, int smearJet, int year, int jetType);
double SmearJetPtLite(double recoPt, double genPt, double scaleFactor);

void normalizeTH2D(TH2D*);
void bestTwoJetsCandidatesPt(vector<jetStruct>, pair<TLorentzVector, TLorentzVector>&);
void bestTwoJetsCandidatesPhi(vector<jetStruct>, pair<TLorentzVector, TLorentzVector>&);

vector<double> buildVecFineBin(int nStdBin, double arrStdBin[], int factChop);

void splitBinsInTwoForTUnfold(int nBinsOriginal, double originalArray[], double splitArray[]);

void welcomeMessage();
void bTagVetoMessage(int doBJets);

#endif
