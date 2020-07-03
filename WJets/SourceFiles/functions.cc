#include <iostream>
#include <sstream>
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <vector>
#include <cstdarg>

#include "getFilesAndHistograms.h"
#include "functions.h"

using namespace std;

bool LepDescendingOrder(leptonStruct l1, leptonStruct l2){
    return (l1.pt > l2.pt);
}

bool JetDescendingOrder(jetStruct j1, jetStruct j2){
    return (j1.pt > j2.pt);
}

bool JetYDescendingOrder(TLorentzVector tj1, TLorentzVector tj2){
    return (tj1.Rapidity() > tj2.Rapidity());
}

double deltaRYPhi(TLorentzVector j1, TLorentzVector j2){
    double dY = j1.Rapidity() - j2.Rapidity();
    double dPhi = deltaPhi(j1, j2);
    return sqrt(dY * dY + dPhi * dPhi);
}

vector<double> makeVector(int num, ...)
{
    va_list list;
    va_start(list, num);
    vector<double> vec;
    for (int i(0); i < num; i++) {
        double next = va_arg(list, double);
        vec.push_back(next);
    }
    va_end(list);
    return vec;
}

void insertVector(vector<double>& veca, int num, ...)
{
    va_list list;
    va_start(list, num);
    vector<double> vecb;
    for (int i(0); i < num; i++) {
        double next = va_arg(list, double);
        vecb.push_back(next);
    }
    va_end(list);
    veca.insert(veca.end(), vecb.begin(), vecb.end());
}

TH1D* newTH1D(string name, string title, string xTitle, int nBins, double *xBins)
{
    TH1D* hist = new TH1D(name.c_str(), title.c_str(), nBins, xBins);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle("# Events");
    return hist;
}

TH1D* newTH1D(string name, string title, string xTitle, vector<double>& xBinsVect)
{
    int nBins = xBinsVect.size()-1;
    double *xBins = new double[xBinsVect.size()];
    std::copy(xBinsVect.begin(), xBinsVect.end(), xBins);
    TH1D* hist = new TH1D(name.c_str(), title.c_str(), nBins, xBins);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle("# Events");
    delete [] xBins;
    return hist;
}

TH1D* newTH1D(string name, string title, string xTitle, int nBins, double xLow, double xUp){
    TH1D* hist = new TH1D(name.c_str(), title.c_str(), nBins, xLow, xUp);
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle("# Events");
    hist->SetOption("HIST");
    return hist;
}

TH2D* newTH2D(string name, string title, int nBinsX, double *xBins, int nBinsY, double *yBinsY){
    TH2D* hist = new TH2D(name.c_str(), title.c_str(), nBinsX, xBins, nBinsY, yBinsY);
    hist->GetZaxis()->SetTitle("# Events");
    return hist;
}

TH2D* newTH2D(string name, string title, int nBinsX, double *xBins, int nBinsY, double yLow, double yUp){
    TH2D* hist = new TH2D(name.c_str(), title.c_str(), nBinsX, xBins, nBinsY, yLow, yUp);
    hist->GetZaxis()->SetTitle("# Events");
    return hist;
}

TH2D* newTH2D(string name, string title, int nBinsX, double xLow, double xUp, int nBinsY, double *yBins){
    TH2D* hist = new TH2D(name.c_str(), title.c_str(), nBinsX, xLow, xUp, nBinsY, yBins);
    hist->GetZaxis()->SetTitle("# Events");
    return hist;
}

TH2D* newTH2D(string name, string title, int nBinsX, double xLow, double xUp, int nBinsY, double yLow, double yUp){
    TH2D* hist = new TH2D(name.c_str(), title.c_str(), nBinsX, xLow, xUp, nBinsY, yLow, yUp);
    hist->GetZaxis()->SetTitle("# Events");
    hist->SetOption("HIST");
    return hist;
}

double phi0to2pi(double phi){
    double pi = 3.141592653589793238;
    while (phi >= 2.*pi) phi -= 2.*pi;
    while (phi < 0.) phi += 2.*pi;
    return phi;
}

double deltaPhi(TLorentzVector v1, TLorentzVector v2){
    // build the delta Phi angle between the two vectors
    double pi = 3.141592653589793238;
    double phi1 = phi0to2pi(v1.Phi());
    double phi2 = phi0to2pi(v2.Phi());
    double dPhi = phi0to2pi(phi1 - phi2);
    dPhi = (dPhi > (2*pi - dPhi)) ? 2*pi - dPhi : dPhi;
    return dPhi;
} 

double deltaPhi(double Phi1, double Phi2){
    // build the delta Phi angle between the two vectors
    double pi = 3.141592653589793238;
    double phi1 = phi0to2pi(Phi1);
    double phi2 = phi0to2pi(Phi2);
    double dPhi = phi0to2pi(phi1 - phi2);
    dPhi = (dPhi > (2*pi - dPhi)) ? 2*pi - dPhi : dPhi;
    //cout << "      DeltaPhi: " << endl;
    //cout << "      phi1 = " << phi1 << "  phi2 = " << phi2 << endl;
    //cout << "      DeltaPhi = " << dPhi << endl;
    return dPhi;
} 

double deltaR(TLorentzVector v1, TLorentzVector v2){
    double dEta = v1.Eta() - v2.Eta();
    double dPhi = deltaPhi(v1, v2);
    return sqrt(dEta * dEta + dPhi * dPhi);
}

double deltaR(double Phi1, double Eta1, double Phi2, double Eta2){
    //cout << "DeltaR:" << endl;
    //cout << "phi1 = " << Phi1 << "  eta1 = " << Eta1 << "  phi2 = " << Phi2 << "  eta2 = " << Eta2 << endl; 
    double dEta = Eta1 - Eta2;
    double dPhi = deltaPhi(Phi1, Phi2);
    //cout << "   deltaR = " << sqrt(dEta * dEta + dPhi * dPhi) << endl;
    return sqrt(dEta * dEta + dPhi * dPhi);
}

double PHI(TLorentzVector l1, TLorentzVector l2, TLorentzVector j1, TLorentzVector j2){
    // build the angle PHI between the two subsytems (l1+l2, j1+j2) vectors 
    double lPx = (l1.Px() + l2.Px());
    double lPy = (l1.Py() + l2.Py());
    double lPz = (l1.Pz() + l2.Pz());
    double lNorm = sqrt(lPx * lPx + lPy * lPy + lPz * lPz);
    double jPx = (j1.Px() + j2.Px());
    double jPy = (j1.Py() + j2.Py());
    double jPz = (j1.Pz() + j2.Pz());
    double jNorm = sqrt(jPx * jPx + jPy * jPy + jPz * jPz);
    return acos((jPx * lPx + jPy * lPy + jPz * lPz) / (jNorm * lNorm));
}

double PHI_T(TLorentzVector l1, TLorentzVector l2, TLorentzVector j1, TLorentzVector j2){
    // build the angle PHI between the two subsytems (l1+l2, j1+j2) vectors in the transverse plane
    double lPx = (l1.Px() + l2.Px());
    double lPy = (l1.Py() + l2.Py());
    double lNorm = sqrt(lPx * lPx + lPy * lPy);
    double jPx = (j1.Px() + j2.Px());
    double jPy = (j1.Py() + j2.Py());
    double jNorm = sqrt(jPx * jPx + jPy * jPy);
    return acos((jPx * lPx + jPy * lPy) / (jNorm * lNorm));
}

double SpTsub(TLorentzVector v1, TLorentzVector v2){
    return sqrt(pow(v1.Px() + v2.Px(), 2) + pow(v1.Py() + v2.Py(), 2)) / (v1.Pt() + v2.Pt());
}

double SpT(TLorentzVector l1, TLorentzVector l2, TLorentzVector j1, TLorentzVector j2){
    return sqrt( pow(SpTsub(l1, l2), 2) + pow(SpTsub(j1, j2), 2)  ) / sqrt(2.);
} 

double SPhi(TLorentzVector l1, TLorentzVector l2, TLorentzVector j1, TLorentzVector j2){
    return sqrt(deltaPhi(l1, l2) * deltaPhi(l1, l2) + deltaPhi(j1, j2) * deltaPhi(j1, j2)) / sqrt(2.);
}

record::record(): 
    ptLow(0.), ptHi(0.), etaLow(0.), etaHi(0.), effi(0.), effiErrorLow(0.), effiErrorHigh(0.)
{
}

record::record(double pt1, double pt2, double eta1, double eta2, double eff, double effLow, double effHigh):
    ptLow(pt1), ptHi(pt2), etaLow(eta1), etaHi(eta2), effi(eff), effiErrorLow(effLow), effiErrorHigh(effHigh)
{
}

bool record::belongTo(double pt, double eta)
{
    return (pt < ptHi && pt >= ptLow) && (eta < etaHi && eta >= etaLow);
}

table::table()
{
}

table::table(string filename)
{
    ifstream file(filename.c_str());
    if (file) std::cout << filename << " has been found!" << std::endl;
    else std::cout << filename << "has NOT been found..." << std::endl;

    double data[7];
    while(file){
        // can print out the table line by line to see what is being input into "recd"
        // std::cout << "\tentering new line of table:\t\t";
        for (int i(0); i < 7; i++)
        {
            file >> data[i];
            // std::cout << data[i] << "\t";
        }
        // std::cout << std::endl;
        recd.push_back(record(data[2],data[3],data[0],data[1],data[4],data[5],data[6]));
    }
    // std::cout << std::endl;
}

//used for Lep SFs
double table::getEfficiency(double pt, double eta, int sysLepSF){
    //use the hiPtBin value to grab the efficiency of the last pt-bin if the object pt exceeds the max pt in the table
    double hiPtBin = 0;
    //nominal
    if (sysLepSF == 0){
        // std::cout << " ===== table::getEfficiency =====" << std::endl;
        for (unsigned int i=0; i != recd.size(); i++){

            // look inside elements of "recd" as it's being accessed by table::getEfficiency
            // std::cout << "new line of recd:  ";
            // std::cout << recd[i].ptLow << "  " << recd[i].ptHi << "  " << recd[i].etaLow << "  " << recd[i].etaHi << "  ";
            // std::cout << recd[i].effi << "  " << recd[i].effiErrorLow << "  " << recd[i].effiErrorHigh << "  " << std::endl;

            if((recd[i]).belongTo(pt, eta)) return recd[i].effi;
            if((recd[i]).belongTo(0.5*(recd[i].ptHi + recd[i].ptLow), eta)) hiPtBin = recd[i].effi;
        }
        return hiPtBin;
    } 
    //LepSF variation up
    else if (sysLepSF == 1){
        for (unsigned int i=0; i != recd.size(); i++){
            if((recd[i]).belongTo(pt, eta)) return recd[i].effi + recd[i].effiErrorHigh;
            if((recd[i]).belongTo(0.5*(recd[i].ptHi + recd[i].ptLow), eta)) hiPtBin = recd[i].effi + recd[i].effiErrorHigh;
        }
        return hiPtBin;
    }
    //LepSF variation down
    else if (sysLepSF == -1){
        for (unsigned int i=0; i != recd.size(); i++){
            if((recd[i]).belongTo(pt, eta)) return recd[i].effi - recd[i].effiErrorLow;
            if((recd[i]).belongTo(0.5*(recd[i].ptHi + recd[i].ptLow), eta)) hiPtBin = recd[i].effi - recd[i].effiErrorLow;
        }
        return hiPtBin;
    }
    else return 1.;
}

//used for the smear factors for JES
double table::getEfficiency(double pt, double eta){
    double hiPtBin= 0;
    for (unsigned int i=0; i != recd.size(); i++) {
        // if finds the proper bin, then return the corresponding efficiency
        if((recd[i]).belongTo(pt, eta)) return recd[i].effi;
        // if pt of object goes beyond table bounds, then record the last efficiency value in the pt range
        if((recd[i]).belongTo(0.5*(recd[i].ptHi + recd[i].ptLow), eta)) hiPtBin = recd[i].effi;
    }
    return hiPtBin;
} 

double table::getEfficiencyLow(double pt, double eta){
    double hiPtBin= 0;
    for (unsigned int i=0; i != recd.size(); i++) {
        if((recd[i]).belongTo(pt, eta)) return recd[i].effi-recd[i].effiErrorLow;
        if((recd[i]).belongTo(190, eta)) hiPtBin = recd[i].effi;
    }
    return hiPtBin;
}

double table::getEfficiencyHigh(double pt, double eta){
    double hiPtBin= 0;
    for (unsigned int i=0; i != recd.size(); i++) {
        if((recd[i]).belongTo(pt, eta)) return recd[i].effi+recd[i].effiErrorHigh;
        if((recd[i]).belongTo(190, eta)) hiPtBin = recd[i].effi;
    }
    return hiPtBin;
}

double SmearJetPt(double recoPt, double genPt, double eta, int smearJet, int year, int jetType){

    double centralSF(1.);
    double upSF(1.);
    double downSF(1.);

    // AK4 jets -------------------
    if (jetType == 0){
        if (year == 2016){
            // 2016 MC, AK4 Jets: Summer16_25nsV1_MC_SF_AK4PFchs
            // these SFs are symmetric in eta

            // centralSF -----
            if      (fabs(eta) < 0.522) centralSF = 1.1595;
            else if (fabs(eta) < 0.783) centralSF = 1.1948;
            else if (fabs(eta) < 1.131) centralSF = 1.1464;
            else if (fabs(eta) < 1.305) centralSF = 1.1609;
            else if (fabs(eta) < 1.740) centralSF = 1.1278;
            else if (fabs(eta) < 1.930) centralSF = 1.1000;
            else if (fabs(eta) < 2.043) centralSF = 1.1426;
            else if (fabs(eta) < 2.322) centralSF = 1.1512;
            else if (fabs(eta) < 2.500) centralSF = 1.2963;
            else if (fabs(eta) < 2.853) centralSF = 1.3418;
            else if (fabs(eta) < 2.964) centralSF = 1.7788;
            else if (fabs(eta) < 3.139) centralSF = 1.1869;
            else if (fabs(eta) < 5.191) centralSF = 1.1922;
            else centralSF = 1.1922;
            
            // upSF -----
            if      (fabs(eta) < 0.522) upSF = 1.224;
            else if (fabs(eta) < 0.783) upSF = 1.26;
            else if (fabs(eta) < 1.131) upSF = 1.2096;
            else if (fabs(eta) < 1.305) upSF = 1.2634;
            else if (fabs(eta) < 1.740) upSF = 1.2264;
            else if (fabs(eta) < 1.930) upSF = 1.2079;
            else if (fabs(eta) < 2.043) upSF = 1.264;
            else if (fabs(eta) < 2.322) upSF = 1.2652;
            else if (fabs(eta) < 2.500) upSF = 1.5334;
            else if (fabs(eta) < 2.853) upSF = 1.5509;
            else if (fabs(eta) < 2.964) upSF = 1.9796;
            else if (fabs(eta) < 3.139) upSF = 1.3112;
            else if (fabs(eta) < 5.191) upSF = 1.341;
            else upSF = 1.341;

            // downSF -----
            if      (fabs(eta) < 0.522) downSF = 1.095;
            else if (fabs(eta) < 0.783) downSF = 1.1296;
            else if (fabs(eta) < 1.131) downSF = 1.0832;
            else if (fabs(eta) < 1.305) downSF = 1.0584;
            else if (fabs(eta) < 1.740) downSF = 1.0292;
            else if (fabs(eta) < 1.930) downSF = 0.9921;
            else if (fabs(eta) < 2.043) downSF = 1.0212;
            else if (fabs(eta) < 2.322) downSF = 1.0372;
            else if (fabs(eta) < 2.500) downSF = 1.0592;
            else if (fabs(eta) < 2.853) downSF = 1.1327;
            else if (fabs(eta) < 2.964) downSF = 1.578;
            else if (fabs(eta) < 3.139) downSF = 1.0626;
            else if (fabs(eta) < 5.191) downSF = 1.0434;
            else downSF = 1.0434;

        }   
        else if (year == 2017){
            // 2017 MC, AK4 Jets: Fall17_V3_MC_SF_AK4PFchs
            // these SFs are symmetric in eta

            // centralSF -----
            if      (fabs(eta) < 0.522) centralSF = 1.1432;
            else if (fabs(eta) < 0.783) centralSF = 1.1815;
            else if (fabs(eta) < 1.131) centralSF = 1.0989;
            else if (fabs(eta) < 1.305) centralSF = 1.1137;
            else if (fabs(eta) < 1.740) centralSF = 1.1307;
            else if (fabs(eta) < 1.930) centralSF = 1.1600;
            else if (fabs(eta) < 2.043) centralSF = 1.2393;
            else if (fabs(eta) < 2.322) centralSF = 1.2604;
            else if (fabs(eta) < 2.500) centralSF = 1.4085;
            else if (fabs(eta) < 2.853) centralSF = 1.9909;
            else if (fabs(eta) < 2.964) centralSF = 2.2923;
            else if (fabs(eta) < 3.139) centralSF = 1.2696;
            else if (fabs(eta) < 5.191) centralSF = 1.1542;
            else centralSF = 1.1542;

            // upSF -----
            if      (fabs(eta) < 0.522) upSF = 1.1654;
            else if (fabs(eta) < 0.783) upSF = 1.2299;
            else if (fabs(eta) < 1.131) upSF = 1.1444;
            else if (fabs(eta) < 1.305) upSF = 1.2533;
            else if (fabs(eta) < 1.740) upSF = 1.2778;
            else if (fabs(eta) < 1.930) upSF = 1.2576;
            else if (fabs(eta) < 2.043) upSF = 1.4301;
            else if (fabs(eta) < 2.322) upSF = 1.4105;
            else if (fabs(eta) < 2.500) upSF = 1.6105;
            else if (fabs(eta) < 2.853) upSF = 2.5593;
            else if (fabs(eta) < 2.964) upSF = 2.6665;
            else if (fabs(eta) < 3.139) upSF = 1.3785;
            else if (fabs(eta) < 5.191) upSF = 1.3066;
            else upSF = 1.3066;

            // downSF -----
            if      (fabs(eta) < 0.522) downSF = 1.1210;
            else if (fabs(eta) < 0.783) downSF = 1.1332;
            else if (fabs(eta) < 1.131) downSF = 1.0533;
            else if (fabs(eta) < 1.305) downSF = 0.9740;
            else if (fabs(eta) < 1.740) downSF = 0.9837;
            else if (fabs(eta) < 1.930) downSF = 1.0623;
            else if (fabs(eta) < 2.043) downSF = 1.0484;
            else if (fabs(eta) < 2.322) downSF = 1.1103;
            else if (fabs(eta) < 2.500) downSF = 1.2066;
            else if (fabs(eta) < 2.853) downSF = 1.4225;
            else if (fabs(eta) < 2.964) downSF = 1.9180;
            else if (fabs(eta) < 3.139) downSF = 1.1607;
            else if (fabs(eta) < 5.191) downSF = 1.0019;
            else downSF = 1.0019;

        }
    }

    // AK8 jets -------------------
    if (jetType == 1){
        if (year == 2016){
            // 2016 MC, AK8 Jets: Summer16_25nsV1_MC_SF_AK8PFPuppi
            // these SFs are symmetric in eta

            // centralSF -----
            if      (fabs(eta) < 0.522) centralSF = 1.1595;
            else if (fabs(eta) < 0.783) centralSF = 1.1948;
            else if (fabs(eta) < 1.131) centralSF = 1.1464;
            else if (fabs(eta) < 1.305) centralSF = 1.1609;
            else if (fabs(eta) < 1.740) centralSF = 1.1278;
            else if (fabs(eta) < 1.930) centralSF = 1.1000;
            else if (fabs(eta) < 2.043) centralSF = 1.1426;
            else if (fabs(eta) < 2.322) centralSF = 1.1512;
            else if (fabs(eta) < 2.500) centralSF = 1.2963;
            else if (fabs(eta) < 2.853) centralSF = 1.3418;
            else if (fabs(eta) < 2.964) centralSF = 1.7788;
            else if (fabs(eta) < 3.139) centralSF = 1.1869;
            else if (fabs(eta) < 5.191) centralSF = 1.1922;
            else centralSF = 1.1922;
            
            // upSF -----
            if      (fabs(eta) < 0.522) upSF = 1.224;
            else if (fabs(eta) < 0.783) upSF = 1.26;
            else if (fabs(eta) < 1.131) upSF = 1.2096;
            else if (fabs(eta) < 1.305) upSF = 1.2634;
            else if (fabs(eta) < 1.740) upSF = 1.2264;
            else if (fabs(eta) < 1.930) upSF = 1.2079;
            else if (fabs(eta) < 2.043) upSF = 1.264;
            else if (fabs(eta) < 2.322) upSF = 1.2652;
            else if (fabs(eta) < 2.500) upSF = 1.5334;
            else if (fabs(eta) < 2.853) upSF = 1.5509;
            else if (fabs(eta) < 2.964) upSF = 1.9796;
            else if (fabs(eta) < 3.139) upSF = 1.3112;
            else if (fabs(eta) < 5.191) upSF = 1.341;
            else upSF = 1.341;

            // downSF -----
            if      (fabs(eta) < 0.522) downSF = 1.095;
            else if (fabs(eta) < 0.783) downSF = 1.1296;
            else if (fabs(eta) < 1.131) downSF = 1.0832;
            else if (fabs(eta) < 1.305) downSF = 1.0584;
            else if (fabs(eta) < 1.740) downSF = 1.0292;
            else if (fabs(eta) < 1.930) downSF = 0.9921;
            else if (fabs(eta) < 2.043) downSF = 1.0212;
            else if (fabs(eta) < 2.322) downSF = 1.0372;
            else if (fabs(eta) < 2.500) downSF = 1.0592;
            else if (fabs(eta) < 2.853) downSF = 1.1327;
            else if (fabs(eta) < 2.964) downSF = 1.578;
            else if (fabs(eta) < 3.139) downSF = 1.0626;
            else if (fabs(eta) < 5.191) downSF = 1.0434;
            else downSF = 1.0434;

        }
        else if (year == 2017){
            // 2017 MC, AK8 Jets: Fall17_V3_MC_SF_AK8PFPuppi
            // these SFs are symmetric in eta

            // centralSF -----
            if      (fabs(eta) < 0.522) centralSF = 1.1432;
            else if (fabs(eta) < 0.783) centralSF = 1.1815;
            else if (fabs(eta) < 1.131) centralSF = 1.0989;
            else if (fabs(eta) < 1.305) centralSF = 1.1137;
            else if (fabs(eta) < 1.740) centralSF = 1.1307;
            else if (fabs(eta) < 1.930) centralSF = 1.1600;
            else if (fabs(eta) < 2.043) centralSF = 1.2393;
            else if (fabs(eta) < 2.322) centralSF = 1.2604;
            else if (fabs(eta) < 2.500) centralSF = 1.4085;
            else if (fabs(eta) < 2.853) centralSF = 1.9909;
            else if (fabs(eta) < 2.964) centralSF = 2.2923;
            else if (fabs(eta) < 3.139) centralSF = 1.2696;
            else if (fabs(eta) < 5.191) centralSF = 1.1542;
            else centralSF = 1.1542;

            // upSF -----
            if      (fabs(eta) < 0.522) upSF = 1.1654;
            else if (fabs(eta) < 0.783) upSF = 1.2299;
            else if (fabs(eta) < 1.131) upSF = 1.1444;
            else if (fabs(eta) < 1.305) upSF = 1.2533;
            else if (fabs(eta) < 1.740) upSF = 1.2778;
            else if (fabs(eta) < 1.930) upSF = 1.2576;
            else if (fabs(eta) < 2.043) upSF = 1.4301;
            else if (fabs(eta) < 2.322) upSF = 1.4105;
            else if (fabs(eta) < 2.500) upSF = 1.6105;
            else if (fabs(eta) < 2.853) upSF = 2.5593;
            else if (fabs(eta) < 2.964) upSF = 2.6665;
            else if (fabs(eta) < 3.139) upSF = 1.3785;
            else if (fabs(eta) < 5.191) upSF = 1.3066;
            else upSF = 1.3066;
            
            // downSF -----
            if      (fabs(eta) < 0.522) downSF = 1.1210;
            else if (fabs(eta) < 0.783) downSF = 1.1332;
            else if (fabs(eta) < 1.131) downSF = 1.0533;
            else if (fabs(eta) < 1.305) downSF = 0.9740;
            else if (fabs(eta) < 1.740) downSF = 0.9837;
            else if (fabs(eta) < 1.930) downSF = 1.0623;
            else if (fabs(eta) < 2.043) downSF = 1.0484;
            else if (fabs(eta) < 2.322) downSF = 1.1103;
            else if (fabs(eta) < 2.500) downSF = 1.2066;
            else if (fabs(eta) < 2.853) downSF = 1.4225;
            else if (fabs(eta) < 2.964) downSF = 1.9180;
            else if (fabs(eta) < 3.139) downSF = 1.1607;
            else if (fabs(eta) < 5.191) downSF = 1.0019;
            else downSF = 1.0019;

        }
    }
    
    double smearedPt(0.);
    // central SF
    if (smearJet == 0)       smearedPt = std::max(0., genPt + centralSF*(recoPt - genPt));
    // SF variation up
    else if (smearJet == 1)  smearedPt = std::max(0., genPt + upSF*(recoPt - genPt));
    // SF variation down
    else if (smearJet == -1) smearedPt = std::max(0., genPt + downSF*(recoPt - genPt));
    
    return smearedPt;
}

double SmearJetPtLite(double recoPt, double genPt, double scaleFactor){

    // used for --
    // 2018 MC, AK4 Jets: Autumn18_V7b_MC_SF_AK4PFchs
    // 2018 MC, AK8 Jets: Autumn18_V7b_MC_SF_AK8PFPuppi

    double smearedPt = std::max(0., genPt + scaleFactor*(recoPt - genPt));

    return smearedPt;
}

void normalizeTH2D(TH2D *h)
{
    int xbin(h->GetNbinsX()), ybin(h->GetNbinsY());
    for (int i(1); i <= ybin; i++){
        double sum(0.);
        for (int j(1); j <= xbin; j++){
            sum += h->GetBinContent(j,i);
        }
        for (int j(1); j <= xbin; j++){
            if (sum > 0) h->SetBinContent(j, i, h->GetBinContent(j, i) / sum );
        }
    }
}

void bestTwoJetsCandidatesPt(vector<jetStruct> jets, pair<TLorentzVector, TLorentzVector>& bestTwoJets)
{
    int nGoodJets(jets.size());
    if (nGoodJets >= 2){
        //cout << "\nMore than 2 jets, selecting best pair" << endl;
        double minPt(999999.);
        for (int i(0); i < nGoodJets - 1; i++) {
            TLorentzVector jeti; jeti.SetPtEtaPhiE(jets[i].pt, jets[i].eta, jets[i].phi, jets[i].energy);
            for (int j(i + 1); j < nGoodJets; j++) {
                TLorentzVector jetj; jetj.SetPtEtaPhiE(jets[j].pt, jets[j].eta, jets[j].phi, jets[j].energy);
                TLorentzVector jetij = jeti + jetj;
                //cout << i << " " << j << ": Pair pt = " << jetij.Pt() << endl;
                if (jetij.Pt() < minPt){
                    bestTwoJets.first = jeti; 
                    bestTwoJets.second = jetj;
                    minPt = jetij.Pt();
                    //cout << "Smallest pt = " << jetij.Pt() << endl;
                }
            }
        } 
    }
}

void bestTwoJetsCandidatesPhi(vector<jetStruct> jets, pair<TLorentzVector, TLorentzVector>& bestTwoJets)
{
    int nGoodJets(jets.size());
    if (nGoodJets >= 2){
        //cout << "\nMore than 2 jets, selecting best pair" << endl;
        double maxdPhi(-0.0001);
        for (int i(0); i < nGoodJets - 1; i++) {
            TLorentzVector jeti; jeti.SetPtEtaPhiE(jets[i].pt, jets[i].eta, jets[i].phi, jets[i].energy);
            for (int j(i + 1); j < nGoodJets; j++) {
                TLorentzVector jetj; jetj.SetPtEtaPhiE(jets[j].pt, jets[j].eta, jets[j].phi, jets[j].energy);
                double dPhi = deltaPhi(jeti, jetj);
                //cout << i << " " << j << ": dPhi = " << dPhi << endl;
                if (dPhi > maxdPhi){
                    bestTwoJets.first = jeti; 
                    bestTwoJets.second = jetj;
                    maxdPhi = dPhi;
                    //cout << "Biggest dPhi = " << dPhi << endl;
                }
            }
        } 
    }
}

vector<double> buildVecFineBin( int nStdBin, double arrStdBin[], int factChop)
{
    vector<double> vecTemp;
    for (int i = 0; i < nStdBin; i++){
        double binWidth = (arrStdBin[i+1] - arrStdBin[i])/5;
        for (int j = 0; j < factChop; j++){
            double element(0.);
            element = arrStdBin[i] + (j * binWidth);
            vecTemp.push_back(element);
        }
    }
    vecTemp.push_back(arrStdBin[nStdBin]);
    return vecTemp;
}

void splitBinsInTwoForTUnfold(int nBinsOriginal, double originalArray[], double splitArray[]){

    int j(0);
    for(int i(0); i < nBinsOriginal; i++){
        splitArray[j] = originalArray[i];
        splitArray[j+1] = originalArray[i] + (originalArray[i+1] - originalArray[i])/2.;
        if (i == nBinsOriginal - 1) splitArray[j+2] = originalArray[i+1];
        j += 2;
    }
    
}

void welcomeMessage(){
    std::cout << R"(

                                                  __    __       __       _                           
                                                 / / /\ \ \ _    \ \  ___| |_ ___                     
                                                 \ \/  \/ /| |_   \ \/ _ \ __/ __|                    
                                                  \  /\  /_   _/\_/ /  __/ |_\__ \                    
                                                   \/  \/  |_| \___/ \___|\__|___/                    
                                    __                 _       __      _           _   _              
                                   /__\_   _____ _ __ | |_    / _\ ___| | ___  ___| |_(_) ___  _ __   
                                  /_\ \ \ / / _ \ '_ \| __|   \ \ / _ \ |/ _ \/ __| __| |/ _ \| '_ \  
                                 //__  \ V /  __/ | | | |_    _\ \  __/ |  __/ (__| |_| | (_) | | | | 
                                 \__/   \_/ \___|_| |_|\__|   \__/\___|_|\___|\___|\__|_|\___/|_| |_| 
                             	
                        
    )" << '\n'; 
}

void bTagVetoMessage(int doBJets){
    if (doBJets > -1){
        std::cout << "\n                       ///////////////////////////////////////////////////////////////////////////////////////" << std::endl;
        std::cout << "                       //                                                                                   //" << std::endl;
        std::cout << "                       //                              ---   WARNING   ---                                  //" << std::endl;
        std::cout << "                       //               You are running event selection without a b-tag veto.               //" << std::endl;
        std::cout << "                       //                       Are you sure you want to do this?                           //" << std::endl;
        std::cout << "                       //                                                                                   //" << std::endl;
        std::cout << "                       ///////////////////////////////////////////////////////////////////////////////////////\n" << std::endl;
    }
}



