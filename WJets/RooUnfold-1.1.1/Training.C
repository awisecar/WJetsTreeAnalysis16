#define Training_cxx
#include "Training.h"
#include "RooUnfoldResponse.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
//#include <TSystem.h>

  TH1D* prettyHistogram(std::string name, std::string title, std::string xtitle, double nbins, double rlow, double rhigh ){
    TH1D* hist = new TH1D(name.c_str(),title.c_str(),nbins,rlow,rhigh);
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    hist->Sumw2();
    return hist;
  }


TH1D* prettyHistogram(std::string name, std::string title, std::string xtitle, double nbins, float* xbins){
    TH1D* hist = new TH1D(name.c_str(),title.c_str(),nbins,xbins);
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    hist->Sumw2();
    return hist;
  }

double Phi_0_2pi(double x) {

  while (x >= 2*M_PI) x -= 2*M_PI;
  while (x <     0.)  x += 2*M_PI;
  return x;
}
double radius(double eta1, double phi1,double eta2, double phi2){

  const double TWOPI= 2.0*M_PI;
  phi1=Phi_0_2pi(phi1);
  phi2=Phi_0_2pi(phi2);
  
  double dphi=Phi_0_2pi(phi1-phi2);
  dphi = TMath::Min(dphi,TWOPI-dphi);
  double deta = eta1-eta2;
  
  return sqrt(deta*deta+dphi*dphi);
}
bool DescendingOrder(const TLorentzVector& l1, const TLorentzVector& l2){
  return l1.Pt()>l2.Pt();
}   
    



void Training::Loop()
{
//#ifdef __CINT__
//gSystem->Load("/home/mekter/ZJet/RooUnfold-1.1.1/libRooUnfold");
//#endif
//gSystem->AddIncludePath(" -I/home/mekter/RooUnfold-1.1.1/src ");

  if (fChain == 0) return;
//  TFile* fs=new TFile("testJES.root", "recreate");

//  float jpt1[18] = {30,40, 50, 60, 70, 80, 100, 120, 140, 160, 180, 210, 240, 270, 300, 350, 400, 450};
  float jptG1[19] = {15, 30,40, 50, 60, 70, 80, 100, 120, 140, 160, 180, 210, 240, 270, 300, 350, 400, 450};
  float jptP1[17] = {30,40, 50, 60, 70, 80, 100, 120, 140, 160, 180, 210, 240, 270, 300, 350, 400};
  float jpt2[11] = {15, 30,50,70,90,120,150,200,270,400, 450};
  float jpt3[8] = {15,30,50,70,90,120,400,450};
  float jpt4[8] = {15,30,50,70,90,120,400,450};

  float jptTest[10] = {20,30,40,55,75,105,150,210,315,500};
  TH1D* testJESerr_nominal = prettyHistogram("testJESerr_nominal","First Leading PatJet (Z+PatJets)","p_{T} [GeV]",9,jptTest);
  TH1D* testJESerr_scaled = prettyHistogram("testJESerr_scaled","First Leading PatJet (Z+PatJets)","p_{T} [GeV]",9,jptTest);
  TH1D* ratio = prettyHistogram("ratio","First Leading PatJet (Z+PatJets)","p_{T} [GeV]",9,jptTest);

//TRAINING  
  TH1F* PatJetPt_uniform = new TH1F("PatJetPt_1_uniform","",300,30,330);
//reconstructed (measured) training distribution
  TH1D* PatJetPt_1=prettyHistogram("PatJetPt_1","First Leading PatJet (Z+PatJets)","p_{T} [GeV]",18,jptG1);
  TH1D* PatJetPt_2=prettyHistogram("PatJetPt_2","Second Leading PatJet (Z+PatJets)","p_{T}(PatJet) [GeV]",10,jpt2);
  TH1D* PatJetPt_3=prettyHistogram("PatJetPt_3","Third Leading (Z+PatJets)","p_{T}(PatJet) [GeV]",7,jpt3);
  TH1D* PatJetPt_4=prettyHistogram("PatJetPt_4","Fourth Leading PatJet (Z+PatJets)","p_{T}(PatJet) [GeV]",7,jpt4);
// generated (true) training distribution
  TH1D* GenJetPt_1=prettyHistogram("GenJetPt_1","First Leading GenJet (Z+GenJets)","p_{T} [GeV]",18,jptG1);
  TH1D* GenJetPt_2=prettyHistogram("GenJetPt_2","Second Leading GenJet (Z+GenJets)","p_{T}(GenJet) [GeV]",10,jpt2);
  TH1D* GenJetPt_3=prettyHistogram("GenJetPt_3","Third Leading (Z+GenJets)","p_{T}(GenJet) [GeV]",7,jpt3);
  TH1D* GenJetPt_4=prettyHistogram("GenJetPt_4","Fourth Leading GenJet (Z+GenJets)","p_{T}(GenJet) [GeV]",7,jpt4);
//reconstructed vs generated 
   TH2D *JetPt_1_Pat_Gen = new TH2D("JetPt_1_Pat_Gen","First Leading Reconstructed vs generated jet pt",18,jptG1,18,jptG1);
JetPt_1_Pat_Gen->Sumw2();
   TH2D *JetPt_2_Pat_Gen = new TH2D("JetPt_2_Pat_Gen","Second Leading Reconstructed vs generated jet pt",10,jpt2,10,jpt2);
JetPt_2_Pat_Gen->Sumw2();
   TH2D *JetPt_3_Pat_Gen = new TH2D("JetPt_3_Pat_Gen","Third Leading Reconstructed vs generated jet pt",7,jpt3,7,jpt3);
JetPt_3_Pat_Gen->Sumw2();
   TH2D *JetPt_4_Pat_Gen = new TH2D("JetPt_4_Pat_Gen","Forth Leading Reconstructed vs generated jet pt",7,jpt4,7,jpt4);
JetPt_4_Pat_Gen->Sumw2(); 

//TEST
//reconstructed (measured) test distribution
  TH1D* PatJetPt_1_test=prettyHistogram("PatJetPt_1_test","First Leading PatJet (Z+PatJets)","p_{T} [GeV]",18,jptG1);
  TH1D* PatJetPt_2_test=prettyHistogram("PatJetPt_2_test","Second Leading PatJet (Z+PatJets)","p_{T}(PatJet) [GeV]",10,jpt2);
  TH1D* PatJetPt_3_test=prettyHistogram("PatJetPt_3_test","Third Leading (Z+PatJets)","p_{T}(PatJet) [GeV]",7,jpt3);
  TH1D* PatJetPt_4_test=prettyHistogram("PatJetPt_4_test","Fourth Leading PatJet (Z+PatJets)","p_{T}(PatJet) [GeV]",7,jpt4);
// generated (true) test distribution
  TH1D* GenJetPt_1_test=prettyHistogram("GenJetPt_1_test","First Leading GenJet (Z+GenJets)","p_{T} [GeV]",18,jptG1);
  TH1D* GenJetPt_2_test=prettyHistogram("GenJetPt_2_test","Second Leading GenJet (Z+GenJets)","p_{T}(GenJet) [GeV]",10,jpt2);
  TH1D* GenJetPt_3_test=prettyHistogram("GenJetPt_3_test","Third Leading (Z+GenJets)","p_{T}(GenJet) [GeV]",7,jpt3);
  TH1D* GenJetPt_4_test=prettyHistogram("GenJetPt_4_test","Fourth Leading GenJet (Z+GenJets)","p_{T}(GenJet) [GeV]",7,jpt4);


RooUnfoldResponse *response = new RooUnfoldResponse(PatJetPt_1, GenJetPt_1);


  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
//  for (Long64_t jentry=0; jentry<nentries;jentry++) {


//  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  for (Long64_t jentry=0; jentry<3500000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
  if (jentry%500000 == 0) std::cout<<jentry<<std::endl;
 
    //Require Events to pass full offline selection chain.
    //if((patZJet_pxGenZ ==-9999)||(patZJet_pxPatZ ==-9999))continue;
    if((patZJet_pxD1GenZ==-9999)||(patZJet_pxD2GenZ==-9999))continue;
    //if((patZJet_pxD2PatZ==-9999)||(patZJet_pxD2PatZ==-9999))continue;
    //if(!(patZJet_passD1PatZ && patZJet_passD2PatZ))continue;
    



    TLorentzVector *ZGen = new TLorentzVector();
    TLorentzVector *ZPat = new TLorentzVector();
    TLorentzVector *gEl1 = new TLorentzVector();    
    TLorentzVector *gEl2 = new TLorentzVector();    
    TLorentzVector *pEl1 = new TLorentzVector();    
    TLorentzVector *pEl2 = new TLorentzVector();    
   
  
    ZGen->SetPxPyPzE(patZJet_pxGenZ,patZJet_pyGenZ,patZJet_pzGenZ,patZJet_enGenZ);
    ZPat->SetPxPyPzE(patZJet_pxPatZ,patZJet_pyPatZ,patZJet_pzPatZ,patZJet_enPatZ);
    gEl1->SetPxPyPzE(patZJet_pxD1GenZ,patZJet_pyD1GenZ,patZJet_pzD1GenZ,patZJet_enD1GenZ);
    gEl2->SetPxPyPzE(patZJet_pxD2GenZ,patZJet_pyD2GenZ,patZJet_pzD2GenZ,patZJet_enD2GenZ);
    pEl1->SetPxPyPzE(patZJet_pxD1PatZ,patZJet_pyD1PatZ,patZJet_pzD1PatZ,patZJet_enD1PatZ);
    pEl2->SetPxPyPzE(patZJet_pxD2PatZ,patZJet_pyD2PatZ,patZJet_pzD2PatZ,patZJet_enD2PatZ);
    
    TLorentzVector z(*gEl1+*gEl2);
    if(!(z.M()>70 && z.M()<110))continue;
/*
    //particle level cuts 
    if(!( fabs(gEl1->Eta())<=2.4 && fabs(gEl2->Eta())<=2.4))continue;
    if( fabs(gEl1->Eta())>1.44 && fabs(gEl1->Eta())<1.56)continue;
    if( fabs(gEl2->Eta())>1.44 && fabs(gEl2->Eta())<1.56)continue;
    if(!(gEl1->Pt()>20 && gEl2->Pt()>20))continue;

     //detector level cuts 
    if(!(fabs(pEl1->Eta())<2.4 &&fabs(pEl2->Eta())<2.4))continue;
    if( fabs(pEl1->Eta())>1.44 && fabs(pEl1->Eta())<1.56)continue;
    if( fabs(pEl2->Eta())>1.44 && fabs(pEl2->Eta())<1.56)continue;
    if(!(pEl1->Pt()>20 && pEl2->Pt()>20))continue;

    TLorentzVector zr(*pEl1+*pEl2);
    if(!(zr.M()>70 && zr.M()<110))continue;
    if(!(patZJet_passD1PatZ && patZJet_passD2PatZ))continue;


    std::vector<TLorentzVector> gJets;
    for(int i=0; i<patZJet_genJetPx->size(); i++){
      TLorentzVector gJ(patZJet_genJetPx->at(i),
			    patZJet_genJetPy->at(i),
			    patZJet_genJetPz->at(i),
			    patZJet_genJetEn->at(i)
			    );
      

      if (gJ.Pt()<15 || gJ.Pt()>450)continue;
      if(fabs(gJ.Eta())> 2.4)continue;
      double dRe1 = radius(gJ.Eta(),gJ.Phi(),gEl1->Eta(),gEl1->Phi());
      double dRe2 = radius(gJ.Eta(),gJ.Phi(),gEl2->Eta(),gEl2->Phi());
      if((dRe1<0.3)||(dRe2<0.3))continue;
      gJets.push_back(gJ);
    }
 

    std::vector<TLorentzVector> pJets;
    for(int i=0; i<patZJet_patJetPx->size(); i++){
      TLorentzVector pJ(
			 patZJet_patJetPx->at(i),
			 patZJet_patJetPy->at(i),
			 patZJet_patJetPz->at(i),
			 patZJet_patJetEn->at(i)
			 );
      
      if (pJ.Pt()<15 || pJ.P()>450)continue;
      if(fabs(pJ.Eta())> 2.4)continue;
      double dRe1 = radius(pJ.Eta(),pJ.Phi(),pEl1->Eta(),pEl1->Phi());
      double dRe2 = radius(pJ.Eta(),pJ.Phi(),pEl2->Eta(),pEl2->Phi());
      if((dRe1<0.3)||(dRe2<0.3))continue;
      pJets.push_back(pJ);
    }

    //sort both the collections according to Pt
    
    std::sort(gJets.begin(),gJets.end(),DescendingOrder);
    std::sort(pJets.begin(),pJets.end(),DescendingOrder);


// match reconstructed jets to generated

 int matchedGenIndex[10][10] = {{0}};
if (gJets.size()>0 && pJets.size()>0){

for (int l=0; l<pJets.size(); l++)
    {
  for(int k=0; k< gJets.size(); k++)
    {
    double dR=1;
    dR = radius(gJets[k].Eta(),gJets[k].Phi(),pJets[l].Eta(),pJets[l].Phi()); 
      if(dR<0.2)
	{
	  matchedGenIndex[l][k] = 1;
	}
    } // end of loop over generator-level jets

}// end of loop over reco level jets

}


    double weight = 1; 
*/


//Fill training
//if (jentry<nentries){

   //first jet
/*
  if(gJets.size()>0){
     for(int l=0; l< gJets.size(); l++)
        {
         if (matchedGenIndex[0][l] == 1)
	    {
            GenJetPt_1->Fill(gJets[l].Pt(),weight);
PatJetPt_uniform->Fill(gJets[l].Pt(),weight);

testJESerr_nominal->Fill(gJets[l].Pt(),weight);
testJESerr_scaled->Fill(1.01*gJets[l].Pt(),weight);
             if (pJets[0].Pt()>30){
            JetPt_1_Pat_Gen->Fill(pJets[0].Pt(),gJets[l].Pt(),weight);
            PatJetPt_1->Fill(pJets[0].Pt(),weight);
                }
	    }
        } // end of loop over generator-level jet
   }

   //second jet

   if(gJets.size()>0){
     for(int l=0; l< gJets.size(); l++)
        {
        if (matchedGenIndex[1][l] == 1)
	   {
           GenJetPt_2->Fill(gJets[l].Pt(),weight);

             if (pJets[1].Pt()>30){
           JetPt_2_Pat_Gen->Fill(pJets[1].Pt(),gJets[l].Pt(),weight);
           PatJetPt_2->Fill(pJets[1].Pt(),weight);
                 }
           }
        }
    }
*/


/*
   //third jet
   if(gJets.size()>0 && pJets.size()>2){
     for(int k=0; k< gJets.size(); k++)
        {
        if (matchedGenIndex[2][k] == 1)
	   {
           GenJetPt_3->Fill(gJets[k].Pt(),weight);
           PatJetPt_3->Fill(pJets[2].Pt(),weight);
           JetPt_3_Pat_Gen->Fill(pJets[2].Pt(),gJets[k].Pt(),weight);
           }
        }
     }

   //forth jet
   if(gJets.size()>0 && pJets.size()>3){
     for(int k=0; k< gJets.size(); k++)
        {
        if (matchedGenIndex[3][k] == 1)
	   {
           GenJetPt_4->Fill(gJets[k].Pt(),weight);
           PatJetPt_4->Fill(pJets[3].Pt(),weight);
           JetPt_4_Pat_Gen->Fill(pJets[3].Pt(),gJets[k].Pt(),weight);
           }
        }
     }
*/

//}//jentry%2==0


//Fill test
//if (jentry%2==1){

/*
   //first jet
 
   if(gJets.size()>0){
     for(int l=0; l< gJets.size(); l++)
        {
        if (matchedGenIndex[0][l] == 1)
	   {
           GenJetPt_1_test->Fill(gJets[l].Pt(),weight);
             if (pJets[0].Pt()>30){
           PatJetPt_1_test->Fill(pJets[0].Pt(),weight);
               }
	   }
       } // end of loop over generator-level jet
    }

   //second jet

   if(gJets.size()>0){
     for( int l=0; l< gJets.size(); l++)
        {
        if (matchedGenIndex[1][l] == 1)
	   {
          GenJetPt_2_test->Fill(gJets[l].Pt(),weight); 
             if (pJets[1].Pt()>30){
          PatJetPt_2_test->Fill(pJets[1].Pt(),weight);
              }
           }
        }
    }
*/

/*
   //third jet
   if(gJets.size()>0 && pJets.size()>2){
     for(int k=0; k< gJets.size(); k++)
        {
        if (matchedGenIndex[2][k] == 1)
	   {
           GenJetPt_3_test->Fill(gJets[k].Pt(),weight);
           PatJetPt_3_test->Fill(pJets[2].Pt(),weight);
           }
        }
     }
*/
/*
   //forth jet
   if(gJets.size()>0 && pJets.size()>3){
     for( int k=0; k< gJets.size(); k++)
        {
        if (matchedGenIndex[3][k] == 1)
	   {
           GenJetPt_4_test->Fill(gJets[k].Pt(),weight);
           PatJetPt_4_test->Fill(pJets[3].Pt(),weight);
           }
        }
     }
*/

//}//jentry%2==1




  }//jentry


 //fs->cd();
 //fs->Write();
 //fs->Close();

//ratio->Divide(testJESerr_scaled,testJESerr_nominal,1.,1.);
//ratio->Draw("hist");


}
