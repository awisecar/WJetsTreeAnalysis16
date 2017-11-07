//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 23 20:08:25 2012 by ROOT version 5.28/00c
// from TTree tree/tree
// found on file: zeeTest.root
//////////////////////////////////////////////////////////

#ifndef Training_h
#define Training_h
#include <fstream>
#include<sstream>
#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TLorentzVector.h"
class Training {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        patZJet_pxGenZ;
   Double_t        patZJet_pyGenZ;
   Double_t        patZJet_pzGenZ;
   Double_t        patZJet_enGenZ;
   Double_t        patZJet_pxD1GenZ;
   Double_t        patZJet_pyD1GenZ;
   Double_t        patZJet_pzD1GenZ;
   Double_t        patZJet_enD1GenZ;
   Double_t        patZJet_pxD2GenZ;
   Double_t        patZJet_pyD2GenZ;
   Double_t        patZJet_pzD2GenZ;
   Double_t        patZJet_enD2GenZ;
   Double_t        patZJet_pxPatZ;
   Double_t        patZJet_pyPatZ;
   Double_t        patZJet_pzPatZ;
   Double_t        patZJet_enPatZ;
   Double_t        patZJet_pxD1PatZ;
   Double_t        patZJet_pyD1PatZ;
   Double_t        patZJet_pzD1PatZ;
   Double_t        patZJet_enD1PatZ;
   Int_t           patZJet_passD1PatZ;
   Double_t        patZJet_pxD2PatZ;
   Double_t        patZJet_pyD2PatZ;
   Double_t        patZJet_pzD2PatZ;
   Double_t        patZJet_enD2PatZ;
   Int_t           patZJet_passD2PatZ;
   std::vector<double>  *patZJet_genJetPx;
   std::vector<double>  *patZJet_genJetPy;
   std::vector<double>  *patZJet_genJetPz;
   std::vector<double>  *patZJet_genJetEn;
   std::vector<double>  *patZJet_patJetPx;
   std::vector<double>  *patZJet_patJetPy;
   std::vector<double>  *patZJet_patJetPz;
   std::vector<double>  *patZJet_patJetEn;

   // List of branches
   TBranch        *b_patZJet_pxGenZ;   //!
   TBranch        *b_patZJet_pyGenZ;   //!
   TBranch        *b_patZJet_pzGenZ;   //!
   TBranch        *b_patZJet_enGenZ;   //!
   TBranch        *b_patZJet_pxD1GenZ;   //!
   TBranch        *b_patZJet_pyD1GenZ;   //!
   TBranch        *b_patZJet_pzD1GenZ;   //!
   TBranch        *b_patZJet_enD1GenZ;   //!
   TBranch        *b_patZJet_pxD2GenZ;   //!
   TBranch        *b_patZJet_pyD2GenZ;   //!
   TBranch        *b_patZJet_pzD2GenZ;   //!
   TBranch        *b_patZJet_enD2GenZ;   //!
   TBranch        *b_patZJet_pxPatZ;   //!
   TBranch        *b_patZJet_pyPatZ;   //!
   TBranch        *b_patZJet_pzPatZ;   //!
   TBranch        *b_patZJet_enPatZ;   //!
   TBranch        *b_patZJet_pxD1PatZ;   //!
   TBranch        *b_patZJet_pyD1PatZ;   //!
   TBranch        *b_patZJet_pzD1PatZ;   //!
   TBranch        *b_patZJet_enD1PatZ;   //!
   TBranch        *b_patZJet_passD1PatZ;   //!
   TBranch        *b_patZJet_pxD2PatZ;   //!
   TBranch        *b_patZJet_pyD2PatZ;   //!
   TBranch        *b_patZJet_pzD2PatZ;   //!
   TBranch        *b_patZJet_enD2PatZ;   //!
   TBranch        *b_patZJet_passD2PatZ;   //!
   TBranch        *b_patZJet_genJetPx;   //!
   TBranch        *b_patZJet_genJetPy;   //!
   TBranch        *b_patZJet_genJetPz;   //!
   TBranch        *b_patZJet_genJetEn;   //!
   TBranch        *b_patZJet_patJetPx;   //!
   TBranch        *b_patZJet_patJetPy;   //!
   TBranch        *b_patZJet_patJetPz;   //!
   TBranch        *b_patZJet_patJetEn;   //!

   Training(TTree *tree=0);
   virtual ~Training();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Training_cxx
Training::Training(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/mekter/ZJet/UnfoldTree.root");
      if (!f) {
         f = new TFile("/home/mekter/ZJet/UnfoldTree.root");
         f->cd("/home/mekter/ZJet/UnfoldTree.root:/tree");
      }
      tree = (TTree*)gDirectory->Get("tree");

   }
   Init(tree);
}

Training::~Training()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Training::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Training::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Training::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   patZJet_genJetPx = 0;
   patZJet_genJetPy = 0;
   patZJet_genJetPz = 0;
   patZJet_genJetEn = 0;
   patZJet_patJetPx = 0;
   patZJet_patJetPy = 0;
   patZJet_patJetPz = 0;
   patZJet_patJetEn = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("patZJet_pxGenZ", &patZJet_pxGenZ, &b_patZJet_pxGenZ);
   fChain->SetBranchAddress("patZJet_pyGenZ", &patZJet_pyGenZ, &b_patZJet_pyGenZ);
   fChain->SetBranchAddress("patZJet_pzGenZ", &patZJet_pzGenZ, &b_patZJet_pzGenZ);
   fChain->SetBranchAddress("patZJet_enGenZ", &patZJet_enGenZ, &b_patZJet_enGenZ);
   fChain->SetBranchAddress("patZJet_pxD1GenZ", &patZJet_pxD1GenZ, &b_patZJet_pxD1GenZ);
   fChain->SetBranchAddress("patZJet_pyD1GenZ", &patZJet_pyD1GenZ, &b_patZJet_pyD1GenZ);
   fChain->SetBranchAddress("patZJet_pzD1GenZ", &patZJet_pzD1GenZ, &b_patZJet_pzD1GenZ);
   fChain->SetBranchAddress("patZJet_enD1GenZ", &patZJet_enD1GenZ, &b_patZJet_enD1GenZ);
   fChain->SetBranchAddress("patZJet_pxD2GenZ", &patZJet_pxD2GenZ, &b_patZJet_pxD2GenZ);
   fChain->SetBranchAddress("patZJet_pyD2GenZ", &patZJet_pyD2GenZ, &b_patZJet_pyD2GenZ);
   fChain->SetBranchAddress("patZJet_pzD2GenZ", &patZJet_pzD2GenZ, &b_patZJet_pzD2GenZ);
   fChain->SetBranchAddress("patZJet_enD2GenZ", &patZJet_enD2GenZ, &b_patZJet_enD2GenZ);
   fChain->SetBranchAddress("patZJet_pxPatZ", &patZJet_pxPatZ, &b_patZJet_pxPatZ);
   fChain->SetBranchAddress("patZJet_pyPatZ", &patZJet_pyPatZ, &b_patZJet_pyPatZ);
   fChain->SetBranchAddress("patZJet_pzPatZ", &patZJet_pzPatZ, &b_patZJet_pzPatZ);
   fChain->SetBranchAddress("patZJet_enPatZ", &patZJet_enPatZ, &b_patZJet_enPatZ);
   fChain->SetBranchAddress("patZJet_pxD1PatZ", &patZJet_pxD1PatZ, &b_patZJet_pxD1PatZ);
   fChain->SetBranchAddress("patZJet_pyD1PatZ", &patZJet_pyD1PatZ, &b_patZJet_pyD1PatZ);
   fChain->SetBranchAddress("patZJet_pzD1PatZ", &patZJet_pzD1PatZ, &b_patZJet_pzD1PatZ);
   fChain->SetBranchAddress("patZJet_enD1PatZ", &patZJet_enD1PatZ, &b_patZJet_enD1PatZ);
   fChain->SetBranchAddress("patZJet_passD1PatZ", &patZJet_passD1PatZ, &b_patZJet_passD1PatZ);
   fChain->SetBranchAddress("patZJet_pxD2PatZ", &patZJet_pxD2PatZ, &b_patZJet_pxD2PatZ);
   fChain->SetBranchAddress("patZJet_pyD2PatZ", &patZJet_pyD2PatZ, &b_patZJet_pyD2PatZ);
   fChain->SetBranchAddress("patZJet_pzD2PatZ", &patZJet_pzD2PatZ, &b_patZJet_pzD2PatZ);
   fChain->SetBranchAddress("patZJet_enD2PatZ", &patZJet_enD2PatZ, &b_patZJet_enD2PatZ);
   fChain->SetBranchAddress("patZJet_passD2PatZ", &patZJet_passD2PatZ, &b_patZJet_passD2PatZ);
   fChain->SetBranchAddress("patZJet_genJetPx", &patZJet_genJetPx, &b_patZJet_genJetPx);
   fChain->SetBranchAddress("patZJet_genJetPy", &patZJet_genJetPy, &b_patZJet_genJetPy);
   fChain->SetBranchAddress("patZJet_genJetPz", &patZJet_genJetPz, &b_patZJet_genJetPz);
   fChain->SetBranchAddress("patZJet_genJetEn", &patZJet_genJetEn, &b_patZJet_genJetEn);
   fChain->SetBranchAddress("patZJet_patJetPx", &patZJet_patJetPx, &b_patZJet_patJetPx);
   fChain->SetBranchAddress("patZJet_patJetPy", &patZJet_patJetPy, &b_patZJet_patJetPy);
   fChain->SetBranchAddress("patZJet_patJetPz", &patZJet_patJetPz, &b_patZJet_patJetPz);
   fChain->SetBranchAddress("patZJet_patJetEn", &patZJet_patJetEn, &b_patZJet_patJetEn);
   Notify();
}

Bool_t Training::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Training::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Training::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Training_cxx
