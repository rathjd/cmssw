//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 27 14:35:14 2019 by ROOT version 6.10/05
// from TTree eventTree/Event tree
// found on file: FinalOutput.root
//////////////////////////////////////////////////////////

#ifndef MuonTrackAnalyzer_h
#define MuonTrackAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class MuonTrackAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *trk_pt;
   vector<float>   *trk_eta;
   vector<float>   *trk_phi;
   vector<float>   *trk_charge;
   vector<float>   *trk_d0;
   vector<float>   *trk_z0;
   vector<float>   *trk_chi2;
   vector<int>     *trk_nstub;
   vector<int>     *trk_seed;
   vector<int>     *trk_genuine;
   vector<int>     *trk_loose;
   vector<int>     *trk_unknown;
   vector<int>     *trk_combinatoric;
   vector<int>     *trk_fake;
   vector<int>     *trk_matchtp_pdgid;
   vector<float>   *trk_matchtp_pt;
   vector<float>   *trk_matchtp_eta;
   vector<float>   *trk_matchtp_phi;
   vector<float>   *trk_matchtp_z0;
   vector<float>   *trk_matchtp_dxy;
   vector<int>     *trk_injet;
   vector<int>     *trk_injet_highpt;
   vector<int>     *trk_injet_vhighpt;
   vector<float>   *tp_pt;
   vector<float>   *tp_eta;
   vector<float>   *tp_phi;
   vector<float>   *tp_dxy;
   vector<float>   *tp_d0;
   vector<float>   *tp_z0;
   vector<float>   *tp_d0_prod;
   vector<float>   *tp_z0_prod;
   vector<int>     *tp_pdgid;
   vector<int>     *tp_nmatch;
   vector<int>     *tp_nloosematch;
   vector<int>     *tp_nstub;
   vector<int>     *tp_eventid;
   vector<int>     *tp_charge;
   vector<int>     *tp_injet;
   vector<int>     *tp_injet_highpt;
   vector<int>     *tp_injet_vhighpt;
   vector<float>   *matchtrk_pt;
   vector<float>   *matchtrk_eta;
   vector<float>   *matchtrk_phi;
   vector<float>   *matchtrk_z0;
   vector<float>   *matchtrk_d0;
   vector<float>   *matchtrk_chi2;
   vector<int>     *matchtrk_nstub;
   vector<int>     *matchtrk_seed;
   vector<int>     *matchtrk_injet;
   vector<int>     *matchtrk_injet_highpt;
   vector<int>     *matchtrk_injet_vhighpt;
   vector<bool>    *matchtrk_charge;
   vector<float>   *loosematchtrk_pt;
   vector<float>   *loosematchtrk_eta;
   vector<float>   *loosematchtrk_phi;
   vector<float>   *loosematchtrk_z0;
   vector<float>   *loosematchtrk_d0;
   vector<float>   *loosematchtrk_chi2;
   vector<int>     *loosematchtrk_nstub;
   vector<int>     *loosematchtrk_seed;
   vector<int>     *loosematchtrk_injet;
   vector<int>     *loosematchtrk_injet_highpt;
   vector<int>     *loosematchtrk_injet_vhighpt;
   vector<float>   *matchmuon_pt;
   vector<float>   *matchmuon_eta;
   vector<float>   *matchmuon_phi;
   vector<int>     *matchmuon_charge;
   vector<int>     *matchmuon_type;
   vector<int>     *matchmuon_quality;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_pt;
   vector<float>   *jet_tp_sumpt;
   vector<float>   *jet_trk_sumpt;
   vector<float>   *jet_matchtrk_sumpt;
   vector<float>   *jet_loosematchtrk_sumpt;
   vector<int>     *EMTF_muon_n;
   vector<float>   *EMTF_muon_pt;
   vector<float>   *EMTF_muon_eta;
   vector<float>   *EMTF_muon_phi;
   vector<int>   *EMTF_muon_c;
   vector<int>     *OMTF_muon_n;
   vector<float>   *OMTF_muon_pt;
   vector<float>   *OMTF_muon_eta;
   vector<float>   *OMTF_muon_phi;
   vector<int>   *OMTF_muon_c;
   vector<int>     *BMTF_muon_n;
   vector<float>   *BMTF_muon_pt;
   vector<float>   *BMTF_muon_eta;
   vector<float>   *BMTF_muon_phi;
   vector<int>   *BMTF_muon_c;

   // List of branches
   TBranch        *b_trk_pt;   //!
   TBranch        *b_trk_eta;   //!
   TBranch        *b_trk_phi;   //!
   TBranch        *b_trk_charge;   //!
   TBranch        *b_trk_d0;   //!
   TBranch        *b_trk_z0;   //!
   TBranch        *b_trk_chi2;   //!
   TBranch        *b_trk_nstub;   //!
   TBranch        *b_trk_seed;   //!
   TBranch        *b_trk_genuine;   //!
   TBranch        *b_trk_loose;   //!
   TBranch        *b_trk_unknown;   //!
   TBranch        *b_trk_combinatoric;   //!
   TBranch        *b_trk_fake;   //!
   TBranch        *b_trk_matchtp_pdgid;   //!
   TBranch        *b_trk_matchtp_pt;   //!
   TBranch        *b_trk_matchtp_eta;   //!
   TBranch        *b_trk_matchtp_phi;   //!
   TBranch        *b_trk_matchtp_z0;   //!
   TBranch        *b_trk_matchtp_dxy;   //!
   TBranch        *b_trk_injet;   //!
   TBranch        *b_trk_injet_highpt;   //!
   TBranch        *b_trk_injet_vhighpt;   //!
   TBranch        *b_tp_pt;   //!
   TBranch        *b_tp_eta;   //!
   TBranch        *b_tp_phi;   //!
   TBranch        *b_tp_dxy;   //!
   TBranch        *b_tp_d0;   //!
   TBranch        *b_tp_z0;   //!
   TBranch        *b_tp_d0_prod;   //!
   TBranch        *b_tp_z0_prod;   //!
   TBranch        *b_tp_pdgid;   //!
   TBranch        *b_tp_nmatch;   //!
   TBranch        *b_tp_nloosematch;   //!
   TBranch        *b_tp_nstub;   //!
   TBranch        *b_tp_eventid;   //!
   TBranch        *b_tp_charge;   //!
   TBranch        *b_tp_injet;   //!
   TBranch        *b_tp_injet_highpt;   //!
   TBranch        *b_tp_injet_vhighpt;   //!
   TBranch        *b_matchtrk_pt;   //!
   TBranch        *b_matchtrk_eta;   //!
   TBranch        *b_matchtrk_phi;   //!
   TBranch        *b_matchtrk_z0;   //!
   TBranch        *b_matchtrk_d0;   //!
   TBranch        *b_matchtrk_chi2;   //!
   TBranch        *b_matchtrk_nstub;   //!
   TBranch        *b_matchtrk_seed;   //!
   TBranch        *b_matchtrk_injet;   //!
   TBranch        *b_matchtrk_injet_highpt;   //!
   TBranch        *b_matchtrk_injet_vhighpt;   //!
   TBranch        *b_matchtrk_charge;
   TBranch        *b_loosematchtrk_pt;   //!
   TBranch        *b_loosematchtrk_eta;   //!
   TBranch        *b_loosematchtrk_phi;   //!
   TBranch        *b_loosematchtrk_z0;   //!
   TBranch        *b_loosematchtrk_d0;   //!
   TBranch        *b_loosematchtrk_chi2;   //!
   TBranch        *b_loosematchtrk_nstub;   //!
   TBranch        *b_loosematchtrk_seed;   //!
   TBranch        *b_loosematchtrk_injet;   //!
   TBranch        *b_loosematchtrk_injet_highpt;   //!
   TBranch        *b_loosematchtrk_injet_vhighpt;   //!
   TBranch	  *b_matchmuon_pt; //!
   TBranch	  *b_matchmuon_eta; //!
   TBranch	  *b_matchmuon_phi; //!
   TBranch	  *b_matchmuon_charge; //!
   TBranch	  *b_matchmuon_type; //!
   TBranch	  *b_matchmuon_quality; //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_tp_sumpt;   //!
   TBranch        *b_jet_trk_sumpt;   //!
   TBranch        *b_jet_matchtrk_sumpt;   //!
   TBranch        *b_jet_loosematchtrk_sumpt;   //!
   TBranch        *b_EMTF_muon_n;   //!
   TBranch        *b_EMTF_muon_pt;   //!
   TBranch        *b_EMTF_muon_eta;   //!
   TBranch        *b_EMTF_muon_phi;   //!
   TBranch        *b_EMTF_muon_c;   //!
   TBranch        *b_OMTF_muon_n;   //!
   TBranch        *b_OMTF_muon_pt;   //!
   TBranch        *b_OMTF_muon_eta;   //!
   TBranch        *b_OMTF_muon_phi;   //!
   TBranch        *b_OMTF_muon_c;   //!
   TBranch        *b_BMTF_muon_n;   //!
   TBranch        *b_BMTF_muon_pt;   //!
   TBranch        *b_BMTF_muon_eta;   //!
   TBranch        *b_BMTF_muon_phi;   //!
   TBranch        *b_BMTF_muon_c;   //!

   MuonTrackAnalyzer(bool Propagation_=false, bool GenPropagation_=false, TTree *tree=0);
   virtual ~MuonTrackAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   bool		    Propagation;
   bool		    GenPropagation;
};

#endif

#ifdef MuonTrackAnalyzer_cxx
MuonTrackAnalyzer::MuonTrackAnalyzer(bool Propagation_, bool GenPropagation_,TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("FinalOutput.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("FinalOutput.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("FinalOutput.root:/L1TrackNtuple");
      dir->GetObject("eventTree",tree);

   }
   Propagation=Propagation_;
   GenPropagation=GenPropagation_;
   Init(tree);
}

MuonTrackAnalyzer::~MuonTrackAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MuonTrackAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MuonTrackAnalyzer::LoadTree(Long64_t entry)
{
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

void MuonTrackAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trk_pt = 0;
   trk_eta = 0;
   trk_phi = 0;
   trk_charge = 0;
   trk_d0 = 0;
   trk_z0 = 0;
   trk_chi2 = 0;
   trk_nstub = 0;
   trk_seed = 0;
   trk_genuine = 0;
   trk_loose = 0;
   trk_unknown = 0;
   trk_combinatoric = 0;
   trk_fake = 0;
   trk_matchtp_pdgid = 0;
   trk_matchtp_pt = 0;
   trk_matchtp_eta = 0;
   trk_matchtp_phi = 0;
   trk_matchtp_z0 = 0;
   trk_matchtp_dxy = 0;
   trk_injet = 0;
   trk_injet_highpt = 0;
   trk_injet_vhighpt = 0;
   tp_pt = 0;
   tp_eta = 0;
   tp_phi = 0;
   tp_dxy = 0;
   tp_d0 = 0;
   tp_z0 = 0;
   tp_d0_prod = 0;
   tp_z0_prod = 0;
   tp_pdgid = 0;
   tp_nmatch = 0;
   tp_nloosematch = 0;
   tp_nstub = 0;
   tp_eventid = 0;
   tp_charge = 0;
   tp_injet = 0;
   tp_injet_highpt = 0;
   tp_injet_vhighpt = 0;
   matchtrk_pt = 0;
   matchtrk_eta = 0;
   matchtrk_phi = 0;
   matchtrk_z0 = 0;
   matchtrk_d0 = 0;
   matchtrk_chi2 = 0;
   matchtrk_nstub = 0;
   matchtrk_seed = 0;
   matchtrk_injet = 0;
   matchtrk_injet_highpt = 0;
   matchtrk_injet_vhighpt = 0;
   matchtrk_charge = 0;
   loosematchtrk_pt = 0;
   loosematchtrk_eta = 0;
   loosematchtrk_phi = 0;
   loosematchtrk_z0 = 0;
   loosematchtrk_d0 = 0;
   loosematchtrk_chi2 = 0;
   loosematchtrk_nstub = 0;
   loosematchtrk_seed = 0;
   loosematchtrk_injet = 0;
   loosematchtrk_injet_highpt = 0;
   loosematchtrk_injet_vhighpt = 0;
   matchmuon_pt = 0;
   matchmuon_eta = 0;
   matchmuon_phi = 0;
   matchmuon_charge = 0;
   matchmuon_type = 0;
   matchmuon_quality = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_pt = 0;
   jet_tp_sumpt = 0;
   jet_trk_sumpt = 0;
   jet_matchtrk_sumpt = 0;
   jet_loosematchtrk_sumpt = 0;
   EMTF_muon_n = 0;
   EMTF_muon_pt = 0;
   EMTF_muon_eta = 0;
   EMTF_muon_phi = 0;
   EMTF_muon_c = 0;
   OMTF_muon_n = 0;
   OMTF_muon_pt = 0;
   OMTF_muon_eta = 0;
   OMTF_muon_phi = 0;
   OMTF_muon_c = 0;
   BMTF_muon_n = 0;
   BMTF_muon_pt = 0;
   BMTF_muon_eta = 0;
   BMTF_muon_phi = 0;
   BMTF_muon_c = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
   fChain->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
   fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
   fChain->SetBranchAddress("trk_charge", &trk_charge, &b_trk_charge);
   fChain->SetBranchAddress("trk_d0", &trk_d0, &b_trk_d0);
   fChain->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
   fChain->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
   fChain->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);
   fChain->SetBranchAddress("trk_seed", &trk_seed, &b_trk_seed);
   fChain->SetBranchAddress("trk_genuine", &trk_genuine, &b_trk_genuine);
   fChain->SetBranchAddress("trk_loose", &trk_loose, &b_trk_loose);
   fChain->SetBranchAddress("trk_unknown", &trk_unknown, &b_trk_unknown);
   fChain->SetBranchAddress("trk_combinatoric", &trk_combinatoric, &b_trk_combinatoric);
   fChain->SetBranchAddress("trk_fake", &trk_fake, &b_trk_fake);
   fChain->SetBranchAddress("trk_matchtp_pdgid", &trk_matchtp_pdgid, &b_trk_matchtp_pdgid);
   fChain->SetBranchAddress("trk_matchtp_pt", &trk_matchtp_pt, &b_trk_matchtp_pt);
   fChain->SetBranchAddress("trk_matchtp_eta", &trk_matchtp_eta, &b_trk_matchtp_eta);
   fChain->SetBranchAddress("trk_matchtp_phi", &trk_matchtp_phi, &b_trk_matchtp_phi);
   fChain->SetBranchAddress("trk_matchtp_z0", &trk_matchtp_z0, &b_trk_matchtp_z0);
   fChain->SetBranchAddress("trk_matchtp_dxy", &trk_matchtp_dxy, &b_trk_matchtp_dxy);
   fChain->SetBranchAddress("trk_injet", &trk_injet, &b_trk_injet);
   fChain->SetBranchAddress("trk_injet_highpt", &trk_injet_highpt, &b_trk_injet_highpt);
   fChain->SetBranchAddress("trk_injet_vhighpt", &trk_injet_vhighpt, &b_trk_injet_vhighpt);
   fChain->SetBranchAddress("tp_pt", &tp_pt, &b_tp_pt);
   fChain->SetBranchAddress("tp_eta", &tp_eta, &b_tp_eta);
   fChain->SetBranchAddress("tp_phi", &tp_phi, &b_tp_phi);
   fChain->SetBranchAddress("tp_dxy", &tp_dxy, &b_tp_dxy);
   fChain->SetBranchAddress("tp_d0", &tp_d0, &b_tp_d0);
   fChain->SetBranchAddress("tp_z0", &tp_z0, &b_tp_z0);
   fChain->SetBranchAddress("tp_d0_prod", &tp_d0_prod, &b_tp_d0_prod);
   fChain->SetBranchAddress("tp_z0_prod", &tp_z0_prod, &b_tp_z0_prod);
   fChain->SetBranchAddress("tp_pdgid", &tp_pdgid, &b_tp_pdgid);
   fChain->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
   fChain->SetBranchAddress("tp_nloosematch", &tp_nloosematch, &b_tp_nloosematch);
   fChain->SetBranchAddress("tp_nstub", &tp_nstub, &b_tp_nstub);
   fChain->SetBranchAddress("tp_eventid", &tp_eventid, &b_tp_eventid);
   fChain->SetBranchAddress("tp_charge", &tp_charge, &b_tp_charge);
   fChain->SetBranchAddress("tp_injet", &tp_injet, &b_tp_injet);
   fChain->SetBranchAddress("tp_injet_highpt", &tp_injet_highpt, &b_tp_injet_highpt);
   fChain->SetBranchAddress("tp_injet_vhighpt", &tp_injet_vhighpt, &b_tp_injet_vhighpt);
   fChain->SetBranchAddress("matchtrk_pt", &matchtrk_pt, &b_matchtrk_pt);
   fChain->SetBranchAddress("matchtrk_eta", &matchtrk_eta, &b_matchtrk_eta);
   fChain->SetBranchAddress("matchtrk_phi", &matchtrk_phi, &b_matchtrk_phi);
   fChain->SetBranchAddress("matchtrk_z0", &matchtrk_z0, &b_matchtrk_z0);
   fChain->SetBranchAddress("matchtrk_d0", &matchtrk_d0, &b_matchtrk_d0);
   fChain->SetBranchAddress("matchtrk_chi2", &matchtrk_chi2, &b_matchtrk_chi2);
   fChain->SetBranchAddress("matchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);
   fChain->SetBranchAddress("matchtrk_seed", &matchtrk_seed, &b_matchtrk_seed);
   fChain->SetBranchAddress("matchtrk_injet", &matchtrk_injet, &b_matchtrk_injet);
   fChain->SetBranchAddress("matchtrk_injet_highpt", &matchtrk_injet_highpt, &b_matchtrk_injet_highpt);
   fChain->SetBranchAddress("matchtrk_injet_vhighpt", &matchtrk_injet_vhighpt, &b_matchtrk_injet_vhighpt);
   fChain->SetBranchAddress("matchtrk_charge", &matchtrk_charge, &b_matchtrk_charge);
   fChain->SetBranchAddress("loosematchtrk_pt", &loosematchtrk_pt, &b_loosematchtrk_pt);
   fChain->SetBranchAddress("loosematchtrk_eta", &loosematchtrk_eta, &b_loosematchtrk_eta);
   fChain->SetBranchAddress("loosematchtrk_phi", &loosematchtrk_phi, &b_loosematchtrk_phi);
   fChain->SetBranchAddress("loosematchtrk_z0", &loosematchtrk_z0, &b_loosematchtrk_z0);
   fChain->SetBranchAddress("loosematchtrk_d0", &loosematchtrk_d0, &b_loosematchtrk_d0);
   fChain->SetBranchAddress("loosematchtrk_chi2", &loosematchtrk_chi2, &b_loosematchtrk_chi2);
   fChain->SetBranchAddress("loosematchtrk_nstub", &loosematchtrk_nstub, &b_loosematchtrk_nstub);
   fChain->SetBranchAddress("loosematchtrk_seed", &loosematchtrk_seed, &b_loosematchtrk_seed);
   fChain->SetBranchAddress("loosematchtrk_injet", &loosematchtrk_injet, &b_loosematchtrk_injet);
   fChain->SetBranchAddress("loosematchtrk_injet_highpt", &loosematchtrk_injet_highpt, &b_loosematchtrk_injet_highpt);
   fChain->SetBranchAddress("loosematchtrk_injet_vhighpt", &loosematchtrk_injet_vhighpt, &b_loosematchtrk_injet_vhighpt);
   fChain->SetBranchAddress("matchmuon_pt", &matchmuon_pt, &b_matchmuon_pt);
   fChain->SetBranchAddress("matchmuon_eta", &matchmuon_eta, &b_matchmuon_eta);
   fChain->SetBranchAddress("matchmuon_phi", &matchmuon_phi, &b_matchmuon_phi);
   fChain->SetBranchAddress("matchmuon_charge", &matchmuon_charge, &b_matchmuon_charge);
   fChain->SetBranchAddress("matchmuon_type", &matchmuon_type, &b_matchmuon_type);
   fChain->SetBranchAddress("matchmuon_quality", &matchmuon_quality, &b_matchmuon_quality);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_tp_sumpt", &jet_tp_sumpt, &b_jet_tp_sumpt);
   fChain->SetBranchAddress("jet_trk_sumpt", &jet_trk_sumpt, &b_jet_trk_sumpt);
   fChain->SetBranchAddress("jet_matchtrk_sumpt", &jet_matchtrk_sumpt, &b_jet_matchtrk_sumpt);
   fChain->SetBranchAddress("jet_loosematchtrk_sumpt", &jet_loosematchtrk_sumpt, &b_jet_loosematchtrk_sumpt);
   fChain->SetBranchAddress("EMTF_muon_n", &EMTF_muon_n, &b_EMTF_muon_n);
   fChain->SetBranchAddress("EMTF_muon_pt", &EMTF_muon_pt, &b_EMTF_muon_pt);
   fChain->SetBranchAddress("EMTF_muon_eta", &EMTF_muon_eta, &b_EMTF_muon_eta);
   fChain->SetBranchAddress("EMTF_muon_phi", &EMTF_muon_phi, &b_EMTF_muon_phi);
   fChain->SetBranchAddress("EMTF_muon_c", &EMTF_muon_c, &b_EMTF_muon_c);
   fChain->SetBranchAddress("OMTF_muon_n", &OMTF_muon_n, &b_OMTF_muon_n);
   fChain->SetBranchAddress("OMTF_muon_pt", &OMTF_muon_pt, &b_OMTF_muon_pt);
   fChain->SetBranchAddress("OMTF_muon_eta", &OMTF_muon_eta, &b_OMTF_muon_eta);
   fChain->SetBranchAddress("OMTF_muon_phi", &OMTF_muon_phi, &b_OMTF_muon_phi);
   fChain->SetBranchAddress("OMTF_muon_c", &OMTF_muon_c, &b_OMTF_muon_c);
   fChain->SetBranchAddress("BMTF_muon_n", &BMTF_muon_n, &b_BMTF_muon_n);
   fChain->SetBranchAddress("BMTF_muon_pt", &BMTF_muon_pt, &b_BMTF_muon_pt);
   fChain->SetBranchAddress("BMTF_muon_eta", &BMTF_muon_eta, &b_BMTF_muon_eta);
   fChain->SetBranchAddress("BMTF_muon_phi", &BMTF_muon_phi, &b_BMTF_muon_phi);
   fChain->SetBranchAddress("BMTF_muon_c", &BMTF_muon_c, &b_BMTF_muon_c);
   Notify();
}

Bool_t MuonTrackAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MuonTrackAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MuonTrackAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MuonTrackAnalyzer_cxx
