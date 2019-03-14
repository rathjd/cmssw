#define MuonTrackAnalyzer_cxx
#include "MuonTrackAnalyzer.h"
#include <TH2.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <iostream>
#include "TH3F.h"
#include "TVector2.h"
#include <math.h>

struct tpMatch{
  unsigned tpAddress;
  TLorentzVector tpVec;
  std::vector<std::pair<unsigned,float> > match; //muon address, DeltaR to tp
};

bool sortbysec(const pair<unsigned,float> &a, 
              const pair<unsigned,float> &b) 
{ 
    if(a.second!=b.second) return (a.second < b.second); 
    else return (a.first<b.first);
} 

void MuonTrackAnalyzer::Loop()
{



//   In a ROOT session, you can do:
//      root> .L MuonTrackAnalyzer.C
//      root> MuonTrackAnalyzer t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   //nentries=10;
   bool verbose=false;

   //histogram declarations   
   TProfile2D *minDeltaR_tpToMuon = new TProfile2D("minDeltaR_tpToMuon","matching tp to reco muon;|#eta|(tp);p_{T}(tp);min(#Delta R)", 25, 0.0, 2.5, 50, 0., 100.);
   TH3F *minDeltaR3D_tpToMuon = new TH3F("minDeltaR3D_tpToMuon","matching tp to reco muon;|#eta|(tp);p_{T}(tp);min(#Delta R)", 25, 0.0, 2.5, 50, 0., 100., 100, 0., 1.);

   TH3F *DeltaPhiWindow = new TH3F("DeltaPhiWindow","matched to same tp #Delta#phi  ;|#eta(tp)|;p_{T}(tp);|#Delta#phi|"  ,25,0.0,2.5,20,0.,100.,100,0.,1.);
   TH3F *DeltaEtaWindow = new TH3F("DeltaEtaWindow","matched to same tp #Delta#eta  ;|#eta(tp)|;p_{T}(tp);|#Delta#eta|"  ,25,0.0,2.5,20,0.,100.,100,0.,1.);
   TH3F *DeltaPtWindow =  new TH3F("DeltaPtWindow" ,"matched to same tp #Delta p_{T};|#eta(tp)|;p_{T}(tp);|#Delta p_{T}|",25,0.0,2.5,20,0.,100.,75,0.,75.);
   TProfile2D *ChargeVeracity=  new TProfile2D("ChargeVeracity","matched to same tp charge veracity;#eta(tp);p_{T}(tp);sign(c_{1}#times c_{2})",25,0.0,2.5,20,0.,100.);

   TH2F *TrkResolutionB=  new TH2F("TrkResolutionB","tp to matchtrk resolution in barrel;p_{T}(tp);p_{T}(trk)",60,0.,120.,60,0.,120.);
   TH2F *TrkResolutionO=  new TH2F("TrkResolutionO","tp to matchtrk resolution in overlap;p_{T}(tp);p_{T}(trk)",60,0.,120.,60,0.,120.);
   TH2F *TrkResolutionE=  new TH2F("TrkResolutionE","tp to matchtrk resolution in endcaps;p_{T}(tp);p_{T}(trk)",60,0.,120.,60,0.,120.);
   TH2F *MuResolutionB=  new TH2F("MuResolutionB","tp to muon resolution in barrel;p_{T}(tp);p_{T}(trk)",60,0.,120.,60,0.,120.);
   TH2F *MuResolutionO=  new TH2F("MuResolutionO","tp to muon resolution in overlap;p_{T}(tp);p_{T}(trk)",60,0.,120.,60,0.,120.);
   TH2F *MuResolutionE=  new TH2F("MuResolutionE","tp to muon resolution in endcaps;p_{T}(tp);p_{T}(trk)",60,0.,120.,60,0.,120.);

   TH2F *TrkResolutionB_phi=  new TH2F("TrkResolutionB_phi","tp to matchtrk #phi resolution in barrel;#phi(tp);#phi(trk)",60,-M_PI,M_PI,60,-M_PI,M_PI);
   TH2F *TrkResolutionO_phi=  new TH2F("TrkResolutionO_phi","tp to matchtrk #phi resolution in overlap;#phi(tp);#phi(trk)",60,-M_PI,M_PI,60,-M_PI,M_PI);
   TH2F *TrkResolutionE_phi=  new TH2F("TrkResolutionE_phi","tp to matchtrk #phi resolution in endcaps;#phi(tp);#phi(trk)",60,-M_PI,M_PI,60,-M_PI,M_PI);
   TH2F *MuResolutionB_phi=  new TH2F("MuResolutionB_phi","tp to muon #phi resolution in barrel;#phi(tp);#phi(trk)",60,-M_PI,M_PI,60,-M_PI,M_PI);
   TH2F *MuResolutionO_phi=  new TH2F("MuResolutionO_phi","tp to muon #phi resolution in overlap;#phi(tp);#phi(trk)",60,-M_PI,M_PI,60,-M_PI,M_PI);
   TH2F *MuResolutionE_phi=  new TH2F("MuResolutionE_phi","tp to muon #phi resolution in endcaps;#phi(tp);#phi(trk)",60,-M_PI,M_PI,60,-M_PI,M_PI);


   TH2F *TrkResolutionB_eta=  new TH2F("TrkResolutionB_eta","tp to matchtrk #eta resolution in barrel;#eta(tp);#eta(trk)",50,-2.5,2.5,50,-2.5,2.5);
   TH2F *TrkResolutionO_eta=  new TH2F("TrkResolutionO_eta","tp to matchtrk #eta resolution in overlap;#eta(tp);#eta(trk)",50,-2.5,2.5,50,-2.5,2.5);
   TH2F *TrkResolutionE_eta=  new TH2F("TrkResolutionE_eta","tp to matchtrk #eta resolution in endcaps;#eta(tp);#eta(trk)",50,-2.5,2.5,50,-2.5,2.5);
   TH2F *MuResolutionB_eta=  new TH2F("MuResolutionB_eta","tp to muon #eta resolution in barrel;#eta(tp);#eta(trk)",50,-2.5,2.5,50,-2.5,2.5);
   TH2F *MuResolutionO_eta=  new TH2F("MuResolutionO_eta","tp to muon #eta resolution in overlap;#eta(tp);#eta(trk)",50,-2.5,2.5,50,-2.5,2.5);
   TH2F *MuResolutionE_eta=  new TH2F("MuResolutionE_eta","tp to muon #eta resolution in endcaps;#eta(tp);#eta(trk)",50,-2.5,2.5,50,-2.5,2.5);

   //actual event loop
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);

      if(jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      std::vector<tpMatch> Matcher;
      if(verbose) std::cout<<"----------event "<<jentry<<"----------"<<std::endl;
      for(unsigned tp=0; tp<tp_pdgid->size(); ++tp){
      	if(abs(tp_pdgid->at(tp))!=13) continue; //skip nonmuons
	TLorentzVector Current;
	if(verbose) std::cout<<"default for "<<tp<<": ("<<tp_pt->at(tp)<<", "<<tp_eta->at(tp)<<", "<<tp_phi->at(tp)<<")"<<std::endl;
	if(GenPropagation){
	  if(tp_pt->at(tp)*cosh(fabs(tp_eta->at(tp)))<3.5) continue;
	  if(fabs(tp_eta->at(tp))<1.1 && tp_pt->at(tp)<3.5) continue;
	  if(fabs(tp_eta->at(tp))>2.5) continue;
	  if(verbose)std::cout<<"accepted for propagation"<<std::endl;
	  float prop_eta=0.;
	  float prop_pt=tp_pt->at(tp);
	  float prop_phi=0.;
	  float dzCorrPhi=1.;
	  float etaProp=fabs(tp_eta->at(tp));
	  if(fabs(tp_eta->at(tp))<1.1){
	    etaProp=1.1;
	    prop_eta=tp_eta->at(tp)+tp_z0->at(tp)/550./cosh(fabs(tp_eta->at(tp)));
	  }
	  else{
	    dzCorrPhi=1.+pow(-1.,1+signbit(tp_eta->at(tp)))*tp_z0->at(tp)/850.;
	    prop_eta=tp_eta->at(tp)+tp_z0->at(tp)/850./(1.+pow(-1.,1+signbit(tp_eta->at(tp)))*tp_z0->at(tp)/850.)*tanh(tp_eta->at(tp));
	  }
	  prop_phi=TVector2::Phi_mpi_pi(tp_phi->at(tp)-1.464*tp_charge->at(tp)*cosh(1.7)/cosh(etaProp)/tp_pt->at(tp)*dzCorrPhi-M_PI/144.);
	  Current.SetPtEtaPhiM(prop_pt,prop_eta,prop_phi,0.106);
	  if(verbose)std::cout<<"propagated for "<<tp<<": ("<<prop_pt<<", "<<prop_eta<<", "<<prop_phi<<")"<<std::endl;
	}
	else Current.SetPtEtaPhiM(tp_pt->at(tp),tp_eta->at(tp),tp_phi->at(tp),0.106);
	std::vector<std::pair<unsigned,float> > dummy;
	tpMatch CurrentTp={tp,Current,dummy};
	Matcher.push_back(CurrentTp);
      }

      //code type of muon into 0xx E, 1xx O, 2xx B
      for(unsigned emu=0; emu<EMTF_muon_pt->size(); ++emu){
	TLorentzVector CurrentMu;
	CurrentMu.SetPtEtaPhiM(EMTF_muon_pt->at(emu),EMTF_muon_eta->at(emu),EMTF_muon_phi->at(emu),0.106);
	for(unsigned tp=0; tp<Matcher.size(); ++tp) Matcher[tp].match.push_back(std::make_pair(emu,Matcher[tp].tpVec.DeltaR(CurrentMu)));
      }
      for(unsigned omu=0; omu<OMTF_muon_pt->size(); ++omu){
	TLorentzVector CurrentMu;
	CurrentMu.SetPtEtaPhiM(OMTF_muon_pt->at(omu),OMTF_muon_eta->at(omu),OMTF_muon_phi->at(omu),0.106);
	for(unsigned tp=0; tp<Matcher.size(); ++tp) Matcher[tp].match.push_back(std::make_pair(100+omu,Matcher[tp].tpVec.DeltaR(CurrentMu)));
      }
      for(unsigned bmu=0; bmu<BMTF_muon_pt->size(); ++bmu){
	TLorentzVector CurrentMu;
	CurrentMu.SetPtEtaPhiM(BMTF_muon_pt->at(bmu),BMTF_muon_eta->at(bmu),BMTF_muon_phi->at(bmu),0.106);
	for(unsigned tp=0; tp<Matcher.size(); ++tp) Matcher[tp].match.push_back(std::make_pair(200+bmu,Matcher[tp].tpVec.DeltaR(CurrentMu)));
      }

      //loop over all tps and sort by lowest DeltaR value
      for(unsigned tp=0; tp<Matcher.size(); ++tp) std::sort(Matcher[tp].match.begin(), Matcher[tp].match.end(),sortbysec);

      //for(unsigned tp=0; tp<Matcher.size(); ++tp) for(unsigned m=0; m<Matcher[tp].match.size(); ++m) std::cout<<"tp "<<tp<<", m "<<m<<": ("<<Matcher[tp].match[m].first<<","<<Matcher[tp].match[m].second<<")"<<std::endl; //sorttest cout
      
      //find best DeltaR match, remove from all others
      std::vector<std::pair<unsigned,unsigned> > FinalMatch; //tp address, muAddress +200 for barrel, +100 for overlap, +0 for endcap
      while(Matcher.size()){
      	float minDeltaR=999;
	std::pair<unsigned,unsigned> currentMatchPair;
	int position=-1;
      	for(unsigned tp=0; tp<Matcher.size(); ++tp){
	  if(Matcher[tp].match.size()) if(Matcher[tp].match[0].second<minDeltaR){
	    minDeltaR=Matcher[tp].match[0].second;
	    currentMatchPair=std::make_pair(Matcher[tp].tpAddress,Matcher[tp].match[0].first);
	    position=tp;
	  }
	}
	minDeltaR_tpToMuon->Fill(tp_eta->at(currentMatchPair.first),tp_pt->at(currentMatchPair.first),minDeltaR);
	minDeltaR3D_tpToMuon->Fill(tp_eta->at(currentMatchPair.first),tp_pt->at(currentMatchPair.first),minDeltaR);
	if(minDeltaR<0.2) FinalMatch.push_back(currentMatchPair);
	if(position==-1) break;
	Matcher.erase(Matcher.begin()+position); //removes tracking particle
	//remove matched muon from all lists
	for(unsigned tp=0; tp<Matcher.size(); ++tp) for(unsigned m=0; m<Matcher[tp].match.size(); ++m) if(Matcher[tp].match[m].first==currentMatchPair.first){
	  Matcher[tp].match.erase(Matcher[tp].match.begin()+m);
	  --m;
	}
      }

      //fill resolution plots
      for(unsigned f=0; f<FinalMatch.size(); ++f){
	const unsigned tp=FinalMatch[f].first;
	unsigned mu=FinalMatch[f].second;
	const float tpEta=fabs(tp_eta->at(tp));
	const float tpEtaS=tp_eta->at(tp);
	const float tppT=tp_pt->at(tp);
	const float tpPhi=tp_phi->at(tp);
        if(matchtrk_pt->at(tp)<0.) continue;

	//propagated variables for matchtrks;
	float prop_eta=0.;
	float prop_pt=matchtrk_pt->at(tp);
	float prop_phi=0.;
	if(Propagation){
	  if(matchtrk_pt->at(tp)*cosh(fabs(matchtrk_eta->at(tp)))<3.5) continue;
	  if(fabs(matchtrk_eta->at(tp))<1.1 && matchtrk_pt->at(tp)<3.5) continue;
	  if(fabs(matchtrk_eta->at(tp))>2.5) continue;
	  float dzCorrPhi=1.;
	  float etaProp=fabs(matchtrk_eta->at(tp));
	  if(fabs(matchtrk_eta->at(tp))<1.1){
	    etaProp=1.1;
	    prop_eta=matchtrk_eta->at(tp)+matchtrk_z0->at(tp)/550./cosh(fabs(matchtrk_eta->at(tp)));
	  }
	  else{
	    dzCorrPhi=1.+pow(-1.,1+signbit(matchtrk_eta->at(tp)))*matchtrk_z0->at(tp)/850.;
	    prop_eta=matchtrk_eta->at(tp)+matchtrk_z0->at(tp)/850./(1.+pow(-1.,1+signbit(matchtrk_eta->at(tp)))*matchtrk_z0->at(tp)/850.)*tanh(matchtrk_eta->at(tp));
	  }
	  prop_phi=TVector2::Phi_mpi_pi(matchtrk_phi->at(tp)-1.464*pow(-1.,matchtrk_charge->at(tp))*cosh(1.7)/cosh(etaProp)/matchtrk_pt->at(tp)*dzCorrPhi-M_PI/144.);
	}


	//std::cout<<matchtrk_pt->at(tp)<<" ?= ";
	float mtrk_phi=matchtrk_phi->at(tp);
	float mtrk_eta=matchtrk_eta->at(tp);
	float mtrk_pt =matchtrk_pt->at(tp);
	if(Propagation){
	  mtrk_phi=prop_phi;
	  mtrk_eta=prop_eta;
	  mtrk_pt=prop_pt;
	}

	if(mu>=200){ //is barrel muon
	  mu-=200;
	  DeltaPhiWindow->Fill(tpEta,tppT,fabs(TVector2::Phi_mpi_pi(mtrk_phi-BMTF_muon_phi->at(mu))));
	  DeltaEtaWindow->Fill(tpEta,tppT,fabs(mtrk_eta-BMTF_muon_eta->at(mu)));
	  DeltaPtWindow-> Fill(tpEta,tppT,fabs(mtrk_pt-BMTF_muon_pt->at(mu)));
	  TrkResolutionB->Fill(tppT,mtrk_pt);
	  MuResolutionB->Fill(tppT,BMTF_muon_pt->at(mu));
	  TrkResolutionB_phi->Fill(tpPhi,mtrk_phi);
	  MuResolutionB_phi->Fill(tpPhi,TVector2::Phi_mpi_pi(BMTF_muon_phi->at(mu)));
	  TrkResolutionB_eta->Fill(tpEtaS,mtrk_eta);
	  MuResolutionB_eta->Fill(tpEtaS,BMTF_muon_eta->at(mu));
	  //std::cout<<matchtrk_charge->at(tp)<<"?="<<BMTF_muon_c->at(mu)<<std::endl;
	  if((matchtrk_charge->at(tp) && BMTF_muon_c->at(mu)==-1) || (!matchtrk_charge->at(tp) && BMTF_muon_c->at(mu)==1)) ChargeVeracity->Fill(tpEta,tppT,1.);
	  else ChargeVeracity->Fill(tpEta,tppT,0.);
	  //std::cout<<BMTF_muon_pt->at(mu)<<std::endl;
	}
	else if(mu>=100){ //is overlap muon
	  mu-=100;
	  DeltaPhiWindow->Fill(tpEta,tppT,fabs(TVector2::Phi_mpi_pi(mtrk_phi-OMTF_muon_phi->at(mu))));
	  DeltaEtaWindow->Fill(tpEta,tppT,fabs(mtrk_eta-OMTF_muon_eta->at(mu)));
	  DeltaPtWindow-> Fill(tpEta,tppT,fabs(mtrk_pt-OMTF_muon_pt->at(mu)));
	  TrkResolutionO->Fill(tppT,mtrk_pt);
	  MuResolutionO->Fill(tppT,OMTF_muon_pt->at(mu));
	  TrkResolutionO_phi->Fill(tpPhi,mtrk_phi);
	  MuResolutionO_phi->Fill(tpPhi,TVector2::Phi_mpi_pi(OMTF_muon_phi->at(mu)));
	  TrkResolutionO_eta->Fill(tpEtaS,mtrk_eta);
	  MuResolutionO_eta->Fill(tpEtaS,OMTF_muon_eta->at(mu));
	  if((matchtrk_charge->at(tp) && OMTF_muon_c->at(mu)==-1) ||(!matchtrk_charge->at(tp) && OMTF_muon_c->at(mu)==1)) ChargeVeracity->Fill(tpEta,tppT,1.);
	  else ChargeVeracity->Fill(tpEta,tppT,0.);
	  //std::cout<<OMTF_muon_pt->at(mu)<<std::endl;
	}
	else{ //is endcap muon
	  DeltaPhiWindow->Fill(tpEta,tppT,fabs(TVector2::Phi_mpi_pi(mtrk_phi-EMTF_muon_phi->at(mu))));
	  DeltaEtaWindow->Fill(tpEta,tppT,fabs(mtrk_eta-EMTF_muon_eta->at(mu)));
	  DeltaPtWindow-> Fill(tpEta,tppT,fabs(mtrk_pt-EMTF_muon_pt->at(mu)));	  
	  TrkResolutionE->Fill(tppT,mtrk_pt);
	  MuResolutionE->Fill(tppT,EMTF_muon_pt->at(mu));
	  TrkResolutionE_phi->Fill(tpPhi,mtrk_phi);
	  MuResolutionE_phi->Fill(tpPhi,TVector2::Phi_mpi_pi(EMTF_muon_phi->at(mu)));
	  TrkResolutionE_eta->Fill(tpEtaS,mtrk_eta);
	  MuResolutionE_eta->Fill(tpEtaS,EMTF_muon_eta->at(mu));
	  if((matchtrk_charge->at(tp) && EMTF_muon_c->at(mu)==-1) ||(!matchtrk_charge->at(tp) && EMTF_muon_c->at(mu)==1)) ChargeVeracity->Fill(tpEta,tppT,1.);
	  else ChargeVeracity->Fill(tpEta,tppT,0.);
	  //std::cout<<EMTF_muon_pt->at(mu)<<std::endl;
	}
      }
   }

   TFile *savefile;
   if(Propagation && !GenPropagation)     savefile=new TFile("Propagated_savedHists.root","RECREATE");
   else if(!Propagation && GenPropagation)savefile=new TFile("GenPropagated_savedHists.root","RECREATE");
   else if(GenPropagation && Propagation) savefile=new TFile("AllPropagated_savedHists.root","RECREATE");
   else		  			  savefile=new TFile("Unpropagated_savedHists.root","RECREATE");
   minDeltaR_tpToMuon->Write();
   minDeltaR3D_tpToMuon->Write();
   DeltaPhiWindow->Write();
   DeltaPtWindow->Write();
   DeltaEtaWindow->Write();

   //derive nth percentile representations of windows
   TH2F *DeltaXbase=(TH2F*)DeltaPhiWindow->Project3D("yx");
   TH2F *DeltaPhi99=(TH2F*)DeltaXbase->Clone("DeltaPhi99");
   DeltaPhi99->SetTitle("#Delta#phi 99th quantile");
   TH2F *DeltaPhi95=(TH2F*)DeltaXbase->Clone("DeltaPhi95");
   DeltaPhi95->SetTitle("#Delta#phi 95th quantile");
   TH2F *DeltaPhi5=(TH2F*)DeltaXbase->Clone("DeltaPhi5");
   DeltaPhi5->SetTitle("#Delta#phi 5th quantile");
   TH2F *DeltaPhi1=(TH2F*)DeltaXbase->Clone("DeltaPhi1");
   DeltaPhi1->SetTitle("#Delta#phi 1st quantile");
   TH2F *DeltaEta99=(TH2F*)DeltaXbase->Clone("DeltaEta99");
   DeltaEta99->SetTitle("#Delta#eta 99th quantile");
   TH2F *DeltaEta95=(TH2F*)DeltaXbase->Clone("DeltaEta95");
   DeltaEta95->SetTitle("#Delta#eta 95th quantile");
   TH2F *DeltaEta5=(TH2F*)DeltaXbase->Clone("DeltaEta5");
   DeltaEta5->SetTitle("#Delta#eta 5th quantile");
   TH2F *DeltaEta1=(TH2F*)DeltaXbase->Clone("DeltaEta1");
   DeltaEta1->SetTitle("#Delta#eta 1st quantile");
   TH2F *DeltaPt99=(TH2F*)DeltaXbase->Clone("DeltaPt99");
   DeltaPt99->SetTitle("#Delta p_{T} 99th quantile");
   TH2F *DeltaPt95=(TH2F*)DeltaXbase->Clone("DeltaPt95");
   DeltaPt95->SetTitle("#Delta p_{T} 95th quantile");
   TH2F *DeltaPt5=(TH2F*)DeltaXbase->Clone("DeltaPt5");
   DeltaPt5->SetTitle("#Delta p_{T} 5th quantile");
   TH2F *DeltaPt1=(TH2F*)DeltaXbase->Clone("DeltaPt1");
   DeltaPt1->SetTitle("#Delta p_{T} 1st quantile");
   for(unsigned x=1; x<DeltaXbase->GetNbinsX()+1; ++x) for(unsigned y=1; y<DeltaXbase->GetNbinsY()+1; ++y){
     Double_t quantileTargets[4]={0.01,0.05,0.95,0.99};
     Double_t quantileResults[4];
     TH1F *Phi=(TH1F*)DeltaPhiWindow->ProjectionZ("Phi",x,x,y,y);
     //std::cout<<"("<<x<<","<<y<<"): ["<<DeltaXbase->GetXaxis()->GetBinLowEdge(x)<<","<<DeltaXbase->GetYaxis()->GetBinLowEdge(y)<<"] ?= ["<<DeltaPhiWindow->GetXaxis()->GetBinLowEdge(x)<<","<<DeltaPhiWindow->GetYaxis()->GetBinLowEdge(y)<<"]"<<std::endl;
     if(!Phi->Integral(0,-1)){
       DeltaPhi1->SetBinContent(x,y,0.);
       DeltaPhi5->SetBinContent(x,y,0.);
       DeltaPhi95->SetBinContent(x,y,0.);
       DeltaPhi99->SetBinContent(x,y,0.);
     }
     else{
       Phi->GetQuantiles(4,quantileResults,quantileTargets);
       DeltaPhi1->SetBinContent(x,y,quantileResults[0]);
       DeltaPhi5->SetBinContent(x,y,quantileResults[1]);
       DeltaPhi95->SetBinContent(x,y,quantileResults[2]);
       DeltaPhi99->SetBinContent(x,y,quantileResults[3]);
     }

     TH1F *Eta=(TH1F*)DeltaEtaWindow->ProjectionZ("Eta",x,x,y,y);
     if(!Eta->Integral(0,-1)){
       DeltaEta1->SetBinContent(x,y,0.);
       DeltaEta5->SetBinContent(x,y,0.);
       DeltaEta95->SetBinContent(x,y,0.);
       DeltaEta99->SetBinContent(x,y,0.);
     }
     else{
       Eta->GetQuantiles(4,quantileResults,quantileTargets);
       DeltaEta1->SetBinContent(x,y,quantileResults[0]);
       DeltaEta5->SetBinContent(x,y,quantileResults[1]);
       DeltaEta95->SetBinContent(x,y,quantileResults[2]);
       DeltaEta99->SetBinContent(x,y,quantileResults[3]);
     }

     TH1F *Pt=(TH1F*)DeltaPtWindow->ProjectionZ("Pt",x,x,y,y);
     if(!Pt->Integral(0,-1)){
       DeltaPt1->SetBinContent(x,y,0.);
       DeltaPt5->SetBinContent(x,y,0.);
       DeltaPt95->SetBinContent(x,y,0.);
       DeltaPt99->SetBinContent(x,y,0.);
     }
     else{
       Pt->GetQuantiles(4,quantileResults,quantileTargets);
       DeltaPt1->SetBinContent(x,y,quantileResults[0]);
       DeltaPt5->SetBinContent(x,y,quantileResults[1]);
       DeltaPt95->SetBinContent(x,y,quantileResults[2]);
       DeltaPt99->SetBinContent(x,y,quantileResults[3]);
     }
   }
   DeltaPhi99->Write();
   DeltaPhi95->Write();
   DeltaPhi5->Write();
   DeltaPhi1->Write();
   DeltaEta99->Write();
   DeltaEta95->Write();
   DeltaEta5->Write();
   DeltaEta1->Write();
   DeltaPt99->Write();
   DeltaPt95->Write();
   DeltaPt5->Write();
   DeltaPt1->Write();
   ChargeVeracity->Write();
   MuResolutionB->Write();
   MuResolutionO->Write();
   MuResolutionE->Write();
   TrkResolutionB->Write();
   TrkResolutionO->Write();
   TrkResolutionE->Write();
   MuResolutionB_phi->Write();
   MuResolutionO_phi->Write();
   MuResolutionE_phi->Write();
   TrkResolutionB_phi->Write();
   TrkResolutionO_phi->Write();
   TrkResolutionE_phi->Write();
   MuResolutionB_eta->Write();
   MuResolutionO_eta->Write();
   MuResolutionE_eta->Write();
   TrkResolutionB_eta->Write();
   TrkResolutionO_eta->Write();
   TrkResolutionE_eta->Write();
   
   savefile->Close();
}
