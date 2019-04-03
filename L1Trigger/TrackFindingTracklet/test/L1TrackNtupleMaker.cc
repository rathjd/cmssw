//////////////////////////////////////////////////////////////////////
//                                                                  //
//  Analyzer for making mini-ntuple for L1 track performance plots  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

///////////////////////
// DATA FORMATS HEADERS
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

////////////////
// PHYSICS TOOLS
#include "CommonTools/UtilAlgos/interface/TFileService.h"

////////////////
// PHYSICS TOOLS
#include "GEMCode/GEMValidation/interface/SimTrackMatchManager.h"

///////////////
// ROOT HEADERS
#include <TROOT.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>

//////////////
// STD HEADERS
#include <memory>
#include <string>
#include <iostream>
#include <math.h>

//////////////
// additional Sven STD headers
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <set>
#include <sstream>
#include <bitset>

//////////////
// additional Sven muon headers
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "L1Trigger/L1TMuon/interface/MuonRawDigiTranslator.h"
#include "L1Trigger/L1TMuon/interface/RegionalMuonRawDigiTranslator.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"

//////////////
// NAMESPACES
using namespace std;
using namespace edm;


//////////////////////////////
//                          //
//     CLASS DEFINITION     //
//                          //
//////////////////////////////

class L1TrackNtupleMaker : public edm::EDAnalyzer
{
public:

  // Constructor/destructor
  explicit L1TrackNtupleMaker(const edm::ParameterSet& iConfig);
  virtual ~L1TrackNtupleMaker();

  // Mandatory methods
  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

protected:

private:

  //-----------------------------------------------------------------------------------------------
  // Containers of parameters passed by python configuration file
  edm::ParameterSet config;

  int MyProcess;        // 11/13/211 for single electrons/muons/pions, 6/15 for pions from ttbar/taus, 1 for inclusive
  bool DebugMode;       // lots of debug printout statements
  bool SaveAllTracks;   // store in ntuples not only truth-matched tracks but ALL tracks
  bool SaveStubs;       // option to save also stubs in the ntuples (makes them large...)
  int L1Tk_nPar;        // use 4 or 5 parameter track fit?
  int TP_minNStub;      // require TPs to have >= minNStub (defining efficiency denominator) (==0 means to only require >= 1 cluster)
  int TP_minNStubLayer; // require TPs to have stubs in >= minNStubLayer layers/disks (defining efficiency denominator)
  double TP_minPt;      // save TPs with pt > minPt
  double TP_maxEta;     // save TPs with |eta| < maxEta
  double TP_maxZ0;      // save TPs with |z0| < maxZ0
  int L1Tk_minNStub;    // require L1 tracks to have >= minNStub (this is mostly for tracklet purposes)

  bool TrackingInJets;  // do tracking in jets?

  bool SaveTracklet;    // save some tracklet-dedicated variables, turned off by default

  //additional muon stuff
  bool SaveMuons;	//save muons in addition


  edm::InputTag L1TrackInputTag;        // L1 track collection
  edm::InputTag MCTruthTrackInputTag;   // MC truth collection
  edm::InputTag MCTruthClusterInputTag;
  edm::InputTag L1StubInputTag;
  edm::InputTag MCTruthStubInputTag;
  edm::InputTag TrackingParticleInputTag;
  edm::InputTag TrackingVertexInputTag;
  edm::InputTag GenJetInputTag;

  //begin sven stuff for muons
  std::string getFloatPointDataWord(const l1t::RegionalMuonCand& l1mu) const;
  std::string getGlobalPhi(const l1t::RegionalMuonCand& l1mu) const;

  edm::EDGetToken m_emtfToken;
  edm::EDGetToken m_bmtfToken;
  edm::EDGetToken m_omtfToken;
  //end sven stuff for muons

  edm::ParameterSet simTrackCfg;

  edm::EDGetTokenT< edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > > > ttClusterToken_;
  edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > ttStubToken_;
  edm::EDGetTokenT< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > > ttClusterMCTruthToken_;
  edm::EDGetTokenT< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > ttStubMCTruthToken_;

  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_;
  edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > ttTrackMCTruthToken_;

  edm::EDGetTokenT< std::vector< TrackingParticle > > TrackingParticleToken_;
  edm::EDGetTokenT< std::vector< TrackingVertex > > TrackingVertexToken_;

  edm::EDGetTokenT< std::vector<reco::GenJet> > GenJetToken_;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticleInput_;
  edm::EDGetTokenT<edm::SimVertexContainer> simVertexInput_;
  edm::EDGetTokenT<edm::SimTrackContainer> simTrackInput_;

  edm::EDGetTokenT<edm::PSimHitContainer> gemSimHitInput_;
  edm::EDGetTokenT<edm::PSimHitContainer> cscSimHitInput_;
  edm::EDGetTokenT<edm::PSimHitContainer> rpcSimHitInput_;
  edm::EDGetTokenT<edm::PSimHitContainer> me0SimHitInput_;
  edm::EDGetTokenT<edm::PSimHitContainer> dtSimHitInput_;

  edm::EDGetTokenT<GEMDigiCollection> gemDigiInput_;
  edm::EDGetTokenT<GEMPadDigiCollection> gemPadDigiInput_;
  edm::EDGetTokenT<GEMCoPadDigiCollection> gemCoPadDigiInput_;
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHitInput_;

  edm::EDGetTokenT<ME0DigiPreRecoCollection> me0DigiInput_;
  edm::EDGetTokenT<ME0RecHitCollection> me0RecHitInput_;
  edm::EDGetTokenT<ME0SegmentCollection> me0SegmentInput_;

  edm::EDGetTokenT<CSCComparatorDigiCollection> cscComparatorDigiInput_;
  edm::EDGetTokenT<CSCWireDigiCollection> cscWireDigiInput_;
  edm::EDGetTokenT<CSCCLCTDigiCollection> clctInputs_;
  edm::EDGetTokenT<CSCALCTDigiCollection> alctInputs_;
  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> lctInputs_;
  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> mplctInputs_;
  edm::EDGetTokenT<CSCRecHit2DCollection> cscRecHit2DInput_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegmentInput_;

  edm::EDGetTokenT<DTDigiCollection> dtDigiInput_;
  edm::EDGetTokenT<DTLocalTriggerCollection> dtStubInput_;
  edm::EDGetTokenT<DTRecHitCollection> dtRecHit1DPairInput_;
  edm::EDGetTokenT<DTRecSegment2DCollection> dtRecSegment2DInput_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dtRecSegment4DInput_;

  edm::EDGetTokenT<RPCDigiCollection> rpcDigiInput_;
  edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitInput_;

  edm::EDGetTokenT<l1t::EMTFTrackCollection> emtfTrackInputLabel_;
  edm::EDGetTokenT<BXVector<l1t::RegionalMuonCand>> regMuonCandInputLabel_;
  edm::EDGetTokenT<BXVector<l1t::Muon>> gmtInputLabel_;

  //edm::EDGetTokenT<L1TTTrackCollectionType> trackInputLabel_;
  //edm::EDGetTokenT<l1t::L1TkMuonParticleCollection> trackMuonInputLabel_;

  edm::EDGetTokenT<reco::TrackExtraCollection> recoTrackExtraInputLabel_;
  edm::EDGetTokenT<reco::TrackCollection> recoTrackInputLabel_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> recoChargedCandidateInputLabel_;


  //-----------------------------------------------------------------------------------------------
  // tree & branches for mini-ntuple

  TTree* eventTree;

  // all L1 tracks
  std::vector<float>* m_trk_pt;
  std::vector<float>* m_trk_eta;
  std::vector<float>* m_trk_phi;
  std::vector<int>*   m_trk_charge;
  std::vector<float>* m_trk_d0;   // (filled if L1Tk_nPar==5, else 999)
  std::vector<float>* m_trk_z0;
  std::vector<float>* m_trk_chi2;
  std::vector<int>*   m_trk_nstub;
  std::vector<int>*   m_trk_seed;
  std::vector<int>*   m_trk_genuine;
  std::vector<int>*   m_trk_loose;
  std::vector<int>*   m_trk_unknown;
  std::vector<int>*   m_trk_combinatoric;
  std::vector<int>*   m_trk_fake; //0 fake, 1 track from primary interaction, 2 secondary track
  std::vector<int>*   m_trk_matchtp_pdgid;
  std::vector<float>* m_trk_matchtp_pt;
  std::vector<float>* m_trk_matchtp_eta;
  std::vector<float>* m_trk_matchtp_phi;
  std::vector<float>* m_trk_matchtp_z0;
  std::vector<float>* m_trk_matchtp_dxy;
  std::vector<int>*   m_trk_injet;         //is the track within dR<0.4 of a genjet with pt > 30 GeV?
  std::vector<int>*   m_trk_injet_highpt;  //is the track within dR<0.4 of a genjet with pt > 100 GeV?
  std::vector<int>*   m_trk_injet_vhighpt; //is the track within dR<0.4 of a genjet with pt > 200 GeV?

  // all tracking particles
  std::vector<float>* m_tp_pt;
  std::vector<float>* m_tp_eta;
  std::vector<float>* m_tp_phi;
  std::vector<float>* m_tp_dxy;
  std::vector<float>* m_tp_d0;
  std::vector<float>* m_tp_z0;
  std::vector<float>* m_tp_d0_prod;
  std::vector<float>* m_tp_z0_prod;
  std::vector<int>*   m_tp_pdgid;
  std::vector<int>*   m_tp_nmatch;
  std::vector<int>*   m_tp_nloosematch;
  std::vector<int>*   m_tp_nstub;
  std::vector<int>*   m_tp_eventid;
  std::vector<int>*   m_tp_charge;
  std::vector<int>*   m_tp_injet;
  std::vector<int>*   m_tp_injet_highpt;
  std::vector<int>*   m_tp_injet_vhighpt;

  // *L1 track* properties if m_tp_nmatch > 0
  std::vector<float>* m_matchtrk_pt;
  std::vector<float>* m_matchtrk_eta;
  std::vector<float>* m_matchtrk_phi;
  std::vector<float>* m_matchtrk_d0; //this variable is only filled if L1Tk_nPar==5
  std::vector<float>* m_matchtrk_z0;
  std::vector<float>* m_matchtrk_chi2;
  std::vector<int>*   m_matchtrk_nstub;
  std::vector<int>*   m_matchtrk_seed;
  std::vector<int>*   m_matchtrk_injet;
  std::vector<int>*   m_matchtrk_injet_highpt;
  std::vector<int>*   m_matchtrk_injet_vhighpt;
  std::vector<bool>*  m_matchtrk_charge;

  // *L1 track* properties if m_tp_nloosematch > 0
  std::vector<float>* m_loosematchtrk_pt;
  std::vector<float>* m_loosematchtrk_eta;
  std::vector<float>* m_loosematchtrk_phi;
  std::vector<float>* m_loosematchtrk_d0; //this variable is only filled if L1Tk_nPar==5
  std::vector<float>* m_loosematchtrk_z0;
  std::vector<float>* m_loosematchtrk_chi2;
  std::vector<int>*   m_loosematchtrk_nstub;
  std::vector<int>*   m_loosematchtrk_seed;
  std::vector<int>*   m_loosematchtrk_injet;
  std::vector<int>*   m_loosematchtrk_injet_highpt;
  std::vector<int>*   m_loosematchtrk_injet_vhighpt;

  // ALL stubs
  std::vector<float>* m_allstub_x;
  std::vector<float>* m_allstub_y;
  std::vector<float>* m_allstub_z;

  std::vector<int>*   m_allstub_isBarrel; // stub is in barrel (1) or in disk (0)
  std::vector<int>*   m_allstub_layer;
  std::vector<int>*   m_allstub_isPSmodule;

  std::vector<float>* m_allstub_trigDisplace;
  std::vector<float>* m_allstub_trigOffset;
  std::vector<float>* m_allstub_trigPos;
  std::vector<float>* m_allstub_trigBend;

  // stub associated with tracking particle ?
  std::vector<int>*   m_allstub_matchTP_pdgid; // -999 if not matched
  std::vector<float>* m_allstub_matchTP_pt;    // -999 if not matched
  std::vector<float>* m_allstub_matchTP_eta;   // -999 if not matched
  std::vector<float>* m_allstub_matchTP_phi;   // -999 if not matched

  std::vector<int>*   m_allstub_genuine;

  // track jet variables (for each gen jet, store the sum of pt of TPs / tracks inside jet cone)
  std::vector<float>* m_jet_eta;
  std::vector<float>* m_jet_phi;
  std::vector<float>* m_jet_pt;
  std::vector<float>* m_jet_tp_sumpt;
  std::vector<float>* m_jet_trk_sumpt;
  std::vector<float>* m_jet_matchtrk_sumpt;
  std::vector<float>* m_jet_loosematchtrk_sumpt;

  // sven muon variables
  std::vector<int>*   m_EMTF_muon_n;
  std::vector<float>* m_EMTF_muon_pt;
  std::vector<float>* m_EMTF_muon_eta;
  std::vector<float>* m_EMTF_muon_phi;
  std::vector<int>*   m_EMTF_muon_c;

  std::vector<int>*   m_OMTF_muon_n;
  std::vector<float>* m_OMTF_muon_pt;
  std::vector<float>* m_OMTF_muon_eta;
  std::vector<float>* m_OMTF_muon_phi;
  std::vector<int>*   m_OMTF_muon_c;

  std::vector<int>*   m_BMTF_muon_n;
  std::vector<float>* m_BMTF_muon_pt;
  std::vector<float>* m_BMTF_muon_eta;
  std::vector<float>* m_BMTF_muon_phi;
  std::vector<int>*   m_BMTF_muon_c;

  std::vector<float>* m_matchmuon_pt;
  std::vector<float>* m_matchmuon_eta;
  std::vector<float>* m_matchmuon_phi;
  std::vector<int>*   m_matchmuon_charge;
  std::vector<int>*   m_matchmuon_type;
  std::vector<int>*   m_matchmuon_quality;

};


//////////////////////////////////
//                              //
//     CLASS IMPLEMENTATION     //
//                              //
//////////////////////////////////

//////////////
// CONSTRUCTOR
L1TrackNtupleMaker::L1TrackNtupleMaker(edm::ParameterSet const& iConfig) :
  config(iConfig)
{

  MyProcess        = iConfig.getParameter< int >("MyProcess");
  DebugMode        = iConfig.getParameter< bool >("DebugMode");
  SaveAllTracks    = iConfig.getParameter< bool >("SaveAllTracks");
  SaveStubs        = iConfig.getParameter< bool >("SaveStubs");
  SaveMuons	   = iConfig.getParameter< bool >("SaveMuons");
  L1Tk_nPar        = iConfig.getParameter< int >("L1Tk_nPar");
  TP_minNStub      = iConfig.getParameter< int >("TP_minNStub");
  TP_minNStubLayer = iConfig.getParameter< int >("TP_minNStubLayer");
  TP_minPt         = iConfig.getParameter< double >("TP_minPt");
  TP_maxEta        = iConfig.getParameter< double >("TP_maxEta");
  TP_maxZ0         = iConfig.getParameter< double >("TP_maxZ0");
  L1TrackInputTag      = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
  MCTruthTrackInputTag = iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag");
  L1Tk_minNStub        = iConfig.getParameter< int >("L1Tk_minNStub");

  TrackingInJets = iConfig.getParameter< bool >("TrackingInJets");

  L1StubInputTag           = iConfig.getParameter<edm::InputTag>("L1StubInputTag");
  MCTruthClusterInputTag   = iConfig.getParameter<edm::InputTag>("MCTruthClusterInputTag");
  MCTruthStubInputTag      = iConfig.getParameter<edm::InputTag>("MCTruthStubInputTag");
  TrackingParticleInputTag = iConfig.getParameter<edm::InputTag>("TrackingParticleInputTag");
  TrackingVertexInputTag   = iConfig.getParameter<edm::InputTag>("TrackingVertexInputTag");
  GenJetInputTag           = iConfig.getParameter<edm::InputTag>("GenJetInputTag");

  ttTrackToken_          = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(L1TrackInputTag);
  ttTrackMCTruthToken_   = consumes< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthTrackInputTag);
  ttStubToken_           = consumes< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >(L1StubInputTag);
  ttClusterMCTruthToken_ = consumes< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthClusterInputTag);
  ttStubMCTruthToken_    = consumes< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthStubInputTag);

  TrackingParticleToken_ = consumes< std::vector< TrackingParticle > >(TrackingParticleInputTag);
  TrackingVertexToken_   = consumes< std::vector< TrackingVertex > >(TrackingVertexInputTag);
  GenJetToken_           = consumes< std::vector< reco::GenJet > >(GenJetInputTag);

  //Svens muon consumes
  m_emtfToken = consumes<l1t::RegionalMuonCandBxCollection>(edm::InputTag("simEmtfDigis","EMTF"));
  m_bmtfToken = consumes<l1t::RegionalMuonCandBxCollection>(edm::InputTag("simBmtfDigis","BMTF"));
  m_omtfToken = consumes<l1t::RegionalMuonCandBxCollection>(edm::InputTag("simOmtfDigis","OMTF"));

  // simtrack matching
  simTrackCfg = (iConfig.getParameterSet("simTrackMatching"));

  const auto& displacedGenMu = simTrackCfg.getParameterSet("displacedGenMu");
  genParticleInput_ = consumes<reco::GenParticleCollection>(displacedGenMu.getParameter<edm::InputTag>("inputTag"));

  const auto& simVertex = simTrackCfg.getParameterSet("simVertex");
  simVertexInput_ = consumes<edm::SimVertexContainer>(simVertex.getParameter<edm::InputTag>("inputTag"));

  const auto& simTrack = simTrackCfg.getParameterSet("simTrack");
  simTrackInput_ = consumes<edm::SimTrackContainer>(simTrack.getParameter<edm::InputTag>("inputTag"));

  const auto& gemSimHit_ = simTrackCfg.getParameterSet("gemSimHit");
  gemSimHitInput_ = consumes<edm::PSimHitContainer>(gemSimHit_.getParameter<edm::InputTag>("inputTag"));

  const auto& cscSimHit_= simTrackCfg.getParameterSet("cscSimHit");
  cscSimHitInput_ = consumes<edm::PSimHitContainer>(cscSimHit_.getParameter<edm::InputTag>("inputTag"));

  const auto& me0SimHit_ = simTrackCfg.getParameterSet("me0SimHit");
  me0SimHitInput_ = consumes<edm::PSimHitContainer>(me0SimHit_.getParameter<edm::InputTag>("inputTag"));

  const auto& rpcSimHit_ = simTrackCfg.getParameterSet("rpcSimHit");
  rpcSimHitInput_ = consumes<edm::PSimHitContainer>(rpcSimHit_.getParameter<edm::InputTag>("inputTag"));

  const auto& dtSimHit_ = simTrackCfg.getParameterSet("dtSimHit");
  dtSimHitInput_ = consumes<edm::PSimHitContainer>(dtSimHit_.getParameter<edm::InputTag>("inputTag"));

  const auto& gemDigi_= simTrackCfg.getParameterSet("gemStripDigi");
  gemDigiInput_ = consumes<GEMDigiCollection>(gemDigi_.getParameter<edm::InputTag>("inputTag"));

  const auto& gemPad_= simTrackCfg.getParameterSet("gemPadDigi");
  gemPadDigiInput_ = consumes<GEMPadDigiCollection>(gemPad_.getParameter<edm::InputTag>("inputTag"));

  const auto& gemCoPad_= simTrackCfg.getParameterSet("gemCoPadDigi");
  gemCoPadDigiInput_ = consumes<GEMCoPadDigiCollection>(gemCoPad_.getParameter<edm::InputTag>("inputTag"));

  const auto& gemRecHit_= simTrackCfg.getParameterSet("gemRecHit");
  gemRecHitInput_ = consumes<GEMRecHitCollection>(gemRecHit_.getParameter<edm::InputTag>("inputTag"));

  const auto& me0Digi_= simTrackCfg.getParameterSet("me0DigiPreReco");
  me0DigiInput_ = consumes<ME0DigiPreRecoCollection>(me0Digi_.getParameter<edm::InputTag>("inputTag"));

  const auto& me0RecHit_ = simTrackCfg.getParameterSet("me0RecHit");
  me0RecHitInput_ = consumes<ME0RecHitCollection>(me0RecHit_.getParameter<edm::InputTag>("inputTag"));

  const auto& me0Segment_ = simTrackCfg.getParameterSet("me0Segment");
  me0SegmentInput_ = consumes<ME0SegmentCollection>(me0Segment_.getParameter<edm::InputTag>("inputTag"));

  const auto& cscComparatorDigi = simTrackCfg.getParameterSet("cscStripDigi");
  cscComparatorDigiInput_ = consumes<CSCComparatorDigiCollection>(cscComparatorDigi.getParameter<edm::InputTag>("inputTag"));

  const auto& cscWireDigi = simTrackCfg.getParameterSet("cscWireDigi");
  cscWireDigiInput_ = consumes<CSCWireDigiCollection>(cscWireDigi.getParameter<edm::InputTag>("inputTag"));

  const auto& cscCLCT = simTrackCfg.getParameterSet("cscCLCT");
  clctInputs_ = consumes<CSCCLCTDigiCollection>(cscCLCT.getParameter<edm::InputTag>("inputTag"));

  const auto& cscALCT = simTrackCfg.getParameterSet("cscALCT");
  alctInputs_ = consumes<CSCALCTDigiCollection>(cscALCT.getParameter<edm::InputTag>("inputTag"));

  const auto& cscLCT = simTrackCfg.getParameterSet("cscLCT");
  lctInputs_ = consumes<CSCCorrelatedLCTDigiCollection>(cscLCT.getParameter<edm::InputTag>("inputTag"));

  const auto& cscMPLCT = simTrackCfg.getParameterSet("cscMPLCT");
  mplctInputs_ = consumes<CSCCorrelatedLCTDigiCollection>(cscMPLCT.getParameter<edm::InputTag>("inputTag"));

  const auto& cscRecHit2D = simTrackCfg.getParameterSet("cscRecHit");
  cscRecHit2DInput_ = consumes<CSCRecHit2DCollection>(cscRecHit2D.getParameter<edm::InputTag>("inputTag"));

  const auto& cscSegment2D = simTrackCfg.getParameterSet("cscSegment");
  cscSegmentInput_ = consumes<CSCSegmentCollection>(cscSegment2D.getParameter<edm::InputTag>("inputTag"));

  const auto& dtDigi_= simTrackCfg.getParameterSet("dtDigi");
  dtDigiInput_ = consumes<DTDigiCollection>(dtDigi_.getParameter<edm::InputTag>("inputTag"));

  const auto& dtStub_= simTrackCfg.getParameterSet("dtLocalTrigger");
  dtStubInput_ = consumes<DTLocalTriggerCollection>(dtStub_.getParameter<edm::InputTag>("inputTag"));

  const auto& dtRecHit1DPair = simTrackCfg.getParameterSet("dtRecHit");
  dtRecHit1DPairInput_ = consumes<DTRecHitCollection>(dtRecHit1DPair.getParameter<edm::InputTag>("inputTag"));

  const auto& dtSegment2D = simTrackCfg.getParameterSet("dtRecSegment2D");
  dtRecSegment2DInput_ = consumes<DTRecSegment2DCollection>(dtSegment2D.getParameter<edm::InputTag>("inputTag"));

  const auto& dtSegment4D = simTrackCfg.getParameterSet("dtRecSegment4D");
  dtRecSegment4DInput_ = consumes<DTRecSegment4DCollection>(dtSegment4D.getParameter<edm::InputTag>("inputTag"));

  const auto& rpcDigi_= simTrackCfg.getParameterSet("rpcStripDigi");
  rpcDigiInput_ = consumes<RPCDigiCollection>(rpcDigi_.getParameter<edm::InputTag>("inputTag"));

  const auto& rpcRecHit_= simTrackCfg.getParameterSet("rpcRecHit");
  rpcRecHitInput_ = consumes<RPCRecHitCollection>(rpcRecHit_.getParameter<edm::InputTag>("inputTag"));

  const auto& emtfTrack = simTrackCfg.getParameterSet("upgradeEmtfTrack");
  emtfTrackInputLabel_ = consumes<l1t::EMTFTrackCollection>(emtfTrack.getParameter<edm::InputTag>("inputTag"));

  const auto& upgradeemtfCand = simTrackCfg.getParameterSet("upgradeEmtfCand");
  regMuonCandInputLabel_ = consumes< BXVector<l1t::RegionalMuonCand> >(upgradeemtfCand.getParameter<edm::InputTag>("inputTag"));

  const auto& upgradegmt = simTrackCfg.getParameterSet("upgradeGMT");
  gmtInputLabel_ = consumes< BXVector<l1t::Muon> >(upgradegmt.getParameter<edm::InputTag>("inputTag"));

  //const auto& l1Track = simTrackCfg.getParameterSet("l1track");
  //trackInputLabel_ = consumes<L1TTTrackCollectionType>(l1Track.getParameter<edm::InputTag>("inputTag"));
  //verboseL1Track_ = l1Track.getParameter<int>("verbose");

  //const auto& l1TrackMuon = simTrackCfg.getParameterSet("l1tkmuon");
  //trackMuonInputLabel_ = consumes<l1t::L1TkMuonParticleCollection>(l1TrackMuon.getParameter<edm::InputTag>("inputTag"));

  const auto& recoTrackExtra = simTrackCfg.getParameterSet("recoTrackExtra");
  recoTrackExtraInputLabel_ = consumes<reco::TrackExtraCollection>(recoTrackExtra.getParameter<edm::InputTag>("inputTag"));

  const auto& recoTrack = simTrackCfg.getParameterSet("recoTrack");
  recoTrackInputLabel_ = consumes<reco::TrackCollection>(recoTrack.getParameter<edm::InputTag>("inputTag"));

  const auto& recoChargedCandidate = simTrackCfg.getParameterSet("recoChargedCandidate");
  recoChargedCandidateInputLabel_ = consumes<reco::RecoChargedCandidateCollection>(recoChargedCandidate.getParameter<edm::InputTag>("inputTag"));

}

/////////////
// DESTRUCTOR
L1TrackNtupleMaker::~L1TrackNtupleMaker()
{
}

//////////
// END JOB
void L1TrackNtupleMaker::endJob()
{
  // things to be done at the exit of the event Loop
  cerr << "L1TrackNtupleMaker::endJob" << endl;

}

////////////
// BEGIN JOB
void L1TrackNtupleMaker::beginJob()
{

  // things to be done before entering the event Loop
  cerr << "L1TrackNtupleMaker::beginJob" << endl;

  //-----------------------------------------------------------------------------------------------
  // book histograms / make ntuple
  edm::Service<TFileService> fs;

  SaveTracklet = true;

  // initilize
  m_trk_pt    = new std::vector<float>;
  m_trk_eta   = new std::vector<float>;
  m_trk_phi   = new std::vector<float>;
  m_trk_charge= new std::vector<int>;
  m_trk_z0    = new std::vector<float>;
  m_trk_d0    = new std::vector<float>;
  m_trk_chi2  = new std::vector<float>;
  m_trk_nstub = new std::vector<int>;
  m_trk_seed    = new std::vector<int>;
  m_trk_genuine       = new std::vector<int>;
  m_trk_loose         = new std::vector<int>;
  m_trk_unknown       = new std::vector<int>;
  m_trk_combinatoric  = new std::vector<int>;
  m_trk_fake          = new std::vector<int>;
  m_trk_matchtp_pdgid = new std::vector<int>;
  m_trk_matchtp_pt    = new std::vector<float>;
  m_trk_matchtp_eta   = new std::vector<float>;
  m_trk_matchtp_phi   = new std::vector<float>;
  m_trk_matchtp_z0    = new std::vector<float>;
  m_trk_matchtp_dxy   = new std::vector<float>;
  m_trk_injet         = new std::vector<int>;
  m_trk_injet_highpt  = new std::vector<int>;
  m_trk_injet_vhighpt  = new std::vector<int>;

  m_tp_pt      = new std::vector<float>;
  m_tp_eta     = new std::vector<float>;
  m_tp_phi     = new std::vector<float>;
  m_tp_dxy     = new std::vector<float>;
  m_tp_d0      = new std::vector<float>;
  m_tp_z0      = new std::vector<float>;
  m_tp_d0_prod = new std::vector<float>;
  m_tp_z0_prod = new std::vector<float>;
  m_tp_pdgid   = new std::vector<int>;
  m_tp_nmatch  = new std::vector<int>;
  m_tp_nloosematch  = new std::vector<int>;
  m_tp_nstub        = new std::vector<int>;
  m_tp_eventid      = new std::vector<int>;
  m_tp_charge       = new std::vector<int>;
  m_tp_injet        = new std::vector<int>;
  m_tp_injet_highpt = new std::vector<int>;
  m_tp_injet_vhighpt = new std::vector<int>;

  m_matchtrk_pt    = new std::vector<float>;
  m_matchtrk_eta   = new std::vector<float>;
  m_matchtrk_phi   = new std::vector<float>;
  m_matchtrk_z0    = new std::vector<float>;
  m_matchtrk_d0    = new std::vector<float>;
  m_matchtrk_chi2  = new std::vector<float>;
  m_matchtrk_nstub = new std::vector<int>;
  m_matchtrk_seed  = new std::vector<int>;
  m_matchtrk_injet = new std::vector<int>;
  m_matchtrk_injet_highpt = new std::vector<int>;
  m_matchtrk_injet_vhighpt = new std::vector<int>;
  m_matchtrk_charge = new std::vector<bool>;

  m_loosematchtrk_pt    = new std::vector<float>;
  m_loosematchtrk_eta   = new std::vector<float>;
  m_loosematchtrk_phi   = new std::vector<float>;
  m_loosematchtrk_z0    = new std::vector<float>;
  m_loosematchtrk_d0    = new std::vector<float>;
  m_loosematchtrk_chi2  = new std::vector<float>;
  m_loosematchtrk_nstub = new std::vector<int>;
  m_loosematchtrk_seed  = new std::vector<int>;
  m_loosematchtrk_injet = new std::vector<int>;
  m_loosematchtrk_injet_highpt = new std::vector<int>;
  m_loosematchtrk_injet_vhighpt = new std::vector<int>;

  m_matchmuon_pt     = new std::vector<float>;
  m_matchmuon_eta    = new std::vector<float>;
  m_matchmuon_phi    = new std::vector<float>;
  m_matchmuon_charge = new std::vector<int>;
  m_matchmuon_type   = new std::vector<int>;
  m_matchmuon_quality= new std::vector<int>;

  m_allstub_x = new std::vector<float>;
  m_allstub_y = new std::vector<float>;
  m_allstub_z = new std::vector<float>;

  m_allstub_isBarrel = new std::vector<int>;
  m_allstub_layer    = new std::vector<int>;
  m_allstub_isPSmodule   = new std::vector<int>;
  m_allstub_trigDisplace = new std::vector<float>;
  m_allstub_trigOffset   = new std::vector<float>;
  m_allstub_trigPos      = new std::vector<float>;
  m_allstub_trigBend     = new std::vector<float>;

  m_allstub_matchTP_pdgid = new std::vector<int>;
  m_allstub_matchTP_pt    = new std::vector<float>;
  m_allstub_matchTP_eta   = new std::vector<float>;
  m_allstub_matchTP_phi   = new std::vector<float>;

  m_allstub_genuine = new std::vector<int>;

  m_jet_eta = new std::vector<float>;
  m_jet_phi = new std::vector<float>;
  m_jet_pt  = new std::vector<float>;
  m_jet_tp_sumpt  = new std::vector<float>;
  m_jet_trk_sumpt = new std::vector<float>;
  m_jet_matchtrk_sumpt      = new std::vector<float>;
  m_jet_loosematchtrk_sumpt = new std::vector<float>;

  //extra sven muon stuff
  m_EMTF_muon_n = new std::vector<int>;
  m_EMTF_muon_pt = new std::vector<float>;
  m_EMTF_muon_eta = new std::vector<float>;
  m_EMTF_muon_phi = new std::vector<float>;
  m_EMTF_muon_c = new std::vector<int>;

  m_OMTF_muon_n = new std::vector<int>;
  m_OMTF_muon_pt = new std::vector<float>;
  m_OMTF_muon_eta = new std::vector<float>;
  m_OMTF_muon_phi = new std::vector<float>;
  m_OMTF_muon_c = new std::vector<int>;

  m_BMTF_muon_n = new std::vector<int>;
  m_BMTF_muon_pt = new std::vector<float>;
  m_BMTF_muon_eta = new std::vector<float>;
  m_BMTF_muon_phi = new std::vector<float>;
  m_BMTF_muon_c = new std::vector<int>;


  // ntuple
  eventTree = fs->make<TTree>("eventTree", "Event tree");

  if (SaveAllTracks) {
    eventTree->Branch("trk_pt",    &m_trk_pt);
    eventTree->Branch("trk_eta",   &m_trk_eta);
    eventTree->Branch("trk_phi",   &m_trk_phi);
    eventTree->Branch("trk_charge",   &m_trk_charge);
    eventTree->Branch("trk_d0",    &m_trk_d0);
    eventTree->Branch("trk_z0",    &m_trk_z0);
    eventTree->Branch("trk_chi2",  &m_trk_chi2);
    eventTree->Branch("trk_nstub", &m_trk_nstub);
    if (SaveTracklet) eventTree->Branch("trk_seed",    &m_trk_seed);
    eventTree->Branch("trk_genuine",      &m_trk_genuine);
    eventTree->Branch("trk_loose",        &m_trk_loose);
    eventTree->Branch("trk_unknown",      &m_trk_unknown);
    eventTree->Branch("trk_combinatoric", &m_trk_combinatoric);
    eventTree->Branch("trk_fake",         &m_trk_fake);
    eventTree->Branch("trk_matchtp_pdgid",&m_trk_matchtp_pdgid);
    eventTree->Branch("trk_matchtp_pt",   &m_trk_matchtp_pt);
    eventTree->Branch("trk_matchtp_eta",  &m_trk_matchtp_eta);
    eventTree->Branch("trk_matchtp_phi",  &m_trk_matchtp_phi);
    eventTree->Branch("trk_matchtp_z0",   &m_trk_matchtp_z0);
    eventTree->Branch("trk_matchtp_dxy",  &m_trk_matchtp_dxy);
    if (TrackingInJets) {
      eventTree->Branch("trk_injet",         &m_trk_injet);
      eventTree->Branch("trk_injet_highpt",  &m_trk_injet_highpt);
      eventTree->Branch("trk_injet_vhighpt", &m_trk_injet_vhighpt);
    }
  }

  eventTree->Branch("tp_pt",     &m_tp_pt);
  eventTree->Branch("tp_eta",    &m_tp_eta);
  eventTree->Branch("tp_phi",    &m_tp_phi);
  eventTree->Branch("tp_dxy",    &m_tp_dxy);
  eventTree->Branch("tp_d0",     &m_tp_d0);
  eventTree->Branch("tp_z0",     &m_tp_z0);
  eventTree->Branch("tp_d0_prod",     &m_tp_d0_prod);
  eventTree->Branch("tp_z0_prod",     &m_tp_z0_prod);
  eventTree->Branch("tp_pdgid",       &m_tp_pdgid);
  eventTree->Branch("tp_nmatch",      &m_tp_nmatch);
  eventTree->Branch("tp_nloosematch", &m_tp_nloosematch);
  eventTree->Branch("tp_nstub",       &m_tp_nstub);
  eventTree->Branch("tp_eventid",     &m_tp_eventid);
  eventTree->Branch("tp_charge",      &m_tp_charge);
  if (TrackingInJets) {
    eventTree->Branch("tp_injet",         &m_tp_injet);
    eventTree->Branch("tp_injet_highpt",  &m_tp_injet_highpt);
    eventTree->Branch("tp_injet_vhighpt", &m_tp_injet_vhighpt);
  }

  eventTree->Branch("matchtrk_pt",    &m_matchtrk_pt);
  eventTree->Branch("matchtrk_eta",   &m_matchtrk_eta);
  eventTree->Branch("matchtrk_phi",   &m_matchtrk_phi);
  eventTree->Branch("matchtrk_z0",    &m_matchtrk_z0);
  eventTree->Branch("matchtrk_d0",    &m_matchtrk_d0);
  eventTree->Branch("matchtrk_chi2",  &m_matchtrk_chi2);
  eventTree->Branch("matchtrk_nstub", &m_matchtrk_nstub);
  if (SaveTracklet) eventTree->Branch("matchtrk_seed",    &m_matchtrk_seed);
  if (TrackingInJets) {
    eventTree->Branch("matchtrk_injet",         &m_matchtrk_injet);
    eventTree->Branch("matchtrk_injet_highpt",  &m_matchtrk_injet_highpt);
    eventTree->Branch("matchtrk_injet_vhighpt", &m_matchtrk_injet_vhighpt);
  }
  eventTree->Branch("matchtrk_charge",&m_matchtrk_charge);

  eventTree->Branch("loosematchtrk_pt",    &m_loosematchtrk_pt);
  eventTree->Branch("loosematchtrk_eta",   &m_loosematchtrk_eta);
  eventTree->Branch("loosematchtrk_phi",   &m_loosematchtrk_phi);
  eventTree->Branch("loosematchtrk_z0",    &m_loosematchtrk_z0);
  eventTree->Branch("loosematchtrk_d0",    &m_loosematchtrk_d0);
  eventTree->Branch("loosematchtrk_chi2",  &m_loosematchtrk_chi2);
  eventTree->Branch("loosematchtrk_nstub", &m_loosematchtrk_nstub);
  if (SaveTracklet) eventTree->Branch("loosematchtrk_seed",    &m_loosematchtrk_seed);
  if (TrackingInJets) {
    eventTree->Branch("loosematchtrk_injet",         &m_loosematchtrk_injet);
    eventTree->Branch("loosematchtrk_injet_highpt",  &m_loosematchtrk_injet_highpt);
    eventTree->Branch("loosematchtrk_injet_vhighpt", &m_loosematchtrk_injet_vhighpt);
  }

  if (SaveMuons){
    eventTree->Branch("matchmuon_pt", &m_matchmuon_pt);
    eventTree->Branch("matchmuon_eta", &m_matchmuon_eta);
    eventTree->Branch("matchmuon_phi", &m_matchmuon_phi);
    eventTree->Branch("matchmuon_charge",&m_matchmuon_charge);
    eventTree->Branch("matchmuon_type",&m_matchmuon_type);
    eventTree->Branch("matchmuon_quality",&m_matchmuon_quality);
  }

  if (SaveStubs) {
    eventTree->Branch("allstub_x", &m_allstub_x);
    eventTree->Branch("allstub_y", &m_allstub_y);
    eventTree->Branch("allstub_z", &m_allstub_z);

    eventTree->Branch("allstub_isBarrel",   &m_allstub_isBarrel);
    eventTree->Branch("allstub_layer",      &m_allstub_layer);
    eventTree->Branch("allstub_isPSmodule", &m_allstub_isPSmodule);

    eventTree->Branch("allstub_trigDisplace", &m_allstub_trigDisplace);
    eventTree->Branch("allstub_trigOffset",   &m_allstub_trigOffset);
    eventTree->Branch("allstub_trigPos",      &m_allstub_trigPos);
    eventTree->Branch("allstub_trigBend",     &m_allstub_trigBend);

    eventTree->Branch("allstub_matchTP_pdgid", &m_allstub_matchTP_pdgid);
    eventTree->Branch("allstub_matchTP_pt",    &m_allstub_matchTP_pt);
    eventTree->Branch("allstub_matchTP_eta",   &m_allstub_matchTP_eta);
    eventTree->Branch("allstub_matchTP_phi",   &m_allstub_matchTP_phi);

    eventTree->Branch("allstub_genuine", &m_allstub_genuine);
  }

  if (TrackingInJets) {
    eventTree->Branch("jet_eta", &m_jet_eta);
    eventTree->Branch("jet_phi", &m_jet_phi);
    eventTree->Branch("jet_pt",  &m_jet_pt);
    eventTree->Branch("jet_tp_sumpt",  &m_jet_tp_sumpt);
    eventTree->Branch("jet_trk_sumpt", &m_jet_trk_sumpt);
    eventTree->Branch("jet_matchtrk_sumpt",      &m_jet_matchtrk_sumpt);
    eventTree->Branch("jet_loosematchtrk_sumpt", &m_jet_loosematchtrk_sumpt);
  }

  //extra sven muon stuff
  if (SaveMuons) {
    eventTree->Branch("EMTF_muon_n",	 &m_EMTF_muon_n);
    eventTree->Branch("EMTF_muon_pt",  	 &m_EMTF_muon_pt);
    eventTree->Branch("EMTF_muon_eta", 	 &m_EMTF_muon_eta);
    eventTree->Branch("EMTF_muon_phi", 	 &m_EMTF_muon_phi);
    eventTree->Branch("EMTF_muon_c", 	 &m_EMTF_muon_c);

    eventTree->Branch("OMTF_muon_n",	 &m_OMTF_muon_n);
    eventTree->Branch("OMTF_muon_pt", 	 &m_OMTF_muon_pt);
    eventTree->Branch("OMTF_muon_eta", 	 &m_OMTF_muon_eta);
    eventTree->Branch("OMTF_muon_phi", 	 &m_OMTF_muon_phi);
    eventTree->Branch("OMTF_muon_c", 	 &m_OMTF_muon_c);

    eventTree->Branch("BMTF_muon_n",	 &m_BMTF_muon_n);
    eventTree->Branch("BMTF_muon_pt", 	 &m_BMTF_muon_pt);
    eventTree->Branch("BMTF_muon_eta", 	 &m_BMTF_muon_eta);
    eventTree->Branch("BMTF_muon_phi", 	 &m_BMTF_muon_phi);
    eventTree->Branch("BMTF_muon_c", 	 &m_BMTF_muon_c);

  }


}


//////////
// ANALYZE
void L1TrackNtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if (!(MyProcess==13 || MyProcess==11 || MyProcess==211 || MyProcess==6 || MyProcess==15 || MyProcess==1)) {
    cout << "The specified MyProcess is invalid! Exiting..." << endl;
    return;
  }

  if ( !(L1Tk_nPar==4 || L1Tk_nPar==5) ) {
    cout << "Invalid number of track parameters, specified L1Tk_nPar == " << L1Tk_nPar << " but only 4/5 are valid options! Exiting..." << endl;
    return;
  }

  // clear variables
  if (SaveAllTracks) {
    m_trk_pt->clear();
    m_trk_eta->clear();
    m_trk_phi->clear();
    m_trk_charge->clear();
    m_trk_d0->clear();
    m_trk_z0->clear();
    m_trk_chi2->clear();
    m_trk_nstub->clear();
    m_trk_seed->clear();
    m_trk_genuine->clear();
    m_trk_loose->clear();
    m_trk_unknown->clear();
    m_trk_combinatoric->clear();
    m_trk_fake->clear();
    m_trk_matchtp_pdgid->clear();
    m_trk_matchtp_pt->clear();
    m_trk_matchtp_eta->clear();
    m_trk_matchtp_phi->clear();
    m_trk_matchtp_z0->clear();
    m_trk_matchtp_dxy->clear();
    m_trk_injet->clear();
    m_trk_injet_highpt->clear();
    m_trk_injet_vhighpt->clear();
  }

  m_tp_pt->clear();
  m_tp_eta->clear();
  m_tp_phi->clear();
  m_tp_dxy->clear();
  m_tp_d0->clear();
  m_tp_z0->clear();
  m_tp_d0_prod->clear();
  m_tp_z0_prod->clear();
  m_tp_pdgid->clear();
  m_tp_nmatch->clear();
  m_tp_nloosematch->clear();
  m_tp_nstub->clear();
  m_tp_eventid->clear();
  m_tp_charge->clear();
  m_tp_injet->clear();
  m_tp_injet_highpt->clear();
  m_tp_injet_vhighpt->clear();

  m_matchtrk_pt->clear();
  m_matchtrk_eta->clear();
  m_matchtrk_phi->clear();
  m_matchtrk_z0->clear();
  m_matchtrk_d0->clear();
  m_matchtrk_chi2->clear();
  m_matchtrk_nstub->clear();
  m_matchtrk_seed->clear();
  m_matchtrk_injet->clear();
  m_matchtrk_injet_highpt->clear();
  m_matchtrk_injet_vhighpt->clear();
  m_matchtrk_charge->clear();

  m_loosematchtrk_pt->clear();
  m_loosematchtrk_eta->clear();
  m_loosematchtrk_phi->clear();
  m_loosematchtrk_z0->clear();
  m_loosematchtrk_d0->clear();
  m_loosematchtrk_chi2->clear();
  m_loosematchtrk_nstub->clear();
  m_loosematchtrk_seed->clear();
  m_loosematchtrk_injet->clear();
  m_loosematchtrk_injet_highpt->clear();
  m_loosematchtrk_injet_vhighpt->clear();

  m_matchmuon_pt->clear();
  m_matchmuon_eta->clear();
  m_matchmuon_phi->clear();
  m_matchmuon_charge->clear();
  m_matchmuon_type->clear();
  m_matchmuon_quality->clear();

  if (SaveStubs) {
    m_allstub_x->clear();
    m_allstub_y->clear();
    m_allstub_z->clear();

    m_allstub_isBarrel->clear();
    m_allstub_layer->clear();
    m_allstub_isPSmodule->clear();

    m_allstub_trigDisplace->clear();
    m_allstub_trigOffset->clear();
    m_allstub_trigPos->clear();
    m_allstub_trigBend->clear();

    m_allstub_matchTP_pdgid->clear();
    m_allstub_matchTP_pt->clear();
    m_allstub_matchTP_eta->clear();
    m_allstub_matchTP_phi->clear();

    m_allstub_genuine->clear();
  }

  m_jet_eta->clear();
  m_jet_phi->clear();
  m_jet_pt->clear();
  m_jet_tp_sumpt->clear();
  m_jet_trk_sumpt->clear();
  m_jet_matchtrk_sumpt->clear();
  m_jet_loosematchtrk_sumpt->clear();

  //extra sven muon stuff
  if (SaveMuons) {
    m_EMTF_muon_n->clear();
    m_EMTF_muon_pt->clear();
    m_EMTF_muon_eta->clear();
    m_EMTF_muon_phi->clear();
    m_EMTF_muon_c->clear();

    m_OMTF_muon_n->clear();
    m_OMTF_muon_pt->clear();
    m_OMTF_muon_eta->clear();
    m_OMTF_muon_phi->clear();
    m_OMTF_muon_c->clear();

    m_BMTF_muon_n->clear();
    m_BMTF_muon_pt->clear();
    m_BMTF_muon_eta->clear();
    m_BMTF_muon_phi->clear();
    m_BMTF_muon_c->clear();
  }



  // -----------------------------------------------------------------------------------------------
  // retrieve various containers
  // -----------------------------------------------------------------------------------------------

  // L1 tracks
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackHandle;
  iEvent.getByToken(ttTrackToken_, TTTrackHandle);

  // L1 stubs
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > TTStubHandle;
  if (SaveStubs) iEvent.getByToken(ttStubToken_, TTStubHandle);


  // MC truth association maps
  edm::Handle< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTClusterHandle;
  iEvent.getByToken(ttClusterMCTruthToken_, MCTruthTTClusterHandle);
  edm::Handle< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTStubHandle;
  iEvent.getByToken(ttStubMCTruthToken_, MCTruthTTStubHandle);
  edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTTrackHandle;
  iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

  // tracking particles
  edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle;
  edm::Handle< std::vector< TrackingVertex > > TrackingVertexHandle;
  iEvent.getByToken(TrackingParticleToken_, TrackingParticleHandle);
  iEvent.getByToken(TrackingVertexToken_, TrackingVertexHandle);


  // -----------------------------------------------------------------------------------------------
  // more for TTStubs
  edm::ESHandle<TrackerGeometry> geometryHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);

  edm::ESHandle<TrackerTopology> tTopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);

  edm::ESHandle<TrackerGeometry> tGeomHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeomHandle);

  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();


  // ----------------------------------------------------------------------------------------------
  // loop over L1 stubs
  // ----------------------------------------------------------------------------------------------

  if (SaveStubs) {

    for (auto gd=theTrackerGeom->dets().begin(); gd != theTrackerGeom->dets().end(); gd++) {

      DetId detid = (*gd)->geographicalId();
      if(detid.subdetId()!=StripSubdetector::TOB && detid.subdetId()!=StripSubdetector::TID ) continue;
      if(!tTopo->isLower(detid) ) continue; // loop on the stacks: choose the lower arbitrarily
      DetId stackDetid = tTopo->stack(detid); // Stub module detid

      if (TTStubHandle->find( stackDetid ) == TTStubHandle->end() ) continue;

      // Get the DetSets of the Clusters
      edmNew::DetSet< TTStub< Ref_Phase2TrackerDigi_ > > stubs = (*TTStubHandle)[ stackDetid ];
      const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit( detid );
      const PixelGeomDetUnit* theGeomDet = dynamic_cast< const PixelGeomDetUnit* >( det0 );
      const PixelTopology* topol = dynamic_cast< const PixelTopology* >( &(theGeomDet->specificTopology()) );

      // loop over stubs
      for ( auto stubIter = stubs.begin();stubIter != stubs.end();++stubIter ) {
	edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_  > >, TTStub< Ref_Phase2TrackerDigi_  > >
	  tempStubPtr = edmNew::makeRefTo( TTStubHandle, stubIter );

	int isBarrel = 0;
	int layer=-999999;
	if ( detid.subdetId()==StripSubdetector::TOB ) {
	  isBarrel = 1;
	  layer  = static_cast<int>(tTopo->layer(detid));
	}
	else if ( detid.subdetId()==StripSubdetector::TID ) {
	  isBarrel = 0;
	  layer  = static_cast<int>(tTopo->layer(detid));
	}
	else {
	  cout << "WARNING -- neither TOB or TID stub, shouldn't happen..." << endl;
	  layer = -1;
	}

	int isPSmodule=0;
	if (topol->nrows() == 960) isPSmodule=1;

	MeasurementPoint coords = tempStubPtr->getClusterRef(0)->findAverageLocalCoordinatesCentered();
	LocalPoint clustlp = topol->localPosition(coords);
	GlobalPoint posStub  =  theGeomDet->surface().toGlobal(clustlp);

	double tmp_stub_x=posStub.x();
	double tmp_stub_y=posStub.y();
	double tmp_stub_z=posStub.z();

	float trigDisplace = tempStubPtr->getTriggerDisplacement();
	float trigOffset = tempStubPtr->getTriggerOffset();
	float trigPos = tempStubPtr->getTriggerPosition();
	float trigBend = tempStubPtr->getTriggerBend();

	m_allstub_x->push_back(tmp_stub_x);
	m_allstub_y->push_back(tmp_stub_y);
	m_allstub_z->push_back(tmp_stub_z);

	m_allstub_isBarrel->push_back(isBarrel);
	m_allstub_layer->push_back(layer);
	m_allstub_isPSmodule->push_back(isPSmodule);

	m_allstub_trigDisplace->push_back(trigDisplace);
	m_allstub_trigOffset->push_back(trigOffset);
	m_allstub_trigPos->push_back(trigPos);
	m_allstub_trigBend->push_back(trigBend);

	// matched to tracking particle?
	edm::Ptr< TrackingParticle > my_tp = MCTruthTTStubHandle->findTrackingParticlePtr(tempStubPtr);

	int myTP_pdgid = -999;
	float myTP_pt  = -999;
	float myTP_eta = -999;
	float myTP_phi = -999;

	if (my_tp.isNull() == false) {
	  int tmp_eventid = my_tp->eventId().event();

	  if (tmp_eventid > 0) continue; // this means stub from pileup track

	  myTP_pdgid = my_tp->pdgId();
	  myTP_pt = my_tp->p4().pt();
	  myTP_eta = my_tp->p4().eta();
	  myTP_phi = my_tp->p4().phi();
	}

	m_allstub_matchTP_pdgid->push_back(myTP_pdgid);
	m_allstub_matchTP_pt->push_back(myTP_pt);
	m_allstub_matchTP_eta->push_back(myTP_eta);
	m_allstub_matchTP_phi->push_back(myTP_phi);

	int tmp_stub_genuine = 0;
	if (MCTruthTTStubHandle->isGenuine(tempStubPtr)) tmp_stub_genuine = 1;

	m_allstub_genuine->push_back(tmp_stub_genuine);

      }

    }

  }


  // ----------------------------------------------------------------------------------------------
  // tracking in jets
  // ----------------------------------------------------------------------------------------------

  std::vector<math::XYZTLorentzVector> v_jets;
  std::vector<int> v_jets_highpt;
  std::vector<int> v_jets_vhighpt;

  if (TrackingInJets) {

    // gen jets
    if (DebugMode) cout << "get genjets" << endl;
    edm::Handle< std::vector< reco::GenJet > > GenJetHandle;
    iEvent.getByToken(GenJetToken_, GenJetHandle);


    if (GenJetHandle.isValid()) {

      if (DebugMode) cout << "loop over genjets" << endl;
      std::vector<reco::GenJet>::const_iterator iterGenJet;
      for ( iterGenJet = GenJetHandle->begin(); iterGenJet != GenJetHandle->end(); ++iterGenJet ) {

	reco::GenJet myJet = reco::GenJet(*iterGenJet);

	if (myJet.pt() < 30.0) continue;
	if (fabs(myJet.eta()) > 2.5) continue;

	if (DebugMode) cout << "genjet pt = " << myJet.pt() << ", eta = " << myJet.eta() << endl;

	bool ishighpt = false;
	bool isveryhighpt = false;
	if (myJet.pt() > 100.0)	ishighpt = true;
	if (myJet.pt() > 200.0)	isveryhighpt = true;

	math::XYZTLorentzVector jetP4 = myJet.p4();
	v_jets.push_back(jetP4);
	if (ishighpt) v_jets_highpt.push_back(1);
	else v_jets_highpt.push_back(0);
	if (isveryhighpt) v_jets_vhighpt.push_back(1);
	else v_jets_vhighpt.push_back(0);

      }// end loop over genjets
    }// end isValid

  }// end TrackingInJets

  const int NJETS = 10;
  float jets_tp_sumpt[NJETS] = {0};       //sum pt of TPs with dR<0.4 of jet
  float jets_matchtrk_sumpt[NJETS] = {0}; //sum pt of tracks matched to TP with dR<0.4 of jet
  float jets_loosematchtrk_sumpt[NJETS] = {0}; //sum pt of tracks matched to TP with dR<0.4 of jet
  float jets_trk_sumpt[NJETS] = {0};      //sum pt of all tracks with dR<0.4 of jet



  // ----------------------------------------------------------------------------------------------
  // loop over L1 tracks
  // ----------------------------------------------------------------------------------------------

  if (SaveAllTracks) {

    if (DebugMode) {
      cout << endl << "Loop over L1 tracks!" << endl;
      cout << endl << "Looking at " << L1Tk_nPar << "-parameter tracks!" << endl;
    }

    int this_l1track = 0;
    std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterL1Track;
    for ( iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {

      edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(TTTrackHandle, this_l1track);
      this_l1track++;

      float tmp_trk_pt   = iterL1Track->getMomentum(L1Tk_nPar).perp();
      float tmp_trk_eta  = iterL1Track->getMomentum(L1Tk_nPar).eta();
      float tmp_trk_phi  = iterL1Track->getMomentum(L1Tk_nPar).phi();
      int   tmp_trk_charge=signbit(iterL1Track->getRInv(L1Tk_nPar));
      float tmp_trk_z0   = iterL1Track->getPOCA(L1Tk_nPar).z(); //cm

      float tmp_trk_d0 = -999;
      if (L1Tk_nPar == 5) {
	float tmp_trk_x0   = iterL1Track->getPOCA(L1Tk_nPar).x();
	float tmp_trk_y0   = iterL1Track->getPOCA(L1Tk_nPar).y();
	tmp_trk_d0 = -tmp_trk_x0*sin(tmp_trk_phi) + tmp_trk_y0*cos(tmp_trk_phi);
      }

      float tmp_trk_chi2 = iterL1Track->getChi2(L1Tk_nPar);

      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubRefs = iterL1Track->getStubRefs();
      int tmp_trk_nstub  = (int) stubRefs.size();

      int tmp_trk_seed = 0;
      if (SaveTracklet) tmp_trk_seed = (int) iterL1Track->getWedge();

      /*
      int tmp_trk_nPSstub = 0;
      if (SaveTracklet) {
	for (int is=0; is<tmp_trk_nstub; is++) {

	  DetId detIdStub = theTrackerGeom->idToDet( (stubRefs.at(is)->getClusterRef(0))->getDetId() )->geographicalId();
	  DetId stackDetid = tTopo->stack(detIdStub);

	  bool isPS = (theTrackerGeom->getDetectorType(stackDetid)==TrackerGeometry::ModuleType::Ph2PSP);
	  if (isPS) tmp_trk_nPSstub++;
	}
      }
      */

      // ----------------------------------------------------------------------------------------------
      // loop over stubs on tracks
      /*
      float tmp_trk_bend_chi2 = 0;
      if (SaveStubs) {

	// loop over stubs
	for (int is=0; is<tmp_trk_nstub; is++) {

	  //detID of stub
	  DetId detIdStub = theTrackerGeom->idToDet( (stubRefs.at(is)->getClusterRef(0))->getDetId() )->geographicalId();

	  MeasurementPoint coords = stubRefs.at(is)->getClusterRef(0)->findAverageLocalCoordinatesCentered();
	  const GeomDet* theGeomDet = theTrackerGeom->idToDet(detIdStub);
	  Global3DPoint posStub = theGeomDet->surface().toGlobal( theGeomDet->topology().localPosition(coords) );

	  double x=posStub.x();
	  double y=posStub.y();
	  double z=posStub.z();

	  int isBarrel=-1;
	  int layer=-999999;
	  if ( detIdStub.subdetId()==StripSubdetector::TOB ) {
	    isBarrel=1;
	    layer  = static_cast<int>(tTopo->layer(detIdStub));
	    if (DebugMode) cout << "   stub in layer " << layer << " at position x y z = " << x << " " << y << " " << z << endl;
	  }
	  else if ( detIdStub.subdetId()==StripSubdetector::TID ) {
	    isBarrel=0;
	    layer  = static_cast<int>(tTopo->layer(detIdStub));
	    if (DebugMode) cout << "   stub in disk " << layer << " at position x y z = " << x << " " << y << " " << z << endl;
	  }


	  DetId stackDetid = tTopo->stack(detIdStub);
	  bool isPS = (theTrackerGeom->getDetectorType(stackDetid)==TrackerGeometry::ModuleType::Ph2PSP);

	  float pitch = 0.089;
	  float sigma_bend = 0.45;

	  if (isPS) pitch = 0.099;
	  double tmp_stub_r = posStub.perp();

	  float signedPt = 0.3*3.811202/100.0/(iterL1Track->getRInv());
	  float trackBend = -(1.8*0.57*tmp_stub_r/100)/(pitch*signedPt);

	  float stubBend = stubRefs.at(is)->getTriggerBend();
	  if ( !isBarrel && z<0 ) stubBend=-stubBend;
	  float tmp_bend_diff = stubBend - trackBend;
	  float bend_chi2 = (tmp_bend_diff)*(tmp_bend_diff)/(sigma_bend*sigma_bend);
	  tmp_trk_bend_chi2 += bend_chi2;

	}//end loop over stubs
      }
      */
      // ----------------------------------------------------------------------------------------------


      int tmp_trk_genuine = 0;
      int tmp_trk_loose = 0;
      int tmp_trk_unknown = 0;
      int tmp_trk_combinatoric = 0;
      if (MCTruthTTTrackHandle->isLooselyGenuine(l1track_ptr)) tmp_trk_loose = 1;
      if (MCTruthTTTrackHandle->isGenuine(l1track_ptr)) tmp_trk_genuine = 1;
      if (MCTruthTTTrackHandle->isUnknown(l1track_ptr)) tmp_trk_unknown = 1;
      if (MCTruthTTTrackHandle->isCombinatoric(l1track_ptr)) tmp_trk_combinatoric = 1;

      if (DebugMode) {
	cout << "L1 track, pt: " << tmp_trk_pt << " eta: " << tmp_trk_eta << " phi: " << tmp_trk_phi
	     << " z0: " << tmp_trk_z0 << " chi2: " << tmp_trk_chi2 << " nstub: " << tmp_trk_nstub;
	if (tmp_trk_genuine) cout << " (is genuine)" << endl;
	if (tmp_trk_unknown) cout << " (is unknown)" << endl;
	if (tmp_trk_combinatoric) cout << " (is combinatoric)" << endl;
      }

      m_trk_pt ->push_back(tmp_trk_pt);
      m_trk_eta->push_back(tmp_trk_eta);
      m_trk_phi->push_back(tmp_trk_phi);
      m_trk_charge->push_back(tmp_trk_charge);
      m_trk_z0 ->push_back(tmp_trk_z0);
      if (L1Tk_nPar==5) m_trk_d0->push_back(tmp_trk_d0);
      else m_trk_d0->push_back(999.);
      m_trk_chi2 ->push_back(tmp_trk_chi2);
      m_trk_nstub->push_back(tmp_trk_nstub);
      if (SaveTracklet) m_trk_seed->push_back(tmp_trk_seed);
      m_trk_genuine->push_back(tmp_trk_genuine);
      m_trk_loose->push_back(tmp_trk_loose);
      m_trk_unknown->push_back(tmp_trk_unknown);
      m_trk_combinatoric->push_back(tmp_trk_combinatoric);


      // ----------------------------------------------------------------------------------------------
      // for studying the fake rate
      // ----------------------------------------------------------------------------------------------

      edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);

      int myFake = 0;

      int myTP_pdgid = -999;
      float myTP_pt = -999;
      float myTP_eta = -999;
      float myTP_phi = -999;
      float myTP_z0 = -999;
      float myTP_dxy = -999;

      if (my_tp.isNull()) myFake = 0;
      else {
	int tmp_eventid = my_tp->eventId().event();

	if (tmp_eventid > 0) myFake = 2;
	else myFake = 1;

	myTP_pdgid = my_tp->pdgId();
	myTP_pt = my_tp->p4().pt();
	myTP_eta = my_tp->p4().eta();
	myTP_phi = my_tp->p4().phi();
	myTP_z0 = my_tp->vertex().z();

	float myTP_x0 = my_tp->vertex().x();
	float myTP_y0 = my_tp->vertex().y();
	myTP_dxy = sqrt(myTP_x0*myTP_x0 + myTP_y0*myTP_y0);

	if (DebugMode) {
	  cout << "TP matched to track has pt = " << my_tp->p4().pt() << " eta = " << my_tp->momentum().eta()
	       << " phi = " << my_tp->momentum().phi() << " z0 = " << my_tp->vertex().z()
	       << " pdgid = " <<  my_tp->pdgId() << " dxy = " << myTP_dxy << endl;
	}
      }

      m_trk_fake->push_back(myFake);

      m_trk_matchtp_pdgid->push_back(myTP_pdgid);
      m_trk_matchtp_pt->push_back(myTP_pt);
      m_trk_matchtp_eta->push_back(myTP_eta);
      m_trk_matchtp_phi->push_back(myTP_phi);
      m_trk_matchtp_z0->push_back(myTP_z0);
      m_trk_matchtp_dxy->push_back(myTP_dxy);


      // ----------------------------------------------------------------------------------------------
      // for tracking in jets
      // ----------------------------------------------------------------------------------------------

      if (TrackingInJets) {

	if (DebugMode) cout << "doing tracking in jets now" << endl;

	int InJet = 0;
	int InJetHighpt = 0;
	int InJetVeryHighpt = 0;

	for (int ij=0; ij<(int)v_jets.size(); ij++) {
	  float deta = tmp_trk_eta - (v_jets.at(ij)).eta();
	  float dphi = tmp_trk_phi - (v_jets.at(ij)).phi();
	  while (dphi > 3.14159) dphi = fabs(2*3.14159 - dphi);
	  float dR = sqrt(deta*deta + dphi*dphi);

	  if (dR < 0.4) {
	    InJet = 1;
	    if (v_jets_highpt.at(ij) == 1) InJetHighpt = 1;
	    if (v_jets_vhighpt.at(ij) == 1) InJetVeryHighpt = 1;
	    if (ij<NJETS) jets_trk_sumpt[ij] += tmp_trk_pt;
	  }
	}

	m_trk_injet->push_back(InJet);
	m_trk_injet_highpt->push_back(InJetHighpt);
	m_trk_injet_vhighpt->push_back(InJetVeryHighpt);

      }//end tracking in jets

    }//end track loop

  }//end if SaveAllTracks



  // ----------------------------------------------------------------------------------------------
  // loop over tracking particles
  // ----------------------------------------------------------------------------------------------

  edm::Handle<edm::SimVertexContainer> sim_vertices;
  iEvent.getByToken(simVertexInput_, sim_vertices);
  const edm::SimVertexContainer & sim_vert = *sim_vertices.product();

  if (DebugMode) cout << endl << "Loop over tracking particles!" << endl;

  int this_tp = 0;
  std::vector< TrackingParticle >::const_iterator iterTP;
  for (iterTP = TrackingParticleHandle->begin(); iterTP != TrackingParticleHandle->end(); ++iterTP) {

    edm::Ptr< TrackingParticle > tp_ptr(TrackingParticleHandle, this_tp);
    this_tp++;

    /*// SimTrack selection
    if (t.noVertex()) continue;
    if (t.noGenpart()) continue;
    // only muons
    if (std::abs(t.type()) != 13) continue;
    // pt selection
    if (t.momentum().pt() < 2) continue;
    // eta selection
    const float eta(std::abs(t.momentum().eta()));
    if (eta > 2.5) continue;


    if (1){
      std::cout << "Processing SimTrack " << this_tp << std::endl;
      std::cout << "pt (GeV) = " << t.momentum().pt() << ", eta = " << t.momentum().eta()
                << ", phi = " << t.momentum().phi() << ", Q = " << t.charge() << std::endl;
    }*/

    int tmp_eventid = iterTP->eventId().event();
    if (MyProcess != 1 && tmp_eventid > 0) continue; //only care about tracking particles from the primary interaction (except for MyProcess==1, i.e. looking at all TPs)

    float tmp_tp_pt  = iterTP->pt();
    float tmp_tp_eta = iterTP->eta();
    float tmp_tp_phi = iterTP->phi();
    float tmp_tp_vz  = iterTP->vz();
    float tmp_tp_vx  = iterTP->vx();
    float tmp_tp_vy  = iterTP->vy();
    int tmp_tp_pdgid = iterTP->pdgId();
    float tmp_tp_z0_prod = tmp_tp_vz;
    float tmp_tp_d0_prod = -tmp_tp_vx*sin(tmp_tp_phi) + tmp_tp_vy*cos(tmp_tp_phi);

    if (MyProcess==13 && abs(tmp_tp_pdgid) != 13) continue;
    if (MyProcess==11 && abs(tmp_tp_pdgid) != 11) continue;
    if ((MyProcess==6 || MyProcess==15 || MyProcess==211) && abs(tmp_tp_pdgid) != 211) continue;

    if (tmp_tp_pt < TP_minPt) continue;
    if (fabs(tmp_tp_eta) > TP_maxEta) continue;


    // ----------------------------------------------------------------------------------------------
    // get d0/z0 propagated back to the IP

    float tmp_tp_t = tan(2.0*atan(1.0)-2.0*atan(exp(-tmp_tp_eta)));

    float delx = -tmp_tp_vx;
    float dely = -tmp_tp_vy;

    float A = 0.01*0.5696;
    float Kmagnitude = A / tmp_tp_pt;

    float tmp_tp_charge = tp_ptr->charge();
    float K = Kmagnitude * tmp_tp_charge;
    float d = 0;

    float tmp_tp_x0p = delx - (d + 1./(2. * K)*sin(tmp_tp_phi));
    float tmp_tp_y0p = dely + (d + 1./(2. * K)*cos(tmp_tp_phi));
    float tmp_tp_rp = sqrt(tmp_tp_x0p*tmp_tp_x0p + tmp_tp_y0p*tmp_tp_y0p);
    float tmp_tp_d0 = tmp_tp_charge*tmp_tp_rp - (1. / (2. * K));

    tmp_tp_d0 = tmp_tp_d0*(-1); //fix d0 sign

    static double pi = 4.0*atan(1.0);
    float delphi = tmp_tp_phi-atan2(-K*tmp_tp_x0p,K*tmp_tp_y0p);
    if (delphi<-pi) delphi+=2.0*pi;
    if (delphi>pi) delphi-=2.0*pi;
    float tmp_tp_z0 = tmp_tp_vz+tmp_tp_t*delphi/(2.0*K);
    // ----------------------------------------------------------------------------------------------

    if (fabs(tmp_tp_z0) > TP_maxZ0) continue;


    // for pions in ttbar, only consider TPs coming from near the IP!
    float dxy = sqrt(tmp_tp_vx*tmp_tp_vx + tmp_tp_vy*tmp_tp_vy);
    float tmp_tp_dxy = dxy;
    if (MyProcess==6 && (dxy > 1.0)) continue;

    if (DebugMode) cout << "Tracking particle, pt: " << tmp_tp_pt << " eta: " << tmp_tp_eta << " phi: " << tmp_tp_phi
			<< " z0: " << tmp_tp_z0 << " d0: " << tmp_tp_d0
			<< " z_prod: " << tmp_tp_z0_prod << " d_prod: " << tmp_tp_d0_prod
			<< " pdgid: " << tmp_tp_pdgid << " eventID: " << iterTP->eventId().event()
			<< " ttclusters " << MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr).size()
			<< " ttstubs " << MCTruthTTStubHandle->findTTStubRefs(tp_ptr).size()
			<< " tttracks " << MCTruthTTTrackHandle->findTTTrackPtrs(tp_ptr).size() << endl;


    // ----------------------------------------------------------------------------------------------
    // only consider TPs associated with >= 1 cluster, or >= X stubs, or have stubs in >= X layers (configurable options)

    if (MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr).size() < 1) {
      if (DebugMode) cout << "No matching TTClusters for TP, continuing..." << endl;
      continue;
    }


    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > theStubRefs = MCTruthTTStubHandle->findTTStubRefs(tp_ptr);
    int nStubTP = (int) theStubRefs.size();


    // how many layers/disks have stubs?
    int hasStubInLayer[11] = {0};
    for (unsigned int is=0; is<theStubRefs.size(); is++) {

      DetId detid( theStubRefs.at(is)->getDetId() );

      int layer = -1;
      if ( detid.subdetId()==StripSubdetector::TOB ) {
	layer = static_cast<int>(tTopo->layer(detid)) - 1; //fill in array as entries 0-5
      }
      else if ( detid.subdetId()==StripSubdetector::TID ) {
	layer = static_cast<int>(tTopo->layer(detid)) + 5; //fill in array as entries 6-10
      }

      //bool isPS = (theTrackerGeom->getDetectorType(detid)==TrackerGeometry::ModuleType::Ph2PSP);

      //treat genuine stubs separately (==2 is genuine, ==1 is not)
      if (MCTruthTTStubHandle->findTrackingParticlePtr(theStubRefs.at(is)).isNull() && hasStubInLayer[layer]<2)
	hasStubInLayer[layer] = 1;
      else
	hasStubInLayer[layer] = 2;
    }

    int nStubLayerTP = 0;
    int nStubLayerTP_g = 0;
    for (int isum=0; isum<11; isum++) {
      if ( hasStubInLayer[isum] >= 1) nStubLayerTP   += 1;
      if ( hasStubInLayer[isum] == 2) nStubLayerTP_g += 1;
    }

    if (DebugMode) cout << "TP is associated with " << nStubTP << " stubs, and has stubs in " << nStubLayerTP
			<< " different layers/disks, and has GENUINE stubs in " << nStubLayerTP_g << " layers " << endl;



    if (TP_minNStub > 0) {
      if (DebugMode) cout << "Only consider TPs with >= " << TP_minNStub << " stubs" << endl;
      if (nStubTP < TP_minNStub) {
	if (DebugMode) cout << "TP fails minimum nbr stubs requirement! Continuing..." << endl;
	continue;
      }
    }
    if (TP_minNStubLayer > 0) {
      if (DebugMode) cout << "Only consider TPs with stubs in >= " << TP_minNStubLayer << " layers/disks" << endl;
      if (nStubLayerTP < TP_minNStubLayer) {
	if (DebugMode) cout << "TP fails stubs in minimum nbr of layers/disks requirement! Continuing..." << endl;
	continue;
      }
    }


    // ----------------------------------------------------------------------------------------------
    // look for L1 regional muon candidates matched to the tracking particle

    // get the associated simtrack
    const SimTrack& t(tp_ptr->g4Tracks()[0]);
    if(abs(t.type())==13){
      //run arcane muon simulation truth matching
      SimTrackMatchManager match(t, sim_vert[t.vertIndex()], simTrackCfg, iEvent, iSetup,
                               genParticleInput_,
                               simVertexInput_,
                               simTrackInput_,
                               gemSimHitInput_,
                               cscSimHitInput_,
                               rpcSimHitInput_,
                               me0SimHitInput_,
                               dtSimHitInput_,
                               gemDigiInput_,
                               gemPadDigiInput_,
                               gemCoPadDigiInput_,
                               gemRecHitInput_,
                               me0DigiInput_,
                               me0RecHitInput_,
                               me0SegmentInput_,
                               cscComparatorDigiInput_,
                               cscWireDigiInput_,
                               clctInputs_,
                               alctInputs_,
                               lctInputs_,
                               mplctInputs_,
                               cscRecHit2DInput_,
                               cscSegmentInput_,
                               dtDigiInput_,
                               dtStubInput_,
                               dtRecHit1DPairInput_,
                               dtRecSegment2DInput_,
                               dtRecSegment4DInput_,
                               rpcDigiInput_,
                               rpcRecHitInput_,
                               emtfTrackInputLabel_,
                               regMuonCandInputLabel_,
                               gmtInputLabel_,
                               //trackInputLabel_,
                               //trackMuonInputLabel_,
                               recoTrackExtraInputLabel_,
                               recoTrackInputLabel_,
                               recoChargedCandidateInputLabel_
                               );

      //access best match
      TFCand *muonCandidate = match.l1Muons().bestGMTCand();
      std::cout<<"found muon, try to get match"<<std::endl;
      //if best regional match exists, dump info
      if(muonCandidate){
	std::cout<<"found match, extract information"<<std::endl;
        m_matchmuon_pt->push_back(muonCandidate->pt());
        m_matchmuon_eta->push_back(muonCandidate->eta());
        m_matchmuon_phi->push_back(muonCandidate->phi());
        m_matchmuon_charge->push_back(muonCandidate->charge());
	m_matchmuon_quality->push_back(muonCandidate->quality());
        int type=-1;
        if(muonCandidate->tracktype()==l1t::tftype::bmtf) type=2;
        else if(muonCandidate->tracktype()==l1t::tftype::omtf_pos || muonCandidate->tracktype()==l1t::tftype::omtf_neg) type=1;
        else if(muonCandidate->tracktype()==l1t::tftype::emtf_pos || muonCandidate->tracktype()==l1t::tftype::emtf_neg) type=0;
        m_matchmuon_type->push_back(type);
      }
      else{
        m_matchmuon_pt->push_back(-999.);
        m_matchmuon_eta->push_back(-999.);
        m_matchmuon_phi->push_back(-999.);
        m_matchmuon_charge->push_back(-999);
        m_matchmuon_type->push_back(-999);
	m_matchmuon_quality->push_back(-999);
      }
    }
    else{
      m_matchmuon_pt->push_back(-999.);
      m_matchmuon_eta->push_back(-999.);
      m_matchmuon_phi->push_back(-999.);
      m_matchmuon_charge->push_back(-999);
      m_matchmuon_type->push_back(-999);
      m_matchmuon_quality->push_back(-999);
    }
    


    // ----------------------------------------------------------------------------------------------
    // look for L1 tracks matched to the tracking particle

    std::vector< edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > > matchedTracks = MCTruthTTTrackHandle->findTTTrackPtrs(tp_ptr);

    int nMatch = 0;
    int nLooseMatch = 0;
    int i_track = -1;
    int i_loosetrack = -1;
    float i_chi2dof = 99999;
    float i_loosechi2dof = 99999;

    if (matchedTracks.size() > 0) {

      if (DebugMode && (matchedTracks.size()>1)) cout << "TrackingParticle has more than one matched L1 track!" << endl;

      // ----------------------------------------------------------------------------------------------
      // loop over matched L1 tracks
      // here, "match" means tracks that can be associated to a TrackingParticle with at least one hit of at least one of its clusters
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SLHCTrackerTriggerSWTools#MC_truth_for_TTTrack

      for (int it=0; it<(int)matchedTracks.size(); it++) {

	bool tmp_trk_genuine = false;
	bool tmp_trk_loosegenuine = false;
	if (MCTruthTTTrackHandle->isGenuine(matchedTracks.at(it))) tmp_trk_genuine = true;
	if (MCTruthTTTrackHandle->isLooselyGenuine(matchedTracks.at(it))) tmp_trk_loosegenuine = true;
	if (!tmp_trk_loosegenuine) continue;


	if (DebugMode) {
	  if (MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it)).isNull()) {
	    cout << "track matched to TP is NOT uniquely matched to a TP" << endl;
	  }
	  else {
	    edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it));
	    cout << "TP matched to track matched to TP ... tp pt = " << my_tp->p4().pt() << " eta = " << my_tp->momentum().eta()
		 << " phi = " << my_tp->momentum().phi() << " z0 = " << my_tp->vertex().z() << endl;
	  }
	  cout << "   ... matched L1 track has pt = " << matchedTracks.at(it)->getMomentum(L1Tk_nPar).perp()
	       << " eta = " << matchedTracks.at(it)->getMomentum(L1Tk_nPar).eta()
	       << " phi = " << matchedTracks.at(it)->getMomentum(L1Tk_nPar).phi()
	       << " chi2 = " << matchedTracks.at(it)->getChi2(L1Tk_nPar)
	       << " consistency = " << matchedTracks.at(it)->getStubPtConsistency(L1Tk_nPar)
	       << " z0 = " << matchedTracks.at(it)->getPOCA(L1Tk_nPar).z()
	       << " nstub = " << matchedTracks.at(it)->getStubRefs().size();
	  if (tmp_trk_genuine) cout << " (genuine!) " << endl;
	  if (tmp_trk_loosegenuine) cout << " (loose genuine!) " << endl;
	}


	// ----------------------------------------------------------------------------------------------
	// further require L1 track to be (loosely) genuine, that there is only one TP matched to the track
	// + have >= L1Tk_minNStub stubs for it to be a valid match (only relevant is your track collection
	// e.g. stores 3-stub tracks but at plot level you require >= 4 stubs (--> tracklet case)

	std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubRefs = matchedTracks.at(it)->getStubRefs();
	int tmp_trk_nstub = stubRefs.size();

	if (tmp_trk_nstub < L1Tk_minNStub) continue;

	/*
	// PS stubs
	int tmp_trk_nPSstub = 0;
	if (SaveTracklet) {
	  for (int is=0; is<tmp_trk_nstub; is++) {
	    DetId detIdStub = theTrackerGeom->idToDet( (stubRefs.at(is)->getClusterRef(0))->getDetId() )->geographicalId();
	    DetId stackDetid = tTopo->stack(detIdStub);
	    bool isPS = (theTrackerGeom->getDetectorType(stackDetid)==TrackerGeometry::ModuleType::Ph2PSP);
	    if (isPS) tmp_trk_nPSstub++;
	  }
	}
	*/

	float dmatch_pt  = 999;
	float dmatch_eta = 999;
	float dmatch_phi = 999;
	int match_id = 999;

	edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it));
	dmatch_pt  = fabs(my_tp->p4().pt() - tmp_tp_pt);
	dmatch_eta = fabs(my_tp->p4().eta() - tmp_tp_eta);
	dmatch_phi = fabs(my_tp->p4().phi() - tmp_tp_phi);
	match_id = my_tp->pdgId();

	float tmp_trk_chi2dof = (matchedTracks.at(it)->getChi2(L1Tk_nPar)) / (2*tmp_trk_nstub - L1Tk_nPar);

	// ensure that track is uniquely matched to the TP we are looking at!
	if (dmatch_pt<0.1 && dmatch_eta<0.1 && dmatch_phi<0.1 && tmp_tp_pdgid==match_id) {
	  nLooseMatch++;
	  if (i_loosetrack < 0 || tmp_trk_chi2dof < i_loosechi2dof) {
	    i_loosetrack = it;
	    i_loosechi2dof = tmp_trk_chi2dof;
	  }

	  if (tmp_trk_genuine) {
	    nMatch++;
	    if (i_track < 0 || tmp_trk_chi2dof < i_chi2dof) {
	      i_track = it;
	      i_chi2dof = tmp_trk_chi2dof;
	    }
	  }
	}


      }// end loop over matched L1 tracks

    }// end has at least 1 matched L1 track
    // ----------------------------------------------------------------------------------------------


    float tmp_matchtrk_pt   = -999;
    float tmp_matchtrk_eta  = -999;
    float tmp_matchtrk_phi  = -999;
    float tmp_matchtrk_z0   = -999;
    float tmp_matchtrk_d0   = -999;
    float tmp_matchtrk_chi2 = -999;
    int tmp_matchtrk_nstub  = -999;
    int tmp_matchtrk_seed   = -999;
    bool tmp_matchtrk_charge= 0;

    float tmp_loosematchtrk_pt   = -999;
    float tmp_loosematchtrk_eta  = -999;
    float tmp_loosematchtrk_phi  = -999;
    float tmp_loosematchtrk_z0   = -999;
    float tmp_loosematchtrk_d0   = -999;
    float tmp_loosematchtrk_chi2 = -999;
    int tmp_loosematchtrk_nstub  = -999;
    int tmp_loosematchtrk_seed   = -999;


    if (nMatch > 1 && DebugMode) cout << "WARNING *** 2 or more matches to genuine L1 tracks ***" << endl;
    if (nLooseMatch > 1 && DebugMode) cout << "WARNING *** 2 or more matches to loosely genuine L1 tracks ***" << endl;

    if (nMatch > 0) {
      tmp_matchtrk_pt   = matchedTracks.at(i_track)->getMomentum(L1Tk_nPar).perp();
      tmp_matchtrk_eta  = matchedTracks.at(i_track)->getMomentum(L1Tk_nPar).eta();
      tmp_matchtrk_phi  = matchedTracks.at(i_track)->getMomentum(L1Tk_nPar).phi();
      tmp_matchtrk_z0   = matchedTracks.at(i_track)->getPOCA(L1Tk_nPar).z();
      tmp_matchtrk_charge=signbit(matchedTracks.at(i_track)->getRInv(L1Tk_nPar));

      if (L1Tk_nPar == 5) {
	float tmp_matchtrk_x0 = matchedTracks.at(i_track)->getPOCA(L1Tk_nPar).x();
	float tmp_matchtrk_y0 = matchedTracks.at(i_track)->getPOCA(L1Tk_nPar).y();
	tmp_matchtrk_d0 = -tmp_matchtrk_x0*sin(tmp_matchtrk_phi) + tmp_matchtrk_y0*cos(tmp_matchtrk_phi);
      }

      tmp_matchtrk_chi2 = matchedTracks.at(i_track)->getChi2(L1Tk_nPar);
      tmp_matchtrk_nstub  = (int) matchedTracks.at(i_track)->getStubRefs().size();
      if (SaveTracklet)	tmp_matchtrk_seed = (int) matchedTracks.at(i_track)->getWedge();




      // ------------------------------------------------------------------------------------------
      /*
      tmp_matchtrk_bend_chi2 = 0;

      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubRefs = matchedTracks.at(i_track)->getStubRefs();
      int tmp_nstub = stubRefs.size();


      for (int is=0; is<tmp_nstub; is++) {

	DetId detIdStub = theTrackerGeom->idToDet( (stubRefs.at(is)->getClusterRef(0))->getDetId() )->geographicalId();

	MeasurementPoint coords = stubRefs.at(is)->getClusterRef(0)->findAverageLocalCoordinatesCentered();
	const GeomDet* theGeomDet = theTrackerGeom->idToDet(detIdStub);
	Global3DPoint posStub = theGeomDet->surface().toGlobal( theGeomDet->topology().localPosition(coords) );

	double x=posStub.x();
	double y=posStub.y();
	double z=posStub.z();

	int isBarrel = 0;
	int layer=-999999;
	if ( detIdStub.subdetId()==StripSubdetector::TOB ) {
	  isBarrel = 1;
	  layer  = static_cast<int>(tTopo->layer(detIdStub));
	}
	else if ( detIdStub.subdetId()==StripSubdetector::TID ) {
	  isBarrel = 0;
	  layer  = static_cast<int>(tTopo->layer(detIdStub));
	}

	DetId stackDetid = tTopo->stack(detIdStub);
	bool isPS = (theTrackerGeom->getDetectorType(stackDetid)==TrackerGeometry::ModuleType::Ph2PSP);

	float pitch = 0.089;
	float sigma_bend = 0.45;
	if (isPS) pitch = 0.099;
	double tmp_stub_r = posStub.perp();

	float signedPt = 0.3*3.811202/100.0/(matchedTracks.at(i_track)->getRInv());
	float trackBend = -(1.8*0.57*tmp_stub_r/100)/(pitch*signedPt);

	float stubBend = stubRefs.at(is)->getTriggerBend();
	if ( !isBarrel && z>0 ) stubBend=-stubBend;
	float tmp_bend_diff = stubBend - trackBend;
	float bend_chi2 = (tmp_bend_diff)*(tmp_bend_diff)/(sigma_bend*sigma_bend);
	tmp_matchtrk_bend_chi2 += bend_chi2;
      }
      */
      // ------------------------------------------------------------------------------------------


    }

    if (nLooseMatch > 0) {
      tmp_loosematchtrk_pt   = matchedTracks.at(i_loosetrack)->getMomentum(L1Tk_nPar).perp();
      tmp_loosematchtrk_eta  = matchedTracks.at(i_loosetrack)->getMomentum(L1Tk_nPar).eta();
      tmp_loosematchtrk_phi  = matchedTracks.at(i_loosetrack)->getMomentum(L1Tk_nPar).phi();
      tmp_loosematchtrk_z0   = matchedTracks.at(i_loosetrack)->getPOCA(L1Tk_nPar).z();

      if (L1Tk_nPar == 5) {
	float tmp_loosematchtrk_x0 = matchedTracks.at(i_loosetrack)->getPOCA(L1Tk_nPar).x();
	float tmp_loosematchtrk_y0 = matchedTracks.at(i_loosetrack)->getPOCA(L1Tk_nPar).y();
	tmp_loosematchtrk_d0 = -tmp_loosematchtrk_x0*sin(tmp_loosematchtrk_phi) + tmp_loosematchtrk_y0*cos(tmp_loosematchtrk_phi);
      }

      tmp_loosematchtrk_chi2 = matchedTracks.at(i_loosetrack)->getChi2(L1Tk_nPar);
      tmp_loosematchtrk_nstub  = (int) matchedTracks.at(i_loosetrack)->getStubRefs().size();
      if (SaveTracklet)	tmp_loosematchtrk_seed = (int) matchedTracks.at(i_loosetrack)->getWedge();
    }


    m_tp_pt->push_back(tmp_tp_pt);
    m_tp_eta->push_back(tmp_tp_eta);
    m_tp_phi->push_back(tmp_tp_phi);
    m_tp_dxy->push_back(tmp_tp_dxy);
    m_tp_z0->push_back(tmp_tp_z0);
    m_tp_d0->push_back(tmp_tp_d0);
    m_tp_z0_prod->push_back(tmp_tp_z0_prod);
    m_tp_d0_prod->push_back(tmp_tp_d0_prod);
    m_tp_pdgid->push_back(tmp_tp_pdgid);
    m_tp_nmatch->push_back(nMatch);
    m_tp_nloosematch->push_back(nLooseMatch);
    m_tp_nstub->push_back(nStubTP);
    m_tp_eventid->push_back(tmp_eventid);
    m_tp_charge->push_back(tmp_tp_charge);

    m_matchtrk_pt ->push_back(tmp_matchtrk_pt);
    m_matchtrk_eta->push_back(tmp_matchtrk_eta);
    m_matchtrk_phi->push_back(tmp_matchtrk_phi);
    m_matchtrk_z0 ->push_back(tmp_matchtrk_z0);
    m_matchtrk_d0 ->push_back(tmp_matchtrk_d0);
    m_matchtrk_chi2 ->push_back(tmp_matchtrk_chi2);
    m_matchtrk_nstub->push_back(tmp_matchtrk_nstub);
    if (SaveTracklet) m_matchtrk_seed->push_back(tmp_matchtrk_seed);
    m_matchtrk_charge->push_back(tmp_matchtrk_charge);

    m_loosematchtrk_pt ->push_back(tmp_loosematchtrk_pt);
    m_loosematchtrk_eta->push_back(tmp_loosematchtrk_eta);
    m_loosematchtrk_phi->push_back(tmp_loosematchtrk_phi);
    m_loosematchtrk_z0 ->push_back(tmp_loosematchtrk_z0);
    m_loosematchtrk_d0 ->push_back(tmp_loosematchtrk_d0);
    m_loosematchtrk_chi2 ->push_back(tmp_loosematchtrk_chi2);
    m_loosematchtrk_nstub->push_back(tmp_loosematchtrk_nstub);
    if (SaveTracklet) m_loosematchtrk_seed->push_back(tmp_loosematchtrk_seed);


    // ----------------------------------------------------------------------------------------------
    // for tracking in jets
    // ----------------------------------------------------------------------------------------------

    if (TrackingInJets) {

      if (DebugMode) cout << "check if TP/matched track is within jet" << endl;

      int tp_InJet = 0;
      int matchtrk_InJet = 0;
      int loosematchtrk_InJet = 0;
      int tp_InJetHighpt = 0;
      int matchtrk_InJetHighpt = 0;
      int loosematchtrk_InJetHighpt = 0;
      int tp_InJetVeryHighpt = 0;
      int matchtrk_InJetVeryHighpt = 0;
      int loosematchtrk_InJetVeryHighpt = 0;

      for (int ij=0; ij<(int)v_jets.size(); ij++) {
	float deta = tmp_tp_eta - (v_jets.at(ij)).eta();
	float dphi = tmp_tp_phi - (v_jets.at(ij)).phi();
	while (dphi > 3.14159) dphi = fabs(2*3.14159 - dphi);
	float dR = sqrt(deta*deta + dphi*dphi);
	if (dR < 0.4) {
	  tp_InJet = 1;
	  if (v_jets_highpt.at(ij) == 1) tp_InJetHighpt = 1;
	  if (v_jets_vhighpt.at(ij) == 1) tp_InJetVeryHighpt = 1;
	  if (ij<NJETS) jets_tp_sumpt[ij] += tmp_tp_pt;
	}

	if (nMatch > 0) {
	  deta = tmp_matchtrk_eta - (v_jets.at(ij)).eta();
	  dphi = tmp_matchtrk_phi - (v_jets.at(ij)).phi();
	  while (dphi > 3.14159) dphi = fabs(2*3.14159 - dphi);
	  dR = sqrt(deta*deta + dphi*dphi);
	  if (dR < 0.4) {
	    matchtrk_InJet = 1;
	    if (v_jets_highpt.at(ij) == 1) matchtrk_InJetHighpt = 1;
	    if (v_jets_vhighpt.at(ij) == 1) matchtrk_InJetVeryHighpt = 1;
	    if (ij<NJETS) jets_matchtrk_sumpt[ij] += tmp_matchtrk_pt;
	  }
	}

	if (nLooseMatch > 0) {
	  deta = tmp_loosematchtrk_eta - (v_jets.at(ij)).eta();
	  dphi = tmp_loosematchtrk_phi - (v_jets.at(ij)).phi();
	  while (dphi > 3.14159) dphi = fabs(2*3.14159 - dphi);
	  dR = sqrt(deta*deta + dphi*dphi);
	  if (dR < 0.4) {
	    loosematchtrk_InJet = 1;
	    if (v_jets_highpt.at(ij) == 1) loosematchtrk_InJetHighpt = 1;
	    if (v_jets_vhighpt.at(ij) == 1) loosematchtrk_InJetVeryHighpt = 1;
	    if (ij<NJETS) jets_loosematchtrk_sumpt[ij] += tmp_loosematchtrk_pt;
	  }
	}
      }

      m_tp_injet->push_back(tp_InJet);
      m_tp_injet_highpt->push_back(tp_InJetHighpt);
      m_tp_injet_vhighpt->push_back(tp_InJetVeryHighpt);
      m_matchtrk_injet->push_back(matchtrk_InJet);
      m_matchtrk_injet_highpt->push_back(matchtrk_InJetHighpt);
      m_matchtrk_injet_vhighpt->push_back(matchtrk_InJetVeryHighpt);
      m_loosematchtrk_injet->push_back(loosematchtrk_InJet);
      m_loosematchtrk_injet_highpt->push_back(loosematchtrk_InJetHighpt);
      m_loosematchtrk_injet_vhighpt->push_back(loosematchtrk_InJetVeryHighpt);

    }//end TrackingInJets

  } //end loop tracking particles


  if (TrackingInJets) {
    for (int ij=0; ij<(int)v_jets.size(); ij++) {
      if (ij<NJETS) {
	m_jet_eta->push_back((v_jets.at(ij)).eta());
	m_jet_phi->push_back((v_jets.at(ij)).phi());
	m_jet_pt->push_back((v_jets.at(ij)).pt());
	m_jet_tp_sumpt->push_back(jets_tp_sumpt[ij]);
	m_jet_trk_sumpt->push_back(jets_trk_sumpt[ij]);
	m_jet_matchtrk_sumpt->push_back(jets_matchtrk_sumpt[ij]);
	m_jet_loosematchtrk_sumpt->push_back(jets_loosematchtrk_sumpt[ij]);
      }
    }
  }

  //finally svens muon stuff
  if (SaveMuons){
    Handle< BXVector<l1t::RegionalMuonCand> > emtfs;
    Handle< BXVector<l1t::RegionalMuonCand> > omtfs;
    Handle< BXVector<l1t::RegionalMuonCand> > bmtfs;

    iEvent.getByToken(m_emtfToken,emtfs);
    iEvent.getByToken(m_omtfToken,omtfs);
    iEvent.getByToken(m_bmtfToken,bmtfs);

    int nEMTF=0;
    for (auto it = emtfs->begin(0); it != emtfs->end(0); it++){
      m_EMTF_muon_eta->push_back(it->hwEta()*0.010875);
      int globPhi=l1t::MicroGMTConfiguration::calcGlobalPhi(it->hwPhi(), it->trackFinderType(), it->processor());
      m_EMTF_muon_phi->push_back(globPhi*2*M_PI/576.);
      m_EMTF_muon_pt->push_back(it->hwPt()*0.5);
      if(!it->hwSignValid()) m_EMTF_muon_c->push_back(0);
      else{
	if(it->hwSign())   m_EMTF_muon_c->push_back(-1);
	else 		   m_EMTF_muon_c->push_back(1);
      }
      ++nEMTF;
    }
    m_EMTF_muon_n->push_back(nEMTF);

    int nOMTF=0;
    for (auto it = omtfs->begin(0); it != omtfs->end(0); it++){
      m_OMTF_muon_eta->push_back(it->hwEta()*0.010875);
      m_OMTF_muon_phi->push_back(l1t::MicroGMTConfiguration::calcGlobalPhi(it->hwPhi(), it->trackFinderType(), it->processor())*2*M_PI/576.);
      m_OMTF_muon_pt->push_back(it->hwPt()*0.5);
      if(!it->hwSignValid()) m_OMTF_muon_c->push_back(0);
      else{
	if(it->hwSign())   m_OMTF_muon_c->push_back(-1);
	else 		   m_OMTF_muon_c->push_back(1);
      }
      ++nOMTF;
    }
    m_OMTF_muon_n->push_back(nOMTF);

    int nBMTF=0;
    for (auto it = bmtfs->begin(0); it != bmtfs->end(0); it++){
      m_BMTF_muon_eta->push_back(it->hwEta()*0.010875);
      m_BMTF_muon_phi->push_back(l1t::MicroGMTConfiguration::calcGlobalPhi(it->hwPhi(), it->trackFinderType(), it->processor())*2*M_PI/576.);
      m_BMTF_muon_pt->push_back(it->hwPt()*0.5);
      if(!it->hwSignValid()) m_BMTF_muon_c->push_back(0);
      else{
	if(it->hwSign())   m_BMTF_muon_c->push_back(-1);
	else 		   m_BMTF_muon_c->push_back(1);
      }
      ++nBMTF;
    }
    m_BMTF_muon_n->push_back(nBMTF);
  }


  eventTree->Fill();


} // end of analyze()


///////////////////////////
// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(L1TrackNtupleMaker);
