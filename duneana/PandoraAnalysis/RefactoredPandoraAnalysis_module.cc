////////////////////////////////////////////////////////////////////////
// Class:       refactoredPandoraAna
// Plugin Type: analyzer (art v3_05_01)
// File:        pandoraAnalysis_module.cc
//
// Generated at Fri Aug  7 15:01:35 2020 by Ryan Cross using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <TTree.h>
#include <vector>
#include <string>
#include <TH1D.h>
#include <TLorentzVector.h>

namespace test {
  class refactoredPandoraAna;
}


class test::refactoredPandoraAna : public art::EDAnalyzer {
public:
  explicit refactoredPandoraAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  refactoredPandoraAna(refactoredPandoraAna const&) = delete;
  refactoredPandoraAna(refactoredPandoraAna&&) = delete;
  refactoredPandoraAna& operator=(refactoredPandoraAna const&) = delete;
  refactoredPandoraAna& operator=(refactoredPandoraAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // High-level functions
  float getPurity(std::vector<art::Ptr<recob::Hit>> const &hits, int const &trackID, int view = -1);
  float getCompleteness(std::vector<art::Ptr<recob::Hit>> const &hits, int const &trackID, int view = -1);

  // Functions per concern
  bool processMCTruth(art::Event const& e);
  bool processMCParticles(art::Event const& e);
  bool parsePFParticles(art::Event const& e);

  void parseTrack(art::Event const& e, const art::Ptr<recob::PFParticle> &pfp, int iPfp);
  void parseShower(art::Event const& e, const art::Ptr<recob::PFParticle> &pfp, int iPfp);
  void recoMCInfo(art::Event const& e, const art::Ptr<recob::PFParticle> &pfp, int iPfp);

  void reset(bool deepClean = false);

private:

  TTree *fTree;

  unsigned int fEventID;
  unsigned int fRunID;
  unsigned int fSubRunID;

  detinfo::DetectorClocksData fClockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

  unsigned int fNMCParticles;
  unsigned int fNPFParticles;

  static const int kNMaxMCParticles = 20000;
  static const int kNMaxPFParticles = 2000;
  static const int kNMaxPFPClusters = 100;
  static const int kNViews = 3;

  int fNuPdgCodeTruth;
  int fNuCCNCTruth;
  int fNuModeTruth;
  double fNuETruth;
  double fNuVertexXTruth;
  double fNuVertexYTruth;
  double fNuVertexZTruth;

  bool        fMCIsPrimary[kNMaxMCParticles];
  int         fMCParticlePdgCode[kNMaxMCParticles];
  double      fMCParticleTrueEnergy [kNMaxMCParticles];
  int         fMCParticleTrackID[kNMaxMCParticles];
  int         fMCParticleParentTrackID[kNMaxMCParticles];
  std::string fMCParticleStartProcess[kNMaxMCParticles];
  std::string fMCParticleEndProcess[kNMaxMCParticles];
  int         fMCParticleNTrajectoryPoint[kNMaxMCParticles];
  double      fMCParticleStartPositionX[kNMaxMCParticles];
  double      fMCParticleStartPositionY[kNMaxMCParticles];
  double      fMCParticleStartPositionZ[kNMaxMCParticles];
  double      fMCParticleStartPositionT[kNMaxMCParticles];
  double      fMCParticleStartMomentumX[kNMaxMCParticles];
  double      fMCParticleStartMomentumY[kNMaxMCParticles];
  double      fMCParticleStartMomentumZ[kNMaxMCParticles];
  double      fMCParticleStartMomentumE[kNMaxMCParticles];
  double      fMCParticleEndPositionX[kNMaxMCParticles];
  double      fMCParticleEndPositionY[kNMaxMCParticles];
  double      fMCParticleEndPositionZ[kNMaxMCParticles];
  double      fMCParticleEndPositionT[kNMaxMCParticles];
  double      fMCParticleEndMomentumX[kNMaxMCParticles];
  double      fMCParticleEndMomentumY[kNMaxMCParticles];
  double      fMCParticleEndMomentumZ[kNMaxMCParticles];
  double      fMCParticleEndMomentumE[kNMaxMCParticles];
  double      fMCParticleVertexTime[kNMaxMCParticles];
  double      fMCParticleEndTime[kNMaxMCParticles];
  int         fMCParticleNHits[kNMaxMCParticles];
  int         fMCParticleNHitsView[kNMaxMCParticles][kNViews];

  int    fPFPID[kNMaxPFParticles];
  bool   fPFPIsPrimary[kNMaxPFParticles];
  int    fPFPTrueParticleMatchedID[kNMaxPFParticles];
  int    fPFPTrueParticleMatchedPosition[kNMaxPFParticles];
  int    fPFPParentID[kNMaxPFParticles];
  int    fPFPPdgCode[kNMaxPFParticles];
  int    fPFPNChildren[kNMaxPFParticles];
  int    fPFPNClusters[kNMaxPFParticles];
  int    fPFPNHits[kNMaxPFParticles];
  int    fPFPNHitsView[kNMaxPFParticles][kNViews];
  int    fPFPNSharedTrueParticleHits[kNMaxPFParticles];
  int    fPFPNSharedTrueParticleHitsView[kNMaxPFParticles][kNViews];
  int    fPFPTrueParticleMatchedIDView[kNMaxPFParticles][kNViews];
  int    fPFPTrueParticleMatchedPositionView[kNMaxPFParticles][kNViews];
  bool   fPFPIsTrack[kNMaxPFParticles];
  bool   fPFPIsShower[kNMaxPFParticles];
  int    fPFPCluPlane[kNMaxPFParticles][kNMaxPFPClusters];
  int    fPFPCluView[kNMaxPFParticles][kNMaxPFPClusters];
  int    fPFPCluNHits[kNMaxPFParticles][kNMaxPFPClusters];
  double fPFPCluIntegral[kNMaxPFParticles][kNMaxPFPClusters];

  int    fPFPTrackID[kNMaxPFParticles];
  double fPFPTrackLength[kNMaxPFParticles];
  double fPFPTrackStartX[kNMaxPFParticles];
  double fPFPTrackStartY[kNMaxPFParticles];
  double fPFPTrackStartZ[kNMaxPFParticles];
  double fPFPTrackVertexX[kNMaxPFParticles];
  double fPFPTrackVertexY[kNMaxPFParticles];
  double fPFPTrackVertexZ[kNMaxPFParticles];
  double fPFPTrackEndX[kNMaxPFParticles];
  double fPFPTrackEndY[kNMaxPFParticles];
  double fPFPTrackEndZ[kNMaxPFParticles];
  double fPFPTrackTheta[kNMaxPFParticles];
  double fPFPTrackPhi[kNMaxPFParticles];
  double fPFPTrackZenithAngle[kNMaxPFParticles];
  double fPFPTrackAzimuthAngle[kNMaxPFParticles];
  double fPFPTrackStartDirectionX[kNMaxPFParticles];
  double fPFPTrackStartDirectionY[kNMaxPFParticles];
  double fPFPTrackStartDirectionZ[kNMaxPFParticles];
  double fPFPTrackVertexDirectionX[kNMaxPFParticles];
  double fPFPTrackVertexDirectionY[kNMaxPFParticles];
  double fPFPTrackVertexDirectionZ[kNMaxPFParticles];
  double fPFPTrackEndDirectionX[kNMaxPFParticles];
  double fPFPTrackEndDirectionY[kNMaxPFParticles];
  double fPFPTrackEndDirectionZ[kNMaxPFParticles];
  float  fPFPTrackChi2[kNMaxPFParticles];
  int    fPFPTrackNdof[kNMaxPFParticles];

  int    fPFPShowerID[kNMaxPFParticles];
  int    fPFPShowerBestPlane[kNMaxPFParticles];
  double fPFPShowerDirectionX[kNMaxPFParticles];
  double fPFPShowerDirectionY[kNMaxPFParticles];
  double fPFPShowerDirectionZ[kNMaxPFParticles];
  double fPFPShowerDirectionErrX[kNMaxPFParticles];
  double fPFPShowerDirectionErrY[kNMaxPFParticles];
  double fPFPShowerDirectionErrZ[kNMaxPFParticles];
  double fPFPShowerStartX[kNMaxPFParticles];
  double fPFPShowerStartY[kNMaxPFParticles];
  double fPFPShowerStartZ[kNMaxPFParticles];
  double fPFPShowerStartErrX[kNMaxPFParticles];
  double fPFPShowerStartErrY[kNMaxPFParticles];
  double fPFPShowerStartErrZ[kNMaxPFParticles];
  double fPFPShowerLength[kNMaxPFParticles];
  double fPFPShowerOpenAngle[kNMaxPFParticles];

  double fPFPMatchedMCPDG[kNMaxMCParticles];
  double fPFPMatchedMCEnergy[kNMaxMCParticles];
  double fPFPMatchedMCNHits[kNMaxMCParticles];

  double fPFPCompleteness[kNMaxMCParticles];
  double fPFPCompletenessView[kNMaxMCParticles][kNViews];
  double fPFPPurity[kNMaxMCParticles];
  double fPFPPurityView[kNMaxMCParticles][kNViews];

  std::string fTruthLabel;
  std::string fSimLabel;
  std::string fHitLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fPFParticleLabel;
  bool fRollUpUnsavedIDs;

  const geo::Geometry* fGeom;
};


test::refactoredPandoraAna::refactoredPandoraAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  fTruthLabel       = p.get<std::string>("TruthLabel");
  fSimLabel         = p.get<std::string>("SimulationLabel");
  fHitLabel         = p.get<std::string>("HitLabel");
  fPFParticleLabel  = p.get<std::string>("PFParticleLabel");
  fTrackLabel       = p.get<std::string>("TrackLabel");
  fShowerLabel      = p.get<std::string>("ShowerLabel");
  fGeom             = &*art::ServiceHandle<geo::Geometry>();
  fRollUpUnsavedIDs = p.get<bool>("RollUpUnsavedIDs");
}


float test::refactoredPandoraAna::getCompleteness(std::vector<art::Ptr<recob::Hit>> const &hits, int const &trackID, int view) {

  std::map<int,int> objectHitsMap;

  for(unsigned int i = 0; i < hits.size(); ++i) {
    int primaryID = TruthMatchUtils::TrueParticleID(fClockData, hits[i], true);
    objectHitsMap[primaryID]++;
  }

  int baseNumOfHits = fMCParticleNHits[fPFPTrueParticleMatchedPosition[trackID]];

  if (view != -1)
    baseNumOfHits = fMCParticleNHitsView[fPFPTrueParticleMatchedPosition[trackID]][0];

  return (baseNumOfHits == 999999) ? 999999 : objectHitsMap[trackID] / static_cast<float>(baseNumOfHits);
}

float test::refactoredPandoraAna::getPurity(std::vector<art::Ptr<recob::Hit>> const &hits, int const &trackID, int view) {

  std::map<int,int> objectHitsMap;

  for(unsigned int i = 0; i < hits.size(); ++i) {
    int primaryID = TruthMatchUtils::TrueParticleID(fClockData, hits[i], true);
    objectHitsMap[primaryID]++;
  }

  return (hits.size() == 0) ? 999999 : objectHitsMap[trackID] / static_cast<float>(hits.size());
}

bool test::refactoredPandoraAna::processMCTruth(art::Event const& e) {
  // Get the MC truth
  const std::vector<art::Ptr<simb::MCTruth>> mcTruthVect = dune_ana::DUNEAnaEventUtils::GetMCTruths(e, fTruthLabel);

  // Find the neutrino and store its information.
  for (unsigned int iTruth = 0; iTruth < mcTruthVect.size(); iTruth++){

    const art::Ptr<simb::MCTruth> truth = mcTruthVect[iTruth];

    if (truth->Origin() != simb::kBeamNeutrino)
      continue;

    fNuPdgCodeTruth = truth->GetNeutrino().Nu().PdgCode();
    fNuCCNCTruth = truth->GetNeutrino().CCNC();
    fNuModeTruth = truth->GetNeutrino().Mode();
    fNuETruth = truth->GetNeutrino().Nu().E();
    fNuVertexXTruth = truth->GetNeutrino().Nu().Vx();
    fNuVertexYTruth = truth->GetNeutrino().Nu().Vy();
    fNuVertexZTruth = truth->GetNeutrino().Nu().Vz();

    break; // Only one beam nu, so don't keep looping
  }

  return true;
}

bool test::refactoredPandoraAna::processMCParticles(art::Event const& e) {
  // Get the MC Particles
  art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fSimLabel);
  if(! mcParticles.isValid())
    return false;

  fNMCParticles = mcParticles->size();

  std::cout << "There is " << fNMCParticles << " MC particles..." << std::endl;

  // Get all hits.
  art::Handle<std::vector<recob::Hit>> hitHandle;
  std::vector<art::Ptr<recob::Hit> > allHits;
  if(e.getByLabel(fHitLabel, hitHandle))
    art::fill_ptr_vector(allHits, hitHandle);

  // Fill MC particle to hits map.
  std::map<int,int> trueParticleHits, trueParticleHitsView0, trueParticleHitsView1, trueParticleHitsView2;
  for (const auto& hit: allHits) {
      TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleID(fClockData, hit, fRollUpUnsavedIDs));

      if (! TruthMatchUtils::Valid(g4ID))
        continue;

      trueParticleHits[g4ID]++;

      if(hit->View()==0)
        trueParticleHitsView0[g4ID]++;
      if(hit->View()==1)
        trueParticleHitsView1[g4ID]++;
      if(hit->View()==2)
        trueParticleHitsView2[g4ID]++;
  }

  if (fNMCParticles > kNMaxMCParticles) {
    std::cout << "WARN: Too many MC Particles: (" << fNMCParticles << ")!" << std::endl;
    std::cout << "WARN: Limiting to first " << kNMaxMCParticles << " particles!" << std::endl;
    fNMCParticles = kNMaxMCParticles;
  }

  for (unsigned int iMc = 0; iMc < mcParticles->size(); iMc++) {

    const simb::MCParticle trueParticle = mcParticles->at(iMc);

    bool isMCPrimary(false);
    trueParticle.Process() == "primary" ? isMCPrimary=true : isMCPrimary=false;
    fMCIsPrimary[iMc] = isMCPrimary;

    fMCParticleTrueEnergy[iMc]       = trueParticle.E();
    fMCParticlePdgCode[iMc]          = trueParticle.PdgCode();
    fMCParticleTrackID[iMc]          = trueParticle.TrackId();
    fMCParticleVertexTime[iMc]       = trueParticle.T();
    fMCParticleEndTime[iMc]          = trueParticle.EndT();
    fMCParticleParentTrackID[iMc]    = trueParticle.Mother();
    fMCParticleStartProcess[iMc]     = trueParticle.Process();
    fMCParticleEndProcess[iMc]       = trueParticle.EndProcess();
    fMCParticleNTrajectoryPoint[iMc] = trueParticle.NumberTrajectoryPoints();
    fMCParticleNHits[iMc]            = trueParticleHits[trueParticle.TrackId()];
    fMCParticleNHitsView[iMc][0]     = trueParticleHitsView0[trueParticle.TrackId()];
    fMCParticleNHitsView[iMc][1]     = trueParticleHitsView1[trueParticle.TrackId()];
    fMCParticleNHitsView[iMc][2]     = trueParticleHitsView2[trueParticle.TrackId()];
    fMCParticleStartPositionX[iMc]   = trueParticle.Position().X();
    fMCParticleStartPositionY[iMc]   = trueParticle.Position().Y();
    fMCParticleStartPositionZ[iMc]   = trueParticle.Position().Z();
    fMCParticleStartPositionT[iMc]   = trueParticle.Position().T();
    fMCParticleEndPositionX[iMc]     = trueParticle.EndPosition().X();
    fMCParticleEndPositionY[iMc]     = trueParticle.EndPosition().Y();
    fMCParticleEndPositionZ[iMc]     = trueParticle.EndPosition().Z();
    fMCParticleEndPositionT[iMc]     = trueParticle.EndPosition().T();
    fMCParticleStartMomentumX[iMc]   = trueParticle.Momentum().X();
    fMCParticleStartMomentumY[iMc]   = trueParticle.Momentum().Y();
    fMCParticleStartMomentumZ[iMc]   = trueParticle.Momentum().Z();
    fMCParticleStartMomentumE[iMc]   = trueParticle.Momentum().E();
    fMCParticleEndMomentumX[iMc]     = trueParticle.EndMomentum().X();
    fMCParticleEndMomentumY[iMc]     = trueParticle.EndMomentum().Y();
    fMCParticleEndMomentumZ[iMc]     = trueParticle.EndMomentum().Z();
    fMCParticleEndMomentumE[iMc]     = trueParticle.EndMomentum().E();
  }

  return true;
}

void test::refactoredPandoraAna::parseTrack(art::Event const& e, const art::Ptr<recob::PFParticle> &pfp, int iPfp) {
  fPFPIsTrack[iPfp] = true;
  art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, e, fPFParticleLabel, fTrackLabel);

  fPFPTrackID[iPfp]               = track->ID();
  fPFPTrackLength[iPfp]           = track->Length();
  fPFPTrackStartX[iPfp]           = track->Start().X();
  fPFPTrackStartY[iPfp]           = track->Start().Y();
  fPFPTrackStartZ[iPfp]           = track->Start().Z();
  fPFPTrackVertexX[iPfp]          = track->Vertex().X();
  fPFPTrackVertexY[iPfp]          = track->Vertex().Y();
  fPFPTrackVertexZ[iPfp]          = track->Vertex().Z();
  fPFPTrackEndX[iPfp]             = track->End().X();
  fPFPTrackEndY[iPfp]             = track->End().Y();
  fPFPTrackEndZ[iPfp]             = track->End().Z();
  fPFPTrackTheta[iPfp]            = track->Theta();
  fPFPTrackPhi[iPfp]              = track->Phi();
  fPFPTrackZenithAngle[iPfp]      = track->ZenithAngle();
  fPFPTrackAzimuthAngle[iPfp]     = track->AzimuthAngle();
  fPFPTrackStartDirectionX[iPfp]  = track->StartDirection().X();
  fPFPTrackStartDirectionY[iPfp]  = track->StartDirection().Y();
  fPFPTrackStartDirectionZ[iPfp]  = track->StartDirection().Z();
  fPFPTrackVertexDirectionX[iPfp] = track->VertexDirection().X();
  fPFPTrackVertexDirectionY[iPfp] = track->VertexDirection().Y();
  fPFPTrackVertexDirectionZ[iPfp] = track->VertexDirection().Z();
  fPFPTrackEndDirectionX[iPfp]    = track->EndDirection().X();
  fPFPTrackEndDirectionY[iPfp]    = track->EndDirection().Y();
  fPFPTrackEndDirectionZ[iPfp]    = track->EndDirection().Z();
  fPFPTrackChi2[iPfp]             = track->Chi2();
  fPFPTrackNdof[iPfp]             = track->Ndof();

  std::vector<art::Ptr<recob::Hit>> pfpHits;
  pfpHits = dune_ana::DUNEAnaTrackUtils::GetHits(track, e, fTrackLabel);
  fPFPNHits[iPfp] = pfpHits.size();

  std::vector<art::Ptr<recob::Hit> > pfpHitsView0 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, e, fPFParticleLabel, 0);
  fPFPNHitsView[iPfp][0] = pfpHitsView0.size();

  std::vector<art::Ptr<recob::Hit> > pfpHitsView1 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, e, fPFParticleLabel, 1);
  fPFPNHitsView[iPfp][1] = pfpHitsView1.size();

  std::vector<art::Ptr<recob::Hit> > pfpHitsView2 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, e, fPFParticleLabel, 2);
  fPFPNHitsView[iPfp][2] = pfpHitsView2.size();

  if(e.isRealData())
    return;

  // Get total hit matching stats...
  TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(fClockData, pfpHits, fRollUpUnsavedIDs));
  if (TruthMatchUtils::Valid(g4ID)) {
    fPFPTrueParticleMatchedID[iPfp] = g4ID;

    int pos(999999);
    for(int unsigned ipos = 0; ipos < fNMCParticles; ipos++)
      if(fMCParticleTrackID[ipos] == g4ID)
        pos = ipos;

    fPFPTrueParticleMatchedPosition[iPfp] = pos;
  }

  // Get view0 hit matching stats...
  TruthMatchUtils::G4ID g4IDView0(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(fClockData, pfpHitsView0, fRollUpUnsavedIDs));
  if (TruthMatchUtils::Valid(g4ID)) {
    fPFPTrueParticleMatchedIDView[iPfp][0] = g4IDView0;

    for(int unsigned ipos = 0; ipos < fNMCParticles; ipos++) {
      if (fMCParticleTrackID[ipos] != g4IDView0)
        continue;

      fPFPTrueParticleMatchedPositionView[iPfp][0] = ipos;
      break;
    }
  }

  // Get view1 hit matching stats...
  TruthMatchUtils::G4ID g4IDView1(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(fClockData, pfpHitsView1, fRollUpUnsavedIDs));
  if (TruthMatchUtils::Valid(g4ID)) {
    fPFPTrueParticleMatchedIDView[iPfp][1] = g4IDView1;

    for(int unsigned ipos = 0; ipos < fNMCParticles; ipos++) {
      if (fMCParticleTrackID[ipos] != g4IDView1)
        continue;

      fPFPTrueParticleMatchedPositionView[iPfp][1] = ipos;
      break;
    }
  }

  // Get view2 hit matching stats...
  TruthMatchUtils::G4ID g4IDView2(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(fClockData, pfpHitsView2, fRollUpUnsavedIDs));
  if (TruthMatchUtils::Valid(g4ID)) {
    fPFPTrueParticleMatchedIDView[iPfp][2] = g4IDView2;

    for(int unsigned ipos = 0; ipos < fNMCParticles; ipos++) {
      if (fMCParticleTrackID[ipos] != g4IDView2)
        continue;

      fPFPTrueParticleMatchedPositionView[iPfp][2] = ipos;
      break;
    }
  }

  return;
}

void test::refactoredPandoraAna::parseShower(art::Event const& e, const art::Ptr<recob::PFParticle> &pfp, int iPfp) {
  fPFPIsShower[iPfp] = true;
  art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp, e, fPFParticleLabel, fShowerLabel);

  fPFPShowerID[iPfp]            =      shower->ID();
  fPFPShowerBestPlane[iPfp]     =      shower->best_plane();
  fPFPShowerDirectionX[iPfp]    =      shower->Direction().X();
  fPFPShowerDirectionY[iPfp]    =      shower->Direction().Y();
  fPFPShowerDirectionZ[iPfp]    =      shower->Direction().Z();
  fPFPShowerDirectionErrX[iPfp] =      shower->DirectionErr().X();
  fPFPShowerDirectionErrY[iPfp] =      shower->DirectionErr().Y();
  fPFPShowerDirectionErrZ[iPfp] =      shower->DirectionErr().Z();
  fPFPShowerStartX[iPfp]        =      shower->ShowerStart().X();
  fPFPShowerStartY[iPfp]        =      shower->ShowerStart().Y();
  fPFPShowerStartZ[iPfp]        =      shower->ShowerStart().Z();
  fPFPShowerStartErrX[iPfp]     =      shower->ShowerStartErr().X();
  fPFPShowerStartErrY[iPfp]     =      shower->ShowerStartErr().Y();
  fPFPShowerStartErrZ[iPfp]     =      shower->ShowerStartErr().Z();
  fPFPShowerLength[iPfp]        =      shower->Length();
  fPFPShowerOpenAngle[iPfp]     =      shower->OpenAngle();

  std::vector<art::Ptr<recob::Hit>> pfpHits;
  pfpHits = dune_ana::DUNEAnaShowerUtils::GetHits(shower,e,fShowerLabel);
  fPFPNHits[iPfp] = pfpHits.size();
  std::vector<art::Ptr<recob::Hit> > pfpHitsView0 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, e, fPFParticleLabel, 0);
  fPFPNHitsView[iPfp][0] = pfpHitsView0.size();
  std::vector<art::Ptr<recob::Hit> > pfpHitsView1 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, e, fPFParticleLabel, 1);
  fPFPNHitsView[iPfp][1] = pfpHitsView1.size();
  std::vector<art::Ptr<recob::Hit> > pfpHitsView2 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, e, fPFParticleLabel, 2);
  fPFPNHitsView[iPfp][2] = pfpHitsView2.size();

  if(e.isRealData())
    return;

  // Get total hit matching stats...
  TruthMatchUtils::G4ID g4ID(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(fClockData, pfpHits, fRollUpUnsavedIDs));
  if (TruthMatchUtils::Valid(g4ID)) {
    fPFPTrueParticleMatchedID[iPfp] = g4ID;

    int pos(999999);
    for(int unsigned ipos = 0; ipos < fNMCParticles; ipos++)
      if(fMCParticleTrackID[ipos] == g4ID)
        pos=ipos;

    fPFPTrueParticleMatchedPosition[iPfp] = pos;
  }

  // Get view0 hit matching stats...
  TruthMatchUtils::G4ID g4IDView0(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(fClockData, pfpHitsView0, fRollUpUnsavedIDs));
  if (TruthMatchUtils::Valid(g4ID)) {
    fPFPTrueParticleMatchedIDView[iPfp][0] = g4IDView0;

    for(int unsigned ipos = 0; ipos < fNMCParticles; ipos++) {
      if (fMCParticleTrackID[ipos] != g4IDView0)
        continue;

      fPFPTrueParticleMatchedPositionView[iPfp][0] = ipos;
      break;
    }
  }

  // Get view1 hit matching stats...
  TruthMatchUtils::G4ID g4IDView1(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(fClockData, pfpHitsView1, fRollUpUnsavedIDs));
  if (TruthMatchUtils::Valid(g4ID)) {
    fPFPTrueParticleMatchedIDView[iPfp][1] = g4IDView1;

    for(int unsigned ipos = 0; ipos < fNMCParticles; ipos++) {
      if (fMCParticleTrackID[ipos] != g4IDView1)
        continue;

      fPFPTrueParticleMatchedPositionView[iPfp][1] = ipos;
      break;
    }
  }

  // Get view2 hit matching stats...
  TruthMatchUtils::G4ID g4IDView2(TruthMatchUtils::TrueParticleIDFromTotalRecoHits(fClockData, pfpHitsView2, fRollUpUnsavedIDs));
  if (TruthMatchUtils::Valid(g4ID)) {
    fPFPTrueParticleMatchedIDView[iPfp][2] = g4IDView2;

    for(int unsigned ipos = 0; ipos < fNMCParticles; ipos++) {
      if (fMCParticleTrackID[ipos] != g4IDView2)
        continue;

      fPFPTrueParticleMatchedPositionView[iPfp][2] = ipos;
      break;
    }
  }

  return;
}

void test::refactoredPandoraAna::recoMCInfo(art::Event const& e, const art::Ptr<recob::PFParticle> &pfp, int iPfp) {

  if(e.isRealData())
    return;

  // Can't find the MC match, can't do anything.
  if (fPFPTrueParticleMatchedPosition[iPfp] == 999999)
    return;

  std::vector<art::Ptr<recob::Hit>> particleHits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, e, fHitLabel);
  std::vector<art::Ptr<recob::Hit>> particleHits0 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, e, fHitLabel, 0);
  std::vector<art::Ptr<recob::Hit>> particleHits1 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, e, fHitLabel, 1);
  std::vector<art::Ptr<recob::Hit>> particleHits2 = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, e, fHitLabel, 2);

  int truthId = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(fClockData, particleHits, fRollUpUnsavedIDs);

  fPFPCompleteness[iPfp]        = this->getCompleteness(particleHits,      truthId);
  fPFPCompletenessView[iPfp][0] = this->getCompleteness(particleHits0, truthId, 0);
  fPFPCompletenessView[iPfp][1] = this->getCompleteness(particleHits1, truthId, 1);
  fPFPCompletenessView[iPfp][2] = this->getCompleteness(particleHits2, truthId, 2);

  fPFPPurity[iPfp]        = this->getPurity(particleHits,      truthId);
  fPFPPurityView[iPfp][0] = this->getPurity(particleHits0, truthId, 0);
  fPFPPurityView[iPfp][1] = this->getPurity(particleHits1, truthId, 1);
  fPFPPurityView[iPfp][2] = this->getPurity(particleHits2, truthId, 2);

  if(fMCParticlePdgCode[fPFPTrueParticleMatchedPosition[iPfp]] > 0 && fMCParticlePdgCode[fPFPTrueParticleMatchedPosition[iPfp]] < 999999)
    fPFPMatchedMCPDG[iPfp] = fMCParticlePdgCode[fPFPTrueParticleMatchedPosition[iPfp]];
  if(fMCParticleTrueEnergy[fPFPTrueParticleMatchedPosition[iPfp]] > 0 && fMCParticleTrueEnergy[fPFPTrueParticleMatchedPosition[iPfp]] < 999999)
    fPFPMatchedMCEnergy[iPfp] = fMCParticleTrueEnergy[fPFPTrueParticleMatchedPosition[iPfp]];
  if(fMCParticleNHits[fPFPTrueParticleMatchedPosition[iPfp]] > 0 && fMCParticleNHits[fPFPTrueParticleMatchedPosition[iPfp]] < 999999)
    fPFPMatchedMCNHits[iPfp] = fMCParticleNHits[fPFPTrueParticleMatchedPosition[iPfp]];

  return;
}

bool test::refactoredPandoraAna::parsePFParticles(art::Event const& e) {

  // Access the PFParticles from Pandora.
  const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect = dune_ana::DUNEAnaEventUtils::GetPFParticles(e, fPFParticleLabel);
  fNPFParticles = pfparticleVect.size();

  if(! fNPFParticles) {
    std::cout << "No PFParticles found!" << std::endl;
    fTree->Fill();
    return false;
  }

  //Access the Clusters
  art::Handle<std::vector<recob::Cluster>> clusterHandle;
  std::vector<art::Ptr<recob::Cluster>> clusterVect;

  if (e.getByLabel(fPFParticleLabel,clusterHandle))
    art::fill_ptr_vector(clusterVect,clusterHandle);

  art::FindManyP<recob::Cluster> clusterParticleAssoc(pfparticleVect, e, fPFParticleLabel);
  art::Handle<std::vector<recob::Track>> trackHandle;

  if (! (e.getByLabel(fTrackLabel, trackHandle))) {
    std::cout << "Unable to find std::vector<recob::Track> with module label: " << fTrackLabel << std::endl;
    fTree->Fill();
    return false;
  }

  std::vector<art::Ptr<recob::Track> > trackList;
  if (e.getByLabel(fTrackLabel,trackHandle))
    art::fill_ptr_vector(trackList, trackHandle);

  int iPfp(0);
  for(const art::Ptr<recob::PFParticle> &pfp: pfparticleVect) {

    fPFPID[iPfp]        = iPfp;
    fPFPIsPrimary[iPfp] = pfp->IsPrimary();
    fPFPPdgCode[iPfp]   = pfp->PdgCode();
    fPFPNChildren[iPfp] = pfp->NumDaughters();
    (pfp->IsPrimary()) ? fPFPParentID[iPfp] = -1 : fPFPParentID[iPfp] = pfp->Parent();

    std::vector<art::Ptr<recob::Cluster>> pfpClusters = clusterParticleAssoc.at(pfp.key());
    fPFPNClusters[iPfp] = pfpClusters.size();

    if(!pfpClusters.empty()) {
      int iClu(0);
      for(const art::Ptr<recob::Cluster> &clu:pfpClusters) {

        fPFPCluPlane[iPfp][iClu]    = clu->Plane().asPlaneID().Plane;
        fPFPCluView[iPfp][iClu]     = clu->View();
        fPFPCluNHits[iPfp][iClu]    = clu->NHits();
        fPFPCluIntegral[iPfp][iClu] = clu->Integral();
        iClu++;

        if (iClu == kNMaxPFPClusters)
          break;
      }
    }

    if (dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, e, fPFParticleLabel, fTrackLabel))
      this->parseTrack(e, pfp, iPfp);

    if (dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp, e, fPFParticleLabel, fShowerLabel))
      this->parseShower(e, pfp, iPfp);

    this->recoMCInfo(e, pfp, iPfp);

    iPfp++;
  }

  return true;
}

void test::refactoredPandoraAna::analyze(art::Event const& e) {
  // Reset to start.
  reset();

  // Get event level information, before processing...
  const art::ServiceHandle<cheat::BackTrackerService> btServ;
  fEventID = e.id().event();
  fRunID = e.id().run();
  fSubRunID = e.id().subRun();

  std::cout << "============== EVENT ID: " << fEventID << " == RUN ID: " << fRunID << " == SUBRUN ID: " << fSubRunID << " ================" << std::endl;

  // Access the truth information
  if(e.isRealData()) {
    std::cout << "Can't plot for real data!" << std::endl;
    return;
  }

  // Get MCTruth information
  bool mcTruthParsed = this->processMCTruth(e);

  if (! mcTruthParsed) {
    std::cout << "No MCTruth found!" << std::endl;
    fTree->Fill();
    return;
  }

  // Get MCParticle information
  bool mcParticlesParsed = this->processMCParticles(e);

  if (! mcParticlesParsed) {
    std::cout << "No MCTruth found!" << std::endl;
    fTree->Fill();
    return;
  }

  bool pfpParticlesParsed = this->parsePFParticles(e);

  if (! pfpParticlesParsed) {
    std::cout << "Failed parsing the PFParticles!" << std::endl;
    fTree->Fill();
    return;
  }

  fTree->Fill();
}

void test::refactoredPandoraAna::beginJob() {

  // Deep clean the variables
  reset(true);

  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("pandoraOutput","Pandora Output Tree");

  // Event branches
  fTree->Branch("eventID",&fEventID,"eventID/i");
  fTree->Branch("runID",&fRunID,"runID/i");
  fTree->Branch("subrunID",&fSubRunID,"subrunID/i");
  fTree->Branch("nMCParticles",&fNMCParticles,"nMCParticles/i");
  fTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");

  // MC truth branches
  fTree->Branch("nuPdgCodeTruth", &fNuPdgCodeTruth, "nuPdgCodeTruth/i");
  fTree->Branch("nuCCNCTruth", &fNuCCNCTruth, "nuCCNCTruth/i");
  fTree->Branch("nuModeTruth", &fNuModeTruth, "nuModeTruth/i");
  fTree->Branch("nuETruth", &fNuETruth, "nuETruth/d");
  fTree->Branch("nuVertexXTruth", &fNuVertexXTruth, "nuVertexXTruth/d");
  fTree->Branch("nuVertexYTruth", &fNuVertexYTruth, "nuVertexYTruth/d");
  fTree->Branch("nuVertexZTruth", &fNuVertexZTruth, "nuVertexZTruth/d");

  // MC particle branches
  fTree->Branch("mcIsMCPrimary",&fMCIsPrimary,"MCIsPrimary[nMCParticles]/O");
  fTree->Branch("mcParticlePdgCode",&fMCParticlePdgCode,"MCParticlePdgCode[nMCParticles]/I");
  fTree->Branch("mcParticleTrueEnergy",&fMCParticleTrueEnergy,"MCParticleTrueEnergy[nMCParticles]/D");
  fTree->Branch("mcParticleTrackID",&fMCParticleTrackID,"MCParticleTrackID[nMCParticles]/I");
  fTree->Branch("mcParticleParentTrackID",&fMCParticleParentTrackID,"MCParticleParentTrackID[nMCParticles]/I");
  fTree->Branch("mcParticleStartProcess",&fMCParticleStartProcess,"MCParticleStartProcess[nMCParticles]/C");
  fTree->Branch("mcParticleEndProcess",&fMCParticleEndProcess, "MCParticleEndProcess[nMCParticles]/C");
  fTree->Branch("mcParticleNTrajectoryPoints",&fMCParticleNTrajectoryPoint, "MCParticleNTrajectoryPoint[nMCParticles]/I");
  fTree->Branch("mcParticleStartPositionX",&fMCParticleStartPositionX,"MCParticleStartPositionX[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionY",&fMCParticleStartPositionY,"MCParticleStartPositionY[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionZ",&fMCParticleStartPositionZ,"MCParticleStartPositionZ[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionT",&fMCParticleStartPositionT,"MCParticleStartPositionT[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumX",&fMCParticleStartMomentumX, "MCParticleStartMomentumX[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumY",&fMCParticleStartMomentumY, "MCParticleStartMomentumY[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumZ",&fMCParticleStartMomentumZ, "MCParticleStartMomentumZ[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumE",&fMCParticleStartMomentumE, "MCParticleStartMomentumE[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionX",&fMCParticleEndPositionX, "MCParticleEndPositionX[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionY",&fMCParticleEndPositionY, "MCParticleEndPositionY[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionZ",&fMCParticleEndPositionZ, "MCParticleEndPositionZ[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionT",&fMCParticleEndPositionT, "MCParticleEndPositionT[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumX",&fMCParticleEndMomentumX, "MCParticleEndMomentumX[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumY",&fMCParticleEndMomentumY, "MCParticleEndMomentumY[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumZ",&fMCParticleEndMomentumZ, "MCParticleEndMomentumZ[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumE",&fMCParticleEndMomentumE, "MCParticleEndMomentumE[nMCParticles]/D");
  fTree->Branch("mcParticleVertexTime",&fMCParticleVertexTime, "MCParticleVertexTime[nMCParticles]/D");
  fTree->Branch("mcParticleEndTime",&fMCParticleEndTime, "MCParticleEndTime[nMCParticles]/D");
  fTree->Branch("mcParticleNHits", &fMCParticleNHits, "MCParticleNHits[nMCParticles]/I");
  fTree->Branch("mcParticleNHitsView", &fMCParticleNHitsView, "MCParticleNHitsView[nMCParticles][3]/I");

  // PFP branches
  fTree->Branch("pfpTrueParticleMatchedID",&fPFPTrueParticleMatchedID,"PFPTrueParticleMatchedID[nPFParticles]/I");
  fTree->Branch("pfpTrueParticleMatchedPosition",&fPFPTrueParticleMatchedPosition,"PFPTrueParticleMatchedPosition[nPFParticles]/I");
  fTree->Branch("pfpIsPrimary",&fPFPIsPrimary,"PFPIsPrimary[nPFParticles]/O");
  fTree->Branch("pfpID",&fPFPID, "PFPID[nPFParticles]/I");
  fTree->Branch("pfpParentID",&fPFPParentID, "PFPParentID[nPFParticles]/I");
  fTree->Branch("pfpPdgCode",&fPFPPdgCode, "PFPPdgCode[nPFParticles]/I");
  fTree->Branch("pfpNChildren",&fPFPNChildren,"PFPNChildren[nPFParticles]/I");
  fTree->Branch("pfpNClusters",&fPFPNClusters,"PFPNClusters[nPFParticles]/I");
  fTree->Branch("pfpNHits",&fPFPNHits,"PFPNHits[nPFParticles]/I");
  fTree->Branch("pfpNHitsView",&fPFPNHitsView,"PFPNHitsView[nPFParticles][3]/I");
  fTree->Branch("pfpNSharedTrueParticleHits",&fPFPNSharedTrueParticleHits,"PFPNSharedTrueParticleHits[nPFParticles]/I");
  fTree->Branch("pfpNSharedTrueParticleHitsView",&fPFPNSharedTrueParticleHitsView,"PFPNSharedTrueParticleHitsView[nPFParticles][3]/I");
  fTree->Branch("pfpTrueParticleMatchedIDView",&fPFPTrueParticleMatchedIDView,"PFPTrueParticleMatchedIDView[nPFParticles][3]/I");
  fTree->Branch("pfpTrueParticleMatchedPositionView",&fPFPTrueParticleMatchedPositionView,"PFPTrueParticleMatchedPositionView[nPFParticles][3]/I");
  fTree->Branch("pfpIsTrack",&fPFPIsTrack,"PFPIsTrack[nPFParticles]/O");
  fTree->Branch("pfpIsShower",&fPFPIsShower,"PFPIsShower[nPFParticles]/O");
  fTree->Branch("pfpTrackID", &fPFPTrackID,"PFPNClusters[nPFParticles]/I");
  fTree->Branch("pfpTrackLength",&fPFPTrackLength,"PFPTrackLength[nPFParticles]/D");
  fTree->Branch("pfpTrackStartX",&fPFPTrackStartX,"PFPTrackStartX[nPFParticles]/D");
  fTree->Branch("pfpTrackStartY",&fPFPTrackStartY,"PFPTrackStartY[nPFParticles]/D");
  fTree->Branch("pfpTrackStartZ",&fPFPTrackStartZ,"PFPTrackStartZ[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexX",&fPFPTrackVertexX,"PFPTrackVertexX[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexY",&fPFPTrackVertexY, "PFPTrackVertexY[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexZ", &fPFPTrackVertexZ,"PFPTrackVertexZ[nPFParticles]/D");
  fTree->Branch("pfpTrackEndX",&fPFPTrackEndX,"PFPTrackEndX[nPFParticles]/D");
  fTree->Branch("pfpTrackEndY",&fPFPTrackEndY,"PFPTrackEndY[nPFParticles]/D");
  fTree->Branch("pfpTrackEndZ",&fPFPTrackEndZ,"PFPTrackEndZ[nPFParticles]/D");
  fTree->Branch("pfpTrackTheta",&fPFPTrackTheta,"PFPTrackTheta[nPFParticles]/D");
  fTree->Branch("pfpTrackPhi",&fPFPTrackPhi,"PFPTrackPhi[nPFParticles]/D");
  fTree->Branch("pfpTrackZenithAngle",&fPFPTrackZenithAngle,"PFPTrackZenithAngle[nPFParticles]/D");
  fTree->Branch("pfpTrackAzimuthAngle",&fPFPTrackAzimuthAngle,"PFPTrackAzimuthAngle[nPFParticles]/D");
  fTree->Branch("pfpTrackStartDirectionX",&fPFPTrackStartDirectionX,"PFPTrackStartDirectionX[nPFParticles]/D");
  fTree->Branch("pfpTrackStartDirectionY",&fPFPTrackStartDirectionY,"PFPTrackStartDirectionY[nPFParticles]/D");
  fTree->Branch("pfpTrackStartDirectionZ",&fPFPTrackStartDirectionZ,"PFPTrackStartDirectionZ[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexDirectionX",&fPFPTrackVertexDirectionX,"PFPTrackVertexDirectionX[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexDirectionY",&fPFPTrackVertexDirectionY,"PFPTrackVertexDirectionY[nPFParticles]/D");
  fTree->Branch("pfpTrackVertexDirectionZ",&fPFPTrackVertexDirectionZ,"PFPTrackVertexDirectionZ[nPFParticles]/D");
  fTree->Branch("pfpTrackEndDirectionX",&fPFPTrackEndDirectionX,"PFPTrackEndDirectionX[nPFParticles]/D");
  fTree->Branch("pfpTrackEndDirectionY",&fPFPTrackEndDirectionY,"PFPTrackEndDirectionY[nPFParticles]/D");
  fTree->Branch("pfpTrackEndDirectionZ",&fPFPTrackEndDirectionZ,"PFPTrackEndDirectionZ[nPFParticles]/D");
  fTree->Branch("pfpTrackChi2",&fPFPTrackChi2,"PFPTrackChi2[nPFParticles]/F");
  fTree->Branch("pfpTrackStartNdof",&fPFPTrackNdof,"PFPTrackNdof[nPFParticles]/I");

  fTree->Branch("pfpCluPlane",fPFPCluPlane,"PFPCluPlane[nPFParticles][100]/I");
  fTree->Branch("pfpCluView",fPFPCluView,"PFPCluView[nPFParticles][100]/I");
  fTree->Branch("pfpCluNHits",fPFPCluNHits,"PFPCluNHits[nPFParticles][100]/I");
  fTree->Branch("pfpCluIntegral",fPFPCluIntegral,"PFPCluIntegral[nPFParticles][100]/D");

  fTree->Branch("pfpShowerID",&fPFPShowerID,"PFPShowerID[nPFParticles]/I");
  fTree->Branch("pfpShowerBestPlane",&fPFPShowerBestPlane,"PFPShowerBestPlane[nPFParticles]/I");
  fTree->Branch("pfpShowerDirectionX",&fPFPShowerDirectionX,"PFPShowerDirectionX[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionY",&fPFPShowerDirectionY,"PFPShowerDirectionY[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionZ",&fPFPShowerDirectionZ,"PFPShowerDirectionZ[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionErrX",&fPFPShowerDirectionErrX,"PFPShowerDirectionErrX[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionErrY",&fPFPShowerDirectionErrY,"PFPShowerDirectionErrY[nPFParticles]/D");
  fTree->Branch("pfpShowerDirectionErrZ",&fPFPShowerDirectionErrZ,"PFPShowerDirectionErrZ[nPFParticles]/D");
  fTree->Branch("pfpShowerStartX",&fPFPShowerStartX,"PFPShowerStartX[nPFParticles]/D");
  fTree->Branch("pfpShowerStartY",&fPFPShowerStartY,"PFPShowerStartY[nPFParticles]/D");
  fTree->Branch("pfpShowerStartZ",&fPFPShowerStartZ,"PFPShowerStartZ[nPFParticles]/D");
  fTree->Branch("pfpShowerStartErrX",&fPFPShowerStartErrX,"PFPShowerStartErrX[nPFParticles]/D");
  fTree->Branch("pfpShowerStartErrY",&fPFPShowerStartErrY,"PFPShowerStartErrY[nPFParticles]/D");
  fTree->Branch("pfpShowerStartErrZ",&fPFPShowerStartErrZ,"PFPShowerStartErrZ[nPFParticles]/D");
  fTree->Branch("pfpShowerLength",&fPFPShowerLength,"PFPShowerLength[nPFParticles]/D");
  fTree->Branch("pfpShowerOpeningAngle",&fPFPShowerOpenAngle,"PFPShowerOpenAngle[nPFParticles]/D");

  fTree->Branch("pfpMatchedMCPDG", &fPFPMatchedMCPDG, "PFPMatchedMCPDG[nPFParticles]/D");
  fTree->Branch("pfpMatchedMCEnergy", &fPFPMatchedMCEnergy, "PFPMatchedMCEnergy[nPFParticles]/D");
  fTree->Branch("pfpMatchedMCNHits", &fPFPMatchedMCNHits, "PFPMatchedMCNHits[nPFParticles]/D");

  fTree->Branch("pfpCompleteness", &fPFPCompleteness, "PFPCompleteness[nPFParticles]/D");
  fTree->Branch("pfpCompletenessView", &fPFPCompletenessView, "PFPCompletenessView[nPFParticles][3]/D");
  fTree->Branch("pfpPurity", &fPFPPurity, "PFPPurity[nPFParticles]/D");
  fTree->Branch("pfpPurityView", &fPFPPurityView, "PFPPurityView[nPFParticles][3]/D");
}

void test::refactoredPandoraAna::endJob()
{
  // Implementation of optional member function here.
}

void test::refactoredPandoraAna::reset(bool deepClean) {

    fNuPdgCodeTruth = 999999;
    fNuModeTruth    = 999999;
    fNuCCNCTruth    = 999999;
    fNuETruth       = 999999;
    fNuVertexXTruth = 999999;
    fNuVertexYTruth = 999999;
    fNuVertexZTruth = 999999;

    for (unsigned int iMc = 0; iMc < (deepClean ? kNMaxMCParticles : fNMCParticles); iMc++) {
      fMCIsPrimary[iMc]                = 0;
      fMCParticlePdgCode[iMc]          = 999999;
      fMCParticleTrueEnergy[iMc]       = 999999;
      fMCParticleTrackID[iMc]          = 999999;
      fMCParticleParentTrackID[iMc]    = 999999;
      fMCParticleStartProcess[iMc]     = "";
      fMCParticleEndProcess[iMc]       = "";
      fMCParticleNTrajectoryPoint[iMc] = 999999;
      fMCParticleStartPositionX[iMc]   = 999999;
      fMCParticleStartPositionY[iMc]   = 999999;
      fMCParticleStartPositionZ[iMc]   = 999999;
      fMCParticleStartPositionT[iMc]   = 999999;
      fMCParticleStartMomentumX[iMc]   = 999999;
      fMCParticleStartMomentumY[iMc]   = 999999;
      fMCParticleStartMomentumZ[iMc]   = 999999;
      fMCParticleStartMomentumE[iMc]   = 999999;
      fMCParticleEndPositionX[iMc]     = 999999;
      fMCParticleEndPositionY[iMc]     = 999999;
      fMCParticleEndPositionZ[iMc]     = 999999;
      fMCParticleEndPositionT[iMc]     = 999999;
      fMCParticleEndMomentumX[iMc]     = 999999;
      fMCParticleEndMomentumY[iMc]     = 999999;
      fMCParticleEndMomentumZ[iMc]     = 999999;
      fMCParticleEndMomentumE[iMc]     = 999999;
      fMCParticleVertexTime[iMc]       = 999999;
      fMCParticleEndTime[iMc]          = 999999;
      fMCParticleNHits[iMc]            = 999999;
      fMCParticleNHitsView[iMc][0]     = 999999;
      fMCParticleNHitsView[iMc][1]     = 999999;
      fMCParticleNHitsView[iMc][2]     = 999999;
    }

    fNMCParticles = 0;

    for(unsigned int iPfp = 0; iPfp < (deepClean ? kNMaxPFParticles : fNPFParticles); iPfp++) {

      fPFPID[iPfp]                          = 999999;
      fPFPTrueParticleMatchedID[iPfp]       = 999999;
      fPFPTrueParticleMatchedPosition[iPfp] = 999999;

      fPFPNHits[iPfp]                       = 999999;
      fPFPNSharedTrueParticleHits[iPfp]     = 999999;

      fPFPNClusters[iPfp]                   = 999999;
      fPFPIsTrack[iPfp]                     = 0;
      fPFPIsShower[iPfp]                    = 0;

      fPFPTrackID[iPfp]                     = 999999;
      fPFPTrackLength[iPfp]                 = 999999;
      fPFPTrackStartX[iPfp]                 = 999999;
      fPFPTrackStartY[iPfp]                 = 999999;
      fPFPTrackStartZ[iPfp]                 = 999999;
      fPFPTrackVertexX[iPfp]                = 999999;
      fPFPTrackVertexY[iPfp]                = 999999;
      fPFPTrackVertexZ[iPfp]                = 999999;
      fPFPTrackEndX[iPfp]                   = 999999;
      fPFPTrackEndY[iPfp]                   = 999999;
      fPFPTrackEndZ[iPfp]                   = 999999;
      fPFPTrackTheta[iPfp]                  = 999999;
      fPFPTrackPhi[iPfp]                    = 999999;
      fPFPTrackZenithAngle[iPfp]            = 999999;
      fPFPTrackAzimuthAngle[iPfp]           = 999999;
      fPFPTrackStartDirectionX[iPfp]        = 999999;
      fPFPTrackStartDirectionY[iPfp]        = 999999;
      fPFPTrackStartDirectionZ[iPfp]        = 999999;
      fPFPTrackVertexDirectionX[iPfp]       = 999999;
      fPFPTrackVertexDirectionY[iPfp]       = 999999;
      fPFPTrackVertexDirectionZ[iPfp]       = 999999;
      fPFPTrackEndDirectionX[iPfp]          = 999999;
      fPFPTrackEndDirectionY[iPfp]          = 999999;
      fPFPTrackEndDirectionZ[iPfp]          = 999999;
      fPFPTrackChi2[iPfp]                   = 999999;
      fPFPTrackNdof[iPfp]                   = 999999;

      fPFPShowerID[iPfp]                    = 999999;
      fPFPShowerBestPlane[iPfp]             = 999999;
      fPFPShowerDirectionX[iPfp]            = 999999;
      fPFPShowerDirectionY[iPfp]            = 999999;
      fPFPShowerDirectionZ[iPfp]            = 999999;
      fPFPShowerDirectionErrX[iPfp]         = 999999;
      fPFPShowerDirectionErrY[iPfp]         = 999999;
      fPFPShowerDirectionErrZ[iPfp]         = 999999;
      fPFPShowerStartX[iPfp]                = 999999;
      fPFPShowerStartY[iPfp]                = 999999;
      fPFPShowerStartZ[iPfp]                = 999999;
      fPFPShowerStartErrX[iPfp]             = 999999;
      fPFPShowerStartErrY[iPfp]             = 999999;
      fPFPShowerStartErrZ[iPfp]             = 999999;
      fPFPShowerLength[iPfp]                = 999999;
      fPFPShowerOpenAngle[iPfp]             = 999999;

      for (int iClu=0; iClu<kNMaxPFPClusters; iClu++) {
        fPFPCluPlane[iPfp][iClu]    = 999999;
        fPFPCluView[iPfp][iClu]     = 999999;
        fPFPCluNHits[iPfp][iClu]    = 999999;
        fPFPCluIntegral[iPfp][iClu] = 999999;
      }

      fPFPMatchedMCPDG[iPfp]    = 999999;
      fPFPMatchedMCEnergy[iPfp] = 999999;
      fPFPMatchedMCNHits[iPfp]  = 999999;

      fPFPCompleteness[iPfp]    = 999999;
      fPFPPurity[iPfp]          = 999999;

      for (unsigned int iView = 0; iView < kNViews; iView++) {
        fPFPNHitsView[iPfp][iView]                       = 999999;
        fPFPNSharedTrueParticleHitsView[iPfp][iView]     = 999999;
        fPFPTrueParticleMatchedIDView[iPfp][iView]       = 999999;
        fPFPTrueParticleMatchedPositionView[iPfp][iView] = 999999;
        fPFPPurityView[iPfp][iView]                      = 999999;
        fPFPCompletenessView[iPfp][iView]                = 999999;
      }
    }

    fNPFParticles = 0;

    return;
}

DEFINE_ART_MODULE(test::refactoredPandoraAna)
