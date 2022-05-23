////////////////////////////////////////////////////////////////////////
// Class:       mcShowerPandoraAna
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
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "larpandoracontent/LArObjects/LArPfoObjects.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include <TTree.h>
#include <vector>
#include <string>
#include <TH1D.h>
#include <TLorentzVector.h>

namespace test {
  class mcShowerPandoraAna;
}

typedef std::map<int, std::vector<pandora::CartesianVector>> TrackIDToIDEsMap;

class test::mcShowerPandoraAna : public art::EDAnalyzer {
public:
  explicit mcShowerPandoraAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  mcShowerPandoraAna(mcShowerPandoraAna const&) = delete;
  mcShowerPandoraAna(mcShowerPandoraAna&&) = delete;
  mcShowerPandoraAna& operator=(mcShowerPandoraAna const&) = delete;
  mcShowerPandoraAna& operator=(mcShowerPandoraAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  // Functions per concern
  bool processMCTruth(art::Event const& e);
  bool processMCParticles(art::Event const& e);

  // Get shower information
  void getMCShowerInformation(art::Event const& e, TrackIDToIDEsMap &idToIDEs, const simb::MCParticle &mcParticle, const int iMc);

  void reset(bool deepClean = false);

private:

  TTree *fTree;

  unsigned int fEventID;
  unsigned int fRunID;
  unsigned int fSubRunID;

  detinfo::DetectorClocksData fClockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

  unsigned int fNMCParticles;


  static const int kNMaxMCParticles = 20000;

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
  double      fMCParticleShowerLength[kNMaxMCParticles];
  double      fMCParticleShowerOpeningAngle[kNMaxMCParticles];

  std::string fTruthLabel;
  std::string fSimLabel;
  std::string fSimChannelLabel;

  const geo::Geometry* fGeom;
};


test::mcShowerPandoraAna::mcShowerPandoraAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  fTruthLabel       = p.get<std::string>("TruthLabel");
  fSimLabel         = p.get<std::string>("SimulationLabel");
  fSimChannelLabel  = p.get<std::string>("SimChannelLabel");
}


bool test::mcShowerPandoraAna::processMCTruth(art::Event const& e) {
  // Get the MC truth
  const std::vector<art::Ptr<simb::MCTruth>> mcTruthVect = dune_ana::DUNEAnaEventUtils::GetMCTruths(e, fTruthLabel);

  // Find the neutrino and store its information.
  for (unsigned int iTruth = 0; iTruth < mcTruthVect.size(); iTruth++) {

    const art::Ptr<simb::MCTruth> truth = mcTruthVect[iTruth];

    if (truth->Origin() != simb::kBeamNeutrino)
      continue;

    fNuPdgCodeTruth = truth->GetNeutrino().Nu().PdgCode();
    fNuCCNCTruth    = truth->GetNeutrino().CCNC();
    fNuModeTruth    = truth->GetNeutrino().Mode();
    fNuETruth       = truth->GetNeutrino().Nu().E();
    fNuVertexXTruth = truth->GetNeutrino().Nu().Vx();
    fNuVertexYTruth = truth->GetNeutrino().Nu().Vy();
    fNuVertexZTruth = truth->GetNeutrino().Nu().Vz();

    break; // Only one beam nu, so don't keep looping
  }

  return true;
}

void test::mcShowerPandoraAna::getMCShowerInformation(
    art::Event const& e, TrackIDToIDEsMap &idToIDEs, const simb::MCParticle &mcParticle, const int iMc
) {

  if ((std::abs(mcParticle.PdgCode()) != 11) && (std::abs(mcParticle.PdgCode()) != 22)) {
      return;
  }

  if (idToIDEs.count(mcParticle.TrackId()) == 0) {
    return;
  }

  const TLorentzVector mcVertex(mcParticle.Position(0));
  const pandora::CartesianVector vertexPosition(mcVertex.X(), mcVertex.Y(), mcVertex.Z());

  std::vector<pandora::CartesianVector> cartesianPointVector(idToIDEs[mcParticle.TrackId()]);

  std::cout << "This shower has " << cartesianPointVector.size() << " IDEs..." << std::endl;

  try {

    const lar_content::LArShowerPCA initialLArShowerPCA(lar_content::LArPfoHelper::GetPrincipalComponents(cartesianPointVector, vertexPosition));

    const pandora::CartesianVector& centroid(initialLArShowerPCA.GetCentroid());
    const pandora::CartesianVector& primaryAxis(initialLArShowerPCA.GetPrimaryAxis());
    const pandora::CartesianVector& secondaryAxis(initialLArShowerPCA.GetSecondaryAxis());
    const pandora::CartesianVector& tertiaryAxis(initialLArShowerPCA.GetTertiaryAxis());
    const pandora::CartesianVector& eigenValues(initialLArShowerPCA.GetEigenValues());

    const pandora::CartesianVector projectedVertexPosition(centroid - (primaryAxis.GetUnitVector() * (centroid - vertexPosition).GetDotProduct(primaryAxis)));

    const float testProjection(primaryAxis.GetDotProduct(projectedVertexPosition - centroid));
    const float directionScaleFactor((testProjection > std::numeric_limits<float>::epsilon()) ? -1.f : 1.f);

    const lar_content::LArShowerPCA larShowerPCA(centroid,
        primaryAxis * directionScaleFactor,
        secondaryAxis * directionScaleFactor,
        tertiaryAxis * directionScaleFactor,
        eigenValues);

    const pandora::CartesianVector& showerLength(larShowerPCA.GetAxisLengths());

    const float length(showerLength.GetX());
    const float openingAngle(larShowerPCA.GetPrimaryLength() > 0.f ?
        std::atan(larShowerPCA.GetSecondaryLength() / larShowerPCA.GetPrimaryLength()):
        0.f);

    fMCParticleShowerLength[iMc] = length;
    fMCParticleShowerOpeningAngle[iMc] = openingAngle;

  } catch (const pandora::StatusCodeException&) {
    std::cout << "WARN: Failed to extract shower PCA..." << std::endl;
  }

  return;
}

bool test::mcShowerPandoraAna::processMCParticles(art::Event const& e) {

  // Get the MC Particles
  art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fSimLabel);

  // Get the SimChannels
  art::ValidHandle<std::vector<sim::SimChannel>> simChannels = e.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelLabel);

  if (! mcParticles.isValid())
    return false;

  if (! simChannels.isValid())
    return false;

  fNMCParticles = mcParticles->size();

  std::cout << "There is " << fNMCParticles << " MC particles..." << std::endl;

  if (fNMCParticles > kNMaxMCParticles) {
    std::cout << "WARN: Too many MC Particles: (" << fNMCParticles << ")!" << std::endl;
    std::cout << "WARN: Limiting to first " << kNMaxMCParticles << " particles!" << std::endl;
    fNMCParticles = kNMaxMCParticles;
  }

  // Pre-parse the sim-channels to build up a map to use later.
  TrackIDToIDEsMap idToIDEs;

  for (const auto &simChannel : (*simChannels)) {
    for (const auto &tdcIDEs : simChannel.TDCIDEMap()) {
      const auto &IDEs = tdcIDEs.second;

      for (const auto &IDE : IDEs) {
        int trackID = IDE.trackID;
        pandora::CartesianVector position(IDE.x, IDE.y, IDE.z);

        if (idToIDEs.count(trackID) != 0) {
          idToIDEs[trackID].push_back(position);
        } else {
          idToIDEs.insert({trackID, {position}});
        }
      }
    }
  }

  std::cout << "Size of ID -> IDE Map : " << idToIDEs.size() << std::endl;

  for (unsigned int iMc = 0; iMc < mcParticles->size(); iMc++) {

    const simb::MCParticle trueParticle = mcParticles->at(iMc);

    bool isMCPrimary(false);
    trueParticle.Process() == "primary" ? isMCPrimary = true : isMCPrimary = false;
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

    this->getMCShowerInformation(e, idToIDEs, trueParticle, iMc);
  }

  return true;
}

void test::mcShowerPandoraAna::analyze(art::Event const& e) {

  // Reset to start.
  reset();

  // Get event level information, before processing...
  const art::ServiceHandle<cheat::BackTrackerService> btServ;
  fEventID = e.id().event();
  fRunID = e.id().run();
  fSubRunID = e.id().subRun();

  std::cout << "============== EVENT ID: " << fEventID << " == RUN ID: " << fRunID << " == SUBRUN ID: " << fSubRunID << " ================" << std::endl;

  // Access the truth information
  if (e.isRealData()) {
    std::cout << "Is real data, this is an MC module!" << std::endl;
    fTree->Fill();
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

  fTree->Fill();
}

void test::mcShowerPandoraAna::beginJob() {

  // Deep clean the variables
  reset(true);

  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("pandoraOutput","Pandora Output Tree");

  // Event branches
  fTree->Branch("eventID",      &fEventID,      "eventID/i");
  fTree->Branch("runID",        &fRunID,        "runID/i");
  fTree->Branch("subrunID",     &fSubRunID,     "subrunID/i");
  fTree->Branch("nMCParticles", &fNMCParticles, "nMCParticles/i");

  // MC truth branches
  fTree->Branch("nuPdgCodeTruth", &fNuPdgCodeTruth, "nuPdgCodeTruth/i");
  fTree->Branch("nuCCNCTruth",    &fNuCCNCTruth,    "nuCCNCTruth/i");
  fTree->Branch("nuModeTruth",    &fNuModeTruth,    "nuModeTruth/i");
  fTree->Branch("nuETruth",       &fNuETruth,       "nuETruth/d");
  fTree->Branch("nuVertexXTruth", &fNuVertexXTruth, "nuVertexXTruth/d");
  fTree->Branch("nuVertexYTruth", &fNuVertexYTruth, "nuVertexYTruth/d");
  fTree->Branch("nuVertexZTruth", &fNuVertexZTruth, "nuVertexZTruth/d");

  // MC particle branches
  fTree->Branch("mcIsMCPrimary",                &fMCIsPrimary,                  "MCIsPrimary[nMCParticles]/O");
  fTree->Branch("mcParticlePdgCode",            &fMCParticlePdgCode,            "MCParticlePdgCode[nMCParticles]/I");
  fTree->Branch("mcParticleTrueEnergy",         &fMCParticleTrueEnergy,         "MCParticleTrueEnergy[nMCParticles]/D");
  fTree->Branch("mcParticleTrackID",            &fMCParticleTrackID,            "MCParticleTrackID[nMCParticles]/I");
  fTree->Branch("mcParticleParentTrackID",      &fMCParticleParentTrackID,      "MCParticleParentTrackID[nMCParticles]/I");
  fTree->Branch("mcParticleStartProcess",       &fMCParticleStartProcess,       "MCParticleStartProcess[nMCParticles]/C");
  fTree->Branch("mcParticleEndProcess",         &fMCParticleEndProcess,         "MCParticleEndProcess[nMCParticles]/C");
  fTree->Branch("mcParticleNTrajectoryPoints",  &fMCParticleNTrajectoryPoint,   "MCParticleNTrajectoryPoint[nMCParticles]/I");
  fTree->Branch("mcParticleStartPositionX",     &fMCParticleStartPositionX,     "MCParticleStartPositionX[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionY",     &fMCParticleStartPositionY,     "MCParticleStartPositionY[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionZ",     &fMCParticleStartPositionZ,     "MCParticleStartPositionZ[nMCParticles]/D");
  fTree->Branch("mcParticleStartPositionT",     &fMCParticleStartPositionT,     "MCParticleStartPositionT[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumX",     &fMCParticleStartMomentumX,     "MCParticleStartMomentumX[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumY",     &fMCParticleStartMomentumY,     "MCParticleStartMomentumY[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumZ",     &fMCParticleStartMomentumZ,     "MCParticleStartMomentumZ[nMCParticles]/D");
  fTree->Branch("mcParticleStartMomentumE",     &fMCParticleStartMomentumE,     "MCParticleStartMomentumE[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionX",       &fMCParticleEndPositionX,       "MCParticleEndPositionX[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionY",       &fMCParticleEndPositionY,       "MCParticleEndPositionY[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionZ",       &fMCParticleEndPositionZ,       "MCParticleEndPositionZ[nMCParticles]/D");
  fTree->Branch("mcParticleEndPositionT",       &fMCParticleEndPositionT,       "MCParticleEndPositionT[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumX",       &fMCParticleEndMomentumX,       "MCParticleEndMomentumX[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumY",       &fMCParticleEndMomentumY,       "MCParticleEndMomentumY[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumZ",       &fMCParticleEndMomentumZ,       "MCParticleEndMomentumZ[nMCParticles]/D");
  fTree->Branch("mcParticleEndMomentumE",       &fMCParticleEndMomentumE,       "MCParticleEndMomentumE[nMCParticles]/D");
  fTree->Branch("mcParticleVertexTime",         &fMCParticleVertexTime,         "MCParticleVertexTime[nMCParticles]/D");
  fTree->Branch("mcParticleEndTime",            &fMCParticleEndTime,            "MCParticleEndTime[nMCParticles]/D");
  fTree->Branch("mcParticleShowerLength",       &fMCParticleShowerLength,       "MCParticleShowerLength[nMCParticles]/D");
  fTree->Branch("mcParticleShowerOpeningAngle", &fMCParticleShowerOpeningAngle, "MCParticleShowerOpeningAngle[nMCParticles]/D");
}

void test::mcShowerPandoraAna::endJob()
{
  // Implementation of optional member function here.
}

void test::mcShowerPandoraAna::reset(bool deepClean) {

    fNuPdgCodeTruth = 999999;
    fNuModeTruth    = 999999;
    fNuCCNCTruth    = 999999;
    fNuETruth       = 999999;
    fNuVertexXTruth = 999999;
    fNuVertexYTruth = 999999;
    fNuVertexZTruth = 999999;

    for (unsigned int iMc = 0; iMc < (deepClean ? kNMaxMCParticles : fNMCParticles); iMc++) {
      fMCIsPrimary[iMc]                  = 0;
      fMCParticlePdgCode[iMc]            = 999999;
      fMCParticleTrueEnergy[iMc]         = 999999;
      fMCParticleTrackID[iMc]            = 999999;
      fMCParticleParentTrackID[iMc]      = 999999;
      fMCParticleStartProcess[iMc]       = "";
      fMCParticleEndProcess[iMc]         = "";
      fMCParticleNTrajectoryPoint[iMc]   = 999999;
      fMCParticleStartPositionX[iMc]     = 999999;
      fMCParticleStartPositionY[iMc]     = 999999;
      fMCParticleStartPositionZ[iMc]     = 999999;
      fMCParticleStartPositionT[iMc]     = 999999;
      fMCParticleStartMomentumX[iMc]     = 999999;
      fMCParticleStartMomentumY[iMc]     = 999999;
      fMCParticleStartMomentumZ[iMc]     = 999999;
      fMCParticleStartMomentumE[iMc]     = 999999;
      fMCParticleEndPositionX[iMc]       = 999999;
      fMCParticleEndPositionY[iMc]       = 999999;
      fMCParticleEndPositionZ[iMc]       = 999999;
      fMCParticleEndPositionT[iMc]       = 999999;
      fMCParticleEndMomentumX[iMc]       = 999999;
      fMCParticleEndMomentumY[iMc]       = 999999;
      fMCParticleEndMomentumZ[iMc]       = 999999;
      fMCParticleEndMomentumE[iMc]       = 999999;
      fMCParticleVertexTime[iMc]         = 999999;
      fMCParticleEndTime[iMc]            = 999999;
      fMCParticleShowerLength[iMc]       = 999999;
      fMCParticleShowerOpeningAngle[iMc] = 999999;
    }

    fNMCParticles = 0;

    return;
}

DEFINE_ART_MODULE(test::mcShowerPandoraAna)
