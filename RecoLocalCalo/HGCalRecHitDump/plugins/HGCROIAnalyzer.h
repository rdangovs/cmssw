#ifndef _HGCROIAnalyzer_h_
#define _HGCROIAnalyzer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "RecoLocalCalo/HGCalRecHitDump/interface/SlimmedRecHit.h"
#include "RecoLocalCalo/HGCalRecHitDump/interface/SlimmedROI.h"
#include "RecoLocalCalo/HGCalRecHitDump/interface/SlimmedVertex.h"
#include "RecoLocalCalo/HGCalRecHitDump/interface/SlimmedCluster.h"

#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include <unordered_map>

/**
   @class HGCROIAnalyzer
   @author P. Silva (CERN)
*/

class HGCROIAnalyzer : public edm::EDAnalyzer 
{  
 public:
  
  explicit HGCROIAnalyzer( const edm::ParameterSet& );
  ~HGCROIAnalyzer();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );

 private:
  
  void slimRecHits(const edm::Event &iEvent, const edm::EventSetup &iSetup);
  void doMCJetMatching(edm::Handle<std::vector<reco::PFJet> > &pfJets,
		       edm::Handle<reco::GenJetCollection> &genJets,
		       edm::Handle<edm::View<reco::Candidate> > &genParticles,
		       std::unordered_map<uint32_t,uint32_t> &reco2genJet,
		       std::unordered_map<uint32_t,uint32_t> &genJet2Parton,
		       std::unordered_map<uint32_t,uint32_t> &genJet2Stable);
  void doMCJetMatching(edm::Handle<reco::SuperClusterCollection> &superClusters,
		       edm::Handle<reco::GenJetCollection> &genJets,
		       edm::Handle<edm::View<reco::Candidate> > &genParticles,
		       std::unordered_map<uint32_t,uint32_t> &reco2genJet,
		       std::unordered_map<uint32_t,uint32_t> &genJet2Parton,
		       std::unordered_map<uint32_t,uint32_t> &genJet2Stable);

  virtual void endJob() ;

  TTree *tree_;
  Int_t run_,event_,lumi_;
  std::vector<SlimmedRecHit> *slimmedRecHits_;
  std::vector<SlimmedCluster> *slimmedClusters_;
  std::vector<SlimmedROI> *slimmedROIs_;
  std::vector<SlimmedVertex> *slimmedVertices_;
  TLorentzVector *genVertex_;
  
  bool useSuperClustersAsROIs_,useStatus3ForGenVertex_;
  std::string eeSimHitsSource_, hefSimHitsSource_;
  std::string eeRecHitsSource_, hefRecHitsSource_;
  std::string g4TracksSource_, g4VerticesSource_;
  std::string recoVertexSource_;
  std::string genSource_, genCandsFromSimTracksSource_, genJetsSource_, pfJetsSource_, superClustersSource_;
};
 

#endif
