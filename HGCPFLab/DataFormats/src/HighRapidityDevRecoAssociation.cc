#include "HGCPFLab/DataFormats/interface/HighRapidityDevRecoAssociation.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

HighRapidityDevRecoAssociation::HighRapidityDevRecoAssociation()
{}

HighRapidityDevRecoAssociation::~HighRapidityDevRecoAssociation()
{}

void HighRapidityDevRecoAssociation::insertTrack(const edm::Ptr<reco::PFRecTrack> & trk) {
    m_trackPtrs.push_back(trk);
}

void HighRapidityDevRecoAssociation::insertRecHit(const edm::Ptr<reco::PFRecHit> & hit) {
    m_recHitPtrs.push_back(hit);
}

void HighRapidityDevRecoAssociation::insertGenParticle(Barcode_t barcode, const edm::Ptr<reco::GenParticle> & gp) 
{
    m_genParticlePtrs.push_back(gp);
    m_genParticleBarcodes.push_back(barcode);
}

void HighRapidityDevRecoAssociation::insertSimTrack(Barcode_t barcode, const edm::Ptr<SimTrack> & st)
{
    m_simTrackPtrs.push_back(st);
    m_simTrackBarcodes.push_back(barcode);
}

void HighRapidityDevRecoAssociation::insertSimVertex(Barcode_t barcode, const edm::Ptr<SimVertex> & sv)
{
    m_simVertexPtrs.push_back(sv);
    m_simVertexBarcodes.push_back(barcode);
}

void HighRapidityDevRecoAssociation::insertSimHit(Barcode_t barcode, const edm::Ptr<PCaloHit> & sh)
{
    m_simHitPtrs.push_back(sh);
    m_simHitBarcodes.push_back(barcode);
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

