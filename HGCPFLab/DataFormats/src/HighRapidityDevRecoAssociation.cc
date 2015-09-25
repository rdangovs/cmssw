#include "HGCPFLab/DataFormats/interface/HighRapidityDevRecoAssociation.h"
#include "FWCore/Utilities/interface/Exception.h"

using namespace reco;

HighRapidityDevRecoAssociation::HighRapidityDevRecoAssociation()
{}

HighRapidityDevRecoAssociation::~HighRapidityDevRecoAssociation()
{}

void HighRapidityDevRecoAssociation::insertTrack(const TrackRef & trk) {
    m_trackRefs.push_back(trk);
}

void HighRapidityDevRecoAssociation::insertRecHit(const PFRecHitRef & hit) {
    m_recHitRefs.push_back(hit);
}

void HighRapidityDevRecoAssociation::insertGenParticle(Barcode_t barcode, const GenParticleRef & gp) 
{
    m_genParticleRefs.push_back(gp);
    m_genParticleBarcodes.push_back(barcode);
}

void HighRapidityDevRecoAssociation::insertSimTrack(Barcode_t barcode, const SimTrackRef & st)
{
    m_simTrackRefs.push_back(st);
    m_simTrackBarcodes.push_back(barcode);
}

void HighRapidityDevRecoAssociation::insertSimVertex(Barcode_t barcode, const SimVertexRef & sv)
{
    m_simVertexRefs.push_back(sv);
    m_simVertexBarcodes.push_back(barcode);
}

void HighRapidityDevRecoAssociation::insertSimHit(Barcode_t barcode, const SimHitRef & sh)
{
    m_simHitRefs.push_back(sh);
    m_simHitBarcodes.push_back(barcode);
}

void HighRapidityDevRecoAssociation::buildGenParticleMap(bool clear_existing) {
    if (m_genParticleBarcodeToIndex.size() > 0) {
        if (clear_existing) {
            m_genParticleBarcodeToIndex.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building a barcode map when one is already built not currently supported";
        }
    }
    for (unsigned i = 0 ; i < m_genParticleBarcodes.size() ; i++) {
        m_genParticleBarcodeToIndex[m_genParticleBarcodes[i]] = i;
    }
}

void HighRapidityDevRecoAssociation::buildSimTrackMap(bool clear_existing) {
    if (m_simTrackBarcodeToIndex.size() > 0) {
        if (clear_existing) {
            m_simTrackBarcodeToIndex.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building a barcode map when one is already built not currently supported";
        }
    }
    for (unsigned i = 0 ; i < m_simTrackBarcodes.size() ; i++) {
        m_simTrackBarcodeToIndex[m_simTrackBarcodes[i]] = i;
    }
}

void HighRapidityDevRecoAssociation::buildSimVertexMap(bool clear_existing) {
    if (m_simVertexBarcodeToIndex.size() > 0) {
        if (clear_existing) {
            m_simVertexBarcodeToIndex.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building a barcode map when one is already built not currently supported";
        }
    }
    for (unsigned i = 0 ; i < m_simVertexBarcodes.size() ; i++) {
        m_simVertexBarcodeToIndex[m_simVertexBarcodes[i]] = i;
    }
}

void HighRapidityDevRecoAssociation::buildSimHitMap(bool clear_existing) {
    if (m_simHitBarcodeToIndex.size() > 0) {
        if (clear_existing) {
            m_simHitBarcodeToIndex.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building a barcode map when one is already built not currently supported";
        }
    }
    for (unsigned i = 0 ; i < m_simHitBarcodes.size() ; i++) {
        m_simHitBarcodeToIndex[m_simHitBarcodes[i]] = i;
    }
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

