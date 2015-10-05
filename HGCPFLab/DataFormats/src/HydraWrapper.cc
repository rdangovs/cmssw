#include "HGCPFLab/DataFormats/interface/HydraWrapper.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Utilities/interface/Exception.h"

HydraWrapper::HydraWrapper()
{}

HydraWrapper::~HydraWrapper()
{}

HydraWrapper::HydraWrapper( const edm::Ptr<Hydra> & h ) {

    m_hydraCore = h;

    buildGenParticleMap( *h );
    buildSimTrackMap( *h );
    buildSimTrackGenMap( *h );
    buildSimVertexMap( *h );
    buildSimHitMap( *h );
    buildSimHitToSimTrackMap( *h );
    buildSimVertexToSimTrackMap( *h );
    buildSimTrackToSimVertexMap( *h );
    buildSimVertexToSimTrackParentMap( *h );
    buildRecoDetIdToSimHitMap( *h );
}

std::size_t HydraWrapper::simHitSize() const {
    std::size_t result = 0;
    for ( unsigned i = 0 ; i < 3 ; i++ ) {
        result += m_hydraCore->m_simHitPtrs[i].size();
    }
    return result;
}

std::size_t HydraWrapper::simHitCollectionIndex( std::size_t hitIndex ) const {
    std::size_t running_total = 0;
    for ( unsigned collIndex = 0 ; collIndex < 3; collIndex++ ) {
        running_total += m_hydraCore->m_simHitPtrs[collIndex].size();
        if ( hitIndex < running_total ) {
            return collIndex;
        }
    }
    throw cms::Exception( "OutOfBounds" ) << " Requested particle is larger than the number of simHits";
}

edm::Ptr<PCaloHit> HydraWrapper::simHit( std::size_t hitIndex ) const {
    std::size_t hitIndexInColl = hitIndex;
    for ( unsigned collIndex = 0 ; collIndex < 3; collIndex++ ) {
        if ( hitIndexInColl < m_hydraCore->m_simHitPtrs[collIndex].size() ) {
            return m_hydraCore->m_simHitPtrs[collIndex][hitIndexInColl];
        } else {
            hitIndexInColl -= m_hydraCore->m_simHitPtrs[collIndex].size();
        }
    }
    throw cms::Exception( "OutOfBounds" ) << " Requested particle is larger than the number of simHits";
}

void HydraWrapper::buildGenParticleMap(const Hydra &h, bool clear_existing) {
    std::cout << "HydraWrapper::buildGenParticleMap" << std::endl;
    if (m_genParticleBarcodeToIndex.size() > 0) {
        if (clear_existing) {
            m_genParticleBarcodeToIndex.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building a barcode map when one is already built not currently supported";
        }
    }
    for (unsigned i = 0 ; i < h.m_genParticleBarcodes.size() ; i++) {
        m_genParticleBarcodeToIndex.emplace(h.m_genParticleBarcodes[i],i);;
    }
}

void HydraWrapper::buildSimTrackMap(const Hydra &h, bool clear_existing) {
    if (m_simTrackBarcodeToIndex.size() > 0) {
        if (clear_existing) {
            m_simTrackBarcodeToIndex.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building a barcode map when one is already built not currently supported";
        }
    }
    for (unsigned i = 0 ; i < h.m_simTrackBarcodes.size() ; i++) {
        m_simTrackBarcodeToIndex.emplace(h.m_simTrackBarcodes[i],i);
    }
}

void HydraWrapper::buildSimTrackGenMap(const Hydra &h, bool clear_existing) {
    if (m_genBarcodeToSimTrackIndex.size() > 0) {
        if (clear_existing) {
            m_genBarcodeToSimTrackIndex.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building a barcode map when one is already built not currently supported";
        }
    }
    for (unsigned i = 0 ; i < h.m_simTrackBarcodes.size() ; i++) {
        m_genBarcodeToSimTrackIndex.emplace(h.m_simTrackGenBarcodes[i],i);
    }
}

void HydraWrapper::buildSimVertexMap(const Hydra &h, bool clear_existing) {
    if (m_simVertexBarcodeToIndex.size() > 0) {
        if (clear_existing) {
            m_simVertexBarcodeToIndex.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building a barcode map when one is already built not currently supported";
        }
    }
    for (unsigned i = 0 ; i < h.m_simVertexBarcodes.size() ; i++) {
        m_simVertexBarcodeToIndex.emplace(h.m_simVertexBarcodes[i],i);;
    }
}

void HydraWrapper::buildSimHitMap(const Hydra &h, bool clear_existing) {
    if (m_simHitBarcodeToIndices.size() > 0) {
        if (clear_existing) {
            m_simHitBarcodeToIndices.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building a barcode map when one is already built not currently supported";
        }
    }
    for (unsigned iCol = 0 ; iCol < 3 ; iCol++) {
        for (unsigned iHit = 0 ; iHit < h.m_simHitBarcodes[iCol].size() ; iHit++) {
            m_simHitBarcodeToIndices.emplace(h.m_simHitBarcodes[iCol][iHit],std::make_pair(iCol,iHit));
        }
    }
}

void HydraWrapper::buildSimHitToSimTrackMap(const Hydra &h, bool clear_existing) {
    if (m_simHitsToSimTracks.size() > 0) {
        if (clear_existing) {
            m_simHitsToSimTracks.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building a barcode map when one is already built not currently supported";
        }
    }
    for (auto it = h.m_simHitsAndSimTracks.begin() ; it != h.m_simHitsAndSimTracks.end() ; it++) {
        m_simHitsToSimTracks.emplace(it->first,it->second);
    }
}

void HydraWrapper::buildSimVertexToSimTrackMap(const Hydra &h, bool clear_existing) {
    if (m_simVertexToSimTracks.size() > 0) {
        if (clear_existing) {
            m_simVertexToSimTracks.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building a barcode map when one is already built not currently supported";
        }
    }
    for(auto it = h.m_simVertexAndSimTracks.begin() ; it != h.m_simVertexAndSimTracks.end() ; it++) {
        m_simVertexToSimTracks.emplace(it->first,it->second);
    }
}

void HydraWrapper::buildSimTrackToSimVertexMap(const Hydra &h, bool clear_existing) {
    if (m_simTrackToSimVertex.size() > 0) {
        if (clear_existing) {
            m_simTrackToSimVertex.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building a barcode map when one is already built not currently supported";
        }
    }
    for(auto it = h.m_simTrackAndSimVertex.begin() ; it != h.m_simTrackAndSimVertex.end() ; it++) {
        m_simTrackToSimVertex.emplace(it->first,it->second);
    }
}

void HydraWrapper::buildSimVertexToSimTrackParentMap(const Hydra &h, bool clear_existing) {
    if (m_simVertexToSimTrackParent.size() > 0) {
        if (clear_existing) {
            m_simVertexToSimTrackParent.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building a barcode map when one is already built not currently supported";
        }
    }
    for(auto it = h.m_simVertexAndSimTrackParent.begin() ; it != h.m_simVertexAndSimTrackParent.end() ; it++) {
        m_simVertexToSimTrackParent.emplace(it->first,it->second);
    }
}

void HydraWrapper::buildRecoDetIdToSimHitMap(const Hydra &h, bool clear_existing) {
    if (m_recoDetIdToSimHits.size() > 0) {
        if (clear_existing) {
            m_recoDetIdToSimHits.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building a barcode map when one is already built not currently supported";
        }
    }
    for(auto it = h.m_recoDetIdAndSimHits.begin() ; it != h.m_recoDetIdAndSimHits.end() ; it++) {
        m_recoDetIdToSimHits.emplace(it->first,it->second);
    }
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

