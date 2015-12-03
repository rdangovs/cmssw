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
    //    buildSimHitToSimTrackMap( *h );
    buildSimVertexToSimTrackMap( *h );
    buildSimTrackToSimVertexMap( *h );
    buildSimVertexToSimTrackParentMap( *h );
    buildRecoDetIdToSimHitMap( *h );
}

std::size_t HydraWrapper::recHitSize() const {
    std::size_t result = 0;
    for ( unsigned i = 0 ; i < 3 ; i++ ) {
        result += m_hydraCore->m_recHitPtrs[i].size();
    }
    return result;
}

std::size_t HydraWrapper::recHitCollectionIndex( std::size_t hitIndex ) const {
    std::size_t running_total = 0;
    for ( unsigned collIndex = 0 ; collIndex < 3; collIndex++ ) {
        running_total += m_hydraCore->m_recHitPtrs[collIndex].size();
        if ( hitIndex < running_total ) {
            return collIndex;
        }
    }
    throw cms::Exception( "OutOfBounds" ) << " Requested particle is larger than the number of recHits";
}

edm::Ptr<reco::PFRecHit> HydraWrapper::recHit( std::size_t hitIndex ) const {
    IndexPair_t indices = recHitInternalIndices( hitIndex );
    return m_hydraCore->m_recHitPtrs[indices.first][indices.second];
}

IndexPair_t HydraWrapper::recHitInternalIndices( std::size_t hitIndex ) const {
    std::size_t hitIndexInColl = hitIndex;
    for ( unsigned collIndex = 0 ; collIndex < 3; collIndex++ ) {
        if ( hitIndexInColl < m_hydraCore->m_recHitPtrs[collIndex].size() ) {
            return std::make_pair(collIndex, hitIndexInColl);
        } else {
            hitIndexInColl -= m_hydraCore->m_recHitPtrs[collIndex].size();
        }
    }
    throw cms::Exception( "OutOfBounds" ) << " Requested recHit is larger than the number of recHits";
}

edm::Ref<reco::PFRecHitCollection> HydraWrapper::recHitRef( std::size_t i ) const {
    edm::Ptr<reco::PFRecHit> recHitPtr = recHit( i );
    IndexPair_t indices = recHitInternalIndices( i );
    assert(!recHitPtr.isTransient());
    edm::EDProductGetter const* getter = m_hydraCore->m_recHitPtrs[indices.first].productGetter();
    assert(getter);
    return edm::Ref<reco::PFRecHitCollection>(recHitPtr.id(), recHitPtr.key(), getter);
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
    IndexPair_t indices = simHitInternalIndices( hitIndex );
    return m_hydraCore->m_simHitPtrs[indices.first][indices.second];
}

Index_t HydraWrapper::simHitExternalIndex( std::size_t collIndex, std::size_t hitIndexInColl ) const {
    Index_t result = hitIndexInColl;
    for ( unsigned collIndex = 0 ; collIndex < collIndex; collIndex++ ) {
        result += m_hydraCore->m_simHitPtrs[collIndex].size();
    }
    return result;
}

IndexPair_t HydraWrapper::simHitInternalIndices( std::size_t hitIndex ) const {
    std::size_t hitIndexInColl = hitIndex;
    for ( unsigned collIndex = 0 ; collIndex < 3; collIndex++ ) {
        if ( hitIndexInColl < m_hydraCore->m_simHitPtrs[collIndex].size() ) {
            return std::make_pair(collIndex, hitIndexInColl);
        } else {
            hitIndexInColl -= m_hydraCore->m_simHitPtrs[collIndex].size();
        }
    }
    throw cms::Exception( "OutOfBounds" ) << " Requested simHit is larger than the number of simHits";
}

bool HydraWrapper::hasSimTrackFromSimHit( std::size_t i ) const {
    Barcode_t trackId = simHit( i )->geantTrackId();
    return m_simTrackBarcodeToIndex.count( trackId );
}

Index_t HydraWrapper::simTrackFromSimHitIndex( std::size_t i ) const {
    if ( ! hasSimTrackFromSimHit( i ) ) {
        throw cms::Exception( "NoSimTrack" ) << " Requested simHit lists a barcode for which a simtrack is missing - pileup?";
    }
    return m_simTrackBarcodeToIndex.at(simHit( i )->geantTrackId());
}


edm::Ptr<SimTrack> HydraWrapper::simTrackFromSimHit( std::size_t i ) const { 
    return m_hydraCore->m_simTrackPtrs[simTrackFromSimHitIndex( i )];
}

std::vector<edm::Ptr<PCaloHit> > HydraWrapper::simHitsFromSimTrack( std::size_t i, bool include_indirect_descendants ) const {
    std::vector<edm::Ptr<PCaloHit> > result;
    Barcode_t barcode = simTrack( i )->trackId();
    auto range = m_simHitBarcodeToIndices.equal_range( barcode );
    for ( auto iter = range.first ; iter != range.second ; iter++ ) {
        result.push_back ( m_hydraCore->m_simHitPtrs[iter->second.first][iter->second.second] );
    }
    if (include_indirect_descendants) {
        auto daughter_indices = daughterSimTrackIndicesFromSimTrack( i );
        for ( auto dau_i : daughter_indices ) {
            auto daughter_result = simHitsFromSimTrack( dau_i, true );
            result.insert( result.end(), daughter_result.begin(), daughter_result.end() );
        }
    }
    return result;
}

std::vector<Index_t> HydraWrapper::daughterSimTrackIndicesFromSimTrack ( std::size_t i ) const {
    std::vector<Index_t> result;
    Barcode_t trackBarcode = simTrack( i )->trackId();
    if ( m_simTrackToSimVertex.count(trackBarcode) ) {
        auto vertex_range = m_simTrackToSimVertex.equal_range( trackBarcode );
        for ( auto vertex_iter = vertex_range.first ; vertex_iter != vertex_range.second ; vertex_iter++ ) {
            Index_t decayVertexIndex = vertex_iter->second;
            Barcode_t decayVertexBarcode = m_hydraCore->m_simVertexBarcodes[decayVertexIndex];
            auto track_range = m_simVertexToSimTracks.equal_range( decayVertexBarcode );
            for ( auto track_iter = track_range.first ; track_iter != track_range.second ; track_iter++ ) {
                result.push_back ( track_iter->second );
            }
        }
    }
    return result;
}

std::vector<edm::Ptr<SimTrack> > HydraWrapper::daughterSimTracksFromSimTrack ( std::size_t i ) const {
    std::vector<edm::Ptr<SimTrack> > result;
    std::vector<Index_t> indices = daughterSimTrackIndicesFromSimTrack(i);
    for ( auto index : indices ) {
        result.push_back ( m_hydraCore->m_simTrackPtrs[index] );
    }
    return result;
}

bool HydraWrapper::hasSimHitFromRecHit( std::size_t i ) const {
    RecoDetId_t currentRecoDetId = recHit( i )->detId();
    return m_recoDetIdToSimHits.count( currentRecoDetId );
}

std::vector<std::pair<Index_t,float> > HydraWrapper::simHitIndexesAndFractionsFromRecHit( std::size_t i ) const {
    std::vector<std::pair<Index_t,float> > result;
    RecoDetId_t currentRecoDetId = recHit( i )->detId();
    auto range = m_recoDetIdToSimHits.equal_range( currentRecoDetId );
    float ftot = 0.;
    for ( auto iter = range.first ; iter != range.second ; iter++ ) {
        Index_t result_i = simHitExternalIndex( iter->second.first.first, iter->second.first.second );
        float result_f = iter->second.second;
        ftot += result_f;
        result.push_back( std::make_pair(result_i, result_f) );
    }
    std::cout << " quick debug: f_tot=" << ftot << std::endl;
    return result;
}

std::vector<std::pair<edm::Ptr<PCaloHit>,float > > HydraWrapper::simHitsAndFractionsFromRecHit( std::size_t i ) const {
    std::vector<std::pair<edm::Ptr<PCaloHit>,float > > result;
    for ( auto ind_and_fraction : simHitIndexesAndFractionsFromRecHit( i ) ) {
        result.push_back( std::make_pair( simHit( ind_and_fraction.first ), ind_and_fraction.second ) );
    }
    return result;
}

std::pair<Index_t,float> HydraWrapper::simHitIndexAndFractionFromRecHit( std::size_t i ) const {
    //     std::vector<std::pair<RecoDetId_t,SimHitInfo_t> > m_recoDetIdAndSimHits;
    //        typedef std::pair<IndexPair_t,float> SimHitInfo_t;                                                                                  
    RecoDetId_t currentRecoDetId = recHit( i )->detId();
    Index_t best_i = std::numeric_limits<Index_t>::max();
    float best_f = -1.;
    auto range = m_recoDetIdToSimHits.equal_range( currentRecoDetId );
    for ( auto iter = range.first ; iter != range.second ; iter++ ) {
        if ( iter->second.second > best_f ) {
            best_f = iter->second.second;
            best_i = simHitExternalIndex( iter->second.first.first, iter->second.first.second );
        }
    }
    return std::make_pair( best_i, best_f );
}

std::pair<edm::Ptr<PCaloHit>,float> HydraWrapper::simHitAndFractionFromRecHit( std::size_t i ) const {
    auto ind_and_fraction = simHitIndexAndFractionFromRecHit( i );
    return std::make_pair( simHit( ind_and_fraction.first ), ind_and_fraction.second );
}


/*
void HighRapidityDevRecoAssociation::insertSimVertex(const edm::Ptr<SimVertex> & sv)
{   
    Index_t i = m_simVertexPtrs.size();
    m_simVertexPtrs.push_back(sv);
    m_simVertexBarcodes.push_back(sv->vertexId());
    if (!sv->noParent()) {
        m_simTrackAndSimVertex.push_back( std::make_pair(sv->parentIndex(), i ) );
        m_simVertexAndSimTrackParent.push_back( std::make_pair( sv->vertexId(), sv->parentIndex() ) );
    }
}
*/

void HydraWrapper::buildGenParticleMap(const Hydra &h, bool clear_existing) {
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

/*
void HydraWrapper::buildSimHitToSimTrackMap(const Hydra &h, bool clear_existing) {
    if (m_simHitsToSimTracks.size() > 0) {
        if (clear_existing) {
            m_simHitsToSimTracks.clear();
        } else {
            throw cms::Exception( "NotImplemented" ) << "Building an index map when one is already built not currently supported";
        }
    }
    for (auto it = h.m_simHitsAndSimTracks.begin() ; it != h.m_simHitsAndSimTracks.end() ; it++) {
        m_simHitsToSimTracks.emplace(it->first,it->second);
        m_simTracksToSimHits.emplace(it->second,it->first);
    }
}
*/

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

