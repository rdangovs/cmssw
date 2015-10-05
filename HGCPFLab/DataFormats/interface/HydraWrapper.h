#ifndef HGCPFLab_HydraWrapper_h
#define HGCPFLab_HydraWrapper_h

#include "HGCPFLab/DataFormats/interface/HighRapidityDevRecoAssociation.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <unordered_map>
#define Map std::unordered_map
#define MultiMap std::unordered_multimap

class HydraWrapper {

public:
    HydraWrapper();
    HydraWrapper( const edm::Ptr<Hydra>& );
    virtual ~HydraWrapper();

    std::size_t genParticleSize() const { return m_hydraCore->m_genParticlePtrs.size(); }
    std::size_t simTrackSize() const { return m_hydraCore->m_simTrackPtrs.size(); }
    std::size_t simVertexSize() const { return m_hydraCore->m_simVertexPtrs.size(); }
    std::size_t recTrackSize() const { return m_hydraCore->m_trackPtrs.size(); }
    std::size_t recHitSize() const { return m_hydraCore->m_recHitPtrs.size(); }

    edm::Ptr<reco::GenParticle> genParticle( std::size_t i ) const { return m_hydraCore->m_genParticlePtrs[i]; }
    edm::Ptr<SimTrack> simTrack( std::size_t i ) const { return m_hydraCore->m_simTrackPtrs[i]; }
    edm::Ptr<SimVertex> simVertex( std::size_t i ) const { return m_hydraCore->m_simVertexPtrs[i]; }
    edm::Ptr<reco::PFRecTrack> recTrack( std::size_t i ) const { return m_hydraCore->m_trackPtrs[i]; }
    edm::Ptr<reco::PFRecHit> recHit( std::size_t i ) const { return m_hydraCore->m_recHitPtrs[i]; }

    std::size_t simHitSize() const; // suppresses 3 separate collections, use simHitCollectionIndex if needed
    std::size_t simHitCollectionIndex( std::size_t ) const; // cf. HydraProducer SimHitCollection vector
    edm::Ptr<PCaloHit> simHit( std::size_t ) const;

    void buildGenParticleMap(const Hydra &, bool clear_existing = false);
    void buildSimTrackMap(const Hydra &, bool clear_existing = false);
    void buildSimTrackGenMap(const Hydra &, bool clear_existing = false);
    void buildSimVertexMap(const Hydra &, bool clear_existing = false);
    void buildSimHitMap(const Hydra &, bool clear_existing = false);
    void buildSimHitToSimTrackMap(const Hydra &, bool clear_existing = false);
    void buildSimVertexToSimTrackMap(const Hydra &, bool clear_existing = false);
    void buildSimTrackToSimVertexMap(const Hydra &, bool clear_existing = false);
    void buildSimVertexToSimTrackParentMap(const Hydra &, bool clear_existing = false);
    void buildRecoDetIdToSimHitMap(const Hydra &, bool clear_existing = false);

private:

    edm::Ptr<Hydra> m_hydraCore;
   
    // Hashmaps: not persistent
    Map<Barcode_t,Index_t> m_genParticleBarcodeToIndex;
    Map<Barcode_t,Index_t> m_simTrackBarcodeToIndex;
    Map<Barcode_t,Index_t> m_genBarcodeToSimTrackIndex;
    Map<Barcode_t,Index_t> m_simVertexBarcodeToIndex;
    Map<Barcode_t,IndexPair_t> m_simHitBarcodeToIndices;
    MultiMap<Index_t,Index_t> m_simHitsToSimTracks; 
    MultiMap<Barcode_t,Index_t> m_simVertexToSimTracks;
    MultiMap<Barcode_t,Index_t> m_simTrackToSimVertex; 
    MultiMap<Barcode_t,Index_t> m_simVertexToSimTrackParent; 
    MultiMap<RecoDetId_t,SimHitInfo_t> m_recoDetIdToSimHits;
};

#endif

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
