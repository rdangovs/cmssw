#ifndef HGCPFLab_Hydra_h
#define HGCPFLab_Hydra_h

#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"

#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include <vector>

// Workaround for reflex not liking C++11 prior to 76X
#if !defined(__CINT__) && !defined(__MAKECINT__) && !defined(__REFLEX__)
#include <unordered_map>
#define Map std::unordered_map
#define MultiMap std::unordered_multimap
#else 
#include <map>
#define Map std::map
#define MultiMap std::multimap
#endif

// Helpful labelling
typedef unsigned Index_t;
typedef int Barcode_t;
typedef std::pair<unsigned,float> IndexAndFraction_t;
typedef std::pair<Index_t,Index_t> IndexPair_t;
typedef uint32_t RecoDetId_t;

class HighRapidityDevRecoAssociation {

public:
    HighRapidityDevRecoAssociation();
    virtual ~HighRapidityDevRecoAssociation();
    
    void insertTrack(const edm::Ptr<reco::PFRecTrack> &);
    void insertRecHit(const edm::Ptr<reco::PFRecHit> &);
    void insertGenParticle(Barcode_t, const edm::Ptr<reco::GenParticle> &);
    void insertSimTrack(Barcode_t, const edm::Ptr<SimTrack> &);
    void insertSimVertex(Barcode_t, const edm::Ptr<SimVertex> &);
    void insertSimHit(Barcode_t, const edm::Ptr<PCaloHit> &);

    std::size_t genParticleMapSize() const { return m_genParticleBarcodeToIndex.size(); }
    std::size_t genParticleSize() const { return m_genParticlePtrs.size(); }
    std::size_t genParticleBarcodeSize() const { return m_genParticleBarcodes.size(); }

    // methods (primarily for iorule) to build transient hashmaps
    // for iorule these apparently have to be public
    void buildGenParticleMap(bool clear_existing = false);
    void buildSimTrackMap(bool clear_existing = false);
    void buildSimVertexMap(bool clear_existing = false);
    void buildSimHitMap(bool clear_existing = false);
    void buildSimHitToSimTrackMap(bool clear_existing = false);
    void buildSimVertexToSimTrackMap(bool clear_existing = false);
    void buildSimTrackToSimVertexMap(bool clear_existing = false);
    void buildSimVertexToSimTrackParentMap(bool clear_existing = false);
    void buildRecoDetIdToSimHitMap(bool clear_existing = false);

private:

    // Hashmaps: not persistent
    Map<Barcode_t,Index_t> m_genParticleBarcodeToIndex;
    Map<Barcode_t,Index_t> m_simTrackBarcodeToIndex;
    Map<Barcode_t,Index_t> m_simVertexBarcodeToIndex;
    Map<Barcode_t,Index_t> m_simHitBarcodeToIndex;
    MultiMap<Index_t,Index_t> m_simHitsToSimTracks; 
    MultiMap<Index_t,Index_t> m_simVertexToSimTracks;
    MultiMap<Index_t,Index_t> m_simTrackToSimVertex; 
    MultiMap<Index_t,Index_t> m_simVertexToSimTrackParent; 
    MultiMap<RecoDetId_t,IndexAndFraction_t> m_recoDetIdToSimHits;

    // Std::Vectors: persistent, potentially slow at runtime
    // These are used to build the Hashmaps above
    // Other methods will not use these directly
    std::vector<Barcode_t> m_genParticleBarcodes; // indexed as m_genParticlePtrs
    std::vector<Barcode_t> m_simTrackBarcodes; // indexed as m_simTrackPtrs
    std::vector<Barcode_t> m_simVertexBarcodes; // indexed as m_simVertexPtrs
    std::vector<Barcode_t> m_simHitBarcodes; // indexed as m_simHitPtrs
    std::vector<IndexPair_t> m_simHitsAndSimTracks;
    std::vector<IndexPair_t> m_simVertexAndSimTracks;
    std::vector<IndexPair_t> m_simTrackAndSimVertex;
    std::vector<IndexPair_t> m_simVertexAndSimTrackParent;
    std::vector<std::pair<RecoDetId_t,IndexAndFraction_t> > m_recoDetIdAndSimHits;
    
    // Object ptrs: persistent and used in runtime methods
    edm::PtrVector<reco::PFRecTrack> m_trackPtrs;
    edm::PtrVector<reco::PFRecHit> m_recHitPtrs;
    edm::PtrVector<reco::GenParticle> m_genParticlePtrs;
    edm::PtrVector<SimTrack> m_simTrackPtrs;
    edm::PtrVector<SimVertex> m_simVertexPtrs;
    edm::PtrVector<PCaloHit> m_simHitPtrs;
};

typedef HighRapidityDevRecoAssociation Hydra;

#endif

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
