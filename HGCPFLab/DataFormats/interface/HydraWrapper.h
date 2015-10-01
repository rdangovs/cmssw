#ifndef HGCPFLab_HydraWrapper_h
#define HGCPFLab_HydraWrapper_h

#include "HGCPFLab/DataFormats/interface/HighRapidityDevRecoAssociation.h"

#include <unordered_map>
#define Map std::unordered_map
#define MultiMap std::unordered_multimap

class HydraWrapper {

public:
    HydraWrapper();
    HydraWrapper( const Hydra& );
    virtual ~HydraWrapper();

    std::size_t genParticleMapSize() const { return m_genParticleBarcodeToIndex.size(); }

    void buildGenParticleMap(const Hydra &, bool clear_existing = false);
    void buildSimTrackMap(const Hydra &, bool clear_existing = false);
    void buildSimVertexMap(const Hydra &, bool clear_existing = false);
    void buildSimHitMap(const Hydra &, bool clear_existing = false);
    void buildSimHitToSimTrackMap(const Hydra &, bool clear_existing = false);
    void buildSimVertexToSimTrackMap(const Hydra &, bool clear_existing = false);
    void buildSimTrackToSimVertexMap(const Hydra &, bool clear_existing = false);
    void buildSimVertexToSimTrackParentMap(const Hydra &, bool clear_existing = false);
    void buildRecoDetIdToSimHitMap(const Hydra &, bool clear_existing = false);

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
};

#endif

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
