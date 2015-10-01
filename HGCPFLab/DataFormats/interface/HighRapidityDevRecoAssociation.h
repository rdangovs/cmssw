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

// Helpful labelling
typedef unsigned Index_t;
typedef int Barcode_t;
typedef std::pair<unsigned,float> IndexAndFraction_t;
typedef std::pair<Index_t,Index_t> IndexPair_t;
typedef uint32_t RecoDetId_t;

class HighRapidityDevRecoAssociation {
    friend class HydraWrapper;
public:
    HighRapidityDevRecoAssociation();
    virtual ~HighRapidityDevRecoAssociation();
    
    void insertTrack(const edm::Ptr<reco::PFRecTrack> &);
    void insertRecHit(const edm::Ptr<reco::PFRecHit> &);
    void insertGenParticle(Barcode_t, const edm::Ptr<reco::GenParticle> &);
    void insertSimTrack(Barcode_t, const edm::Ptr<SimTrack> &);
    void insertSimVertex(Barcode_t, const edm::Ptr<SimVertex> &);
    void insertSimHit(Barcode_t, const edm::Ptr<PCaloHit> &);

    std::size_t genParticleSize() const { return m_genParticlePtrs.size(); }
    std::size_t genParticleBarcodeSize() const { return m_genParticleBarcodes.size(); }

private:

    // Std::Vectors: persistent, potentially slow at runtime
    // These are used to build the Hashmaps above
    // Other methods will not use these directly
    std::vector<int> m_genParticleBarcodes; // indexed as m_genParticlePtrs
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
