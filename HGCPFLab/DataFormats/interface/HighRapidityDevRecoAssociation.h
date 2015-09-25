#ifndef HGCPFLab_HyDRA_h
#define HGCPFLab_HyDRA_h

#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include <vector>

// Workaround for reflex not liking C++11 prior to 76X
#if !defined(__CINT__) && !defined(__MAKECINT__) && !defined(__REFLEX__)
#include <unordered_map>
#define Map unordered_map
#define MultiMap unordered_multimap
#else 
#include <map>
#define Map map
#define MultiMap multimap
#endif

using namespace edm;
using namespace std;
using namespace reco;

// Usual CMS conventions, but these seem not to be defined anywhere else
typedef vector<SimTrack> SimTrackCollection;
typedef Ref<SimTrackCollection> SimTrackRef;
typedef RefVector<SimTrackCollection> SimTrackRefVector;
typedef vector<SimVertex> SimVertexCollection;
typedef Ref<SimVertexCollection> SimVertexRef;
typedef RefVector<SimVertexCollection> SimVertexRefVector;
typedef vector<PCaloHit> SimHitCollection;
typedef Ref<SimHitCollection> SimHitRef;
typedef RefVector<SimHitCollection> SimHitRefVector;

// Helpful labelling
typedef unsigned Index_t;
typedef unsigned Barcode_t;
typedef pair<unsigned,float> IndexAndFraction_t;
typedef pair<Index_t,Index_t> IndexPair_t;
typedef uint32_t RecoDetId_t;

class HighRapidityDevRecoAssociation {

public:
    HighRapidityDevRecoAssociation();
    virtual ~HighRapidityDevRecoAssociation();
    
    void insertTrack(const PFRecTrackRef &);
    void insertRecHit(const PFRecHitRef &);
    void insertGenParticle(Barcode_t, const GenParticleRef &);
    void insertSimTrack(Barcode_t, const SimTrackRef &);
    void insertSimVertex(Barcode_t, const SimVertexRef &);
    void insertSimHit(Barcode_t, const SimHitRef &);

    // methods (primarily for iorule) to build transient maps
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
    MultiMap<RecoDetId_t,IndexAndFraction_t> m_recoDetIdToSimHit;

    // Vectors: persistent, potentially slow at runtime
    // These are used to build the Hashmaps above
    // Other methods will not use these directly
    vector<Barcode_t> m_genParticleBarcodes; // indexed as m_genParticleRefs
    vector<Barcode_t> m_simTrackBarcodes; // indexed as m_simTrackRefs
    vector<Barcode_t> m_simVertexBarcodes; // indexed as m_simVertexRefs
    vector<Barcode_t> m_simHitBarcodes; // indexed as m_simHitRefs
    vector<IndexPair_t> m_simHitsAndSimTracks;
    vector<IndexPair_t> m_simVertexAndSimTracks;
    vector<IndexPair_t> m_simTrackAndSimVertex;
    vector<IndexPair_t> m_simVertexAndSimTrackParent;
    vector<pair<RecoDetId_t,IndexAndFraction_t> > m_recoDetIdAndSimHit;
    
    // Object refs: persistent and used in runtime methods
    PFRecTrackRefVector m_trackRefs;
    PFRecHitRefVector m_recHitRefs;
    GenParticleRefVector m_genParticleRefs;
    SimTrackRefVector m_simTrackRefs;
    SimVertexRefVector m_simVertexRefs;
    SimHitRefVector m_simHitRefs;
};

#endif

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
