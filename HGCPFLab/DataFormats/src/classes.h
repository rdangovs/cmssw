#include "HGCPFLab/DataFormats/interface/HighRapidityDevRecoAssociation.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <vector>

using namespace std;

namespace {
    struct dictionary {
        HighRapidityDevRecoAssociation hydra;
        std::vector<HighRapidityDevRecoAssociation> vec_hydra;
        edm::Wrapper<std::vector<HighRapidityDevRecoAssociation> > wrp_vec_hydra;
        std::vector<edm::Ptr<reco::PFRecTrack> > trackVecPtrs;
        std::vector<edm::Ptr<reco::PFRecHit> > recHitVecPtrs;
        std::vector<edm::Ptr<reco::GenParticle> > genParticleVecPtrs;
        std::vector<edm::Ptr<SimTrack> > simTrackVecPtrs;
        std::vector<edm::Ptr<SimVertex> > simVertexVecPtrs;
        std::vector<edm::Ptr<PCaloHit> > simHitVecPtrs;
        edm::PtrVector<reco::PFRecTrack> trackPtrVec;
        edm::PtrVector<reco::PFRecHit> recHitPtrVec;
        edm::PtrVector<reco::GenParticle> genParticlePtrVec;
        edm::PtrVector<SimTrack> simTrackPtrVec;
        edm::PtrVector<SimVertex> simVertexPtrVec;
        edm::PtrVector<PCaloHit> simHitPtrVec;
    };
}


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
