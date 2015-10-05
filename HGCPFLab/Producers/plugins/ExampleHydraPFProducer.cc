#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HGCPFLab/DataFormats/interface/HighRapidityDevRecoAssociation.h"
#include "HGCPFLab/DataFormats/interface/HydraWrapper.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"

using namespace std;
using namespace edm;

class ExampleHydraPFProducer : public EDProducer
{
public:
    ExampleHydraPFProducer( const ParameterSet & );
    
private:
    void produce( Event &, const EventSetup & ) override;

    EDGetTokenT<View<Hydra> > tokenHydra_;
    unique_ptr<HydraWrapper> hydraObj;
};

ExampleHydraPFProducer::ExampleHydraPFProducer( const ParameterSet &iConfig ) :
    tokenHydra_( consumes<View<Hydra> >( iConfig.getParameter<InputTag> ( "HydraTag" ) ) )
{
    produces<reco::PFClusterCollection>();
    produces<reco::PFBlockCollection>();
    produces<reco::PFCandidateCollection>();
}

void ExampleHydraPFProducer::produce( Event &iEvent, const EventSetup & )
{
    std::cout << " ExampleHydraPFProducer::produce " << std::endl;
    Handle<View<Hydra> > HydraHandle;
    iEvent.getByToken(tokenHydra_, HydraHandle);
    assert ( HydraHandle->size() == 1 );
    hydraObj.reset( new HydraWrapper( HydraHandle->ptrAt(0)) );

    auto_ptr<reco::PFClusterCollection> pfClusterCol ( new reco::PFClusterCollection );
    auto_ptr<reco::PFBlockCollection> pfBlockCol ( new reco::PFBlockCollection );
    auto_ptr<reco::PFCandidateCollection> pfCandidateCol ( new reco::PFCandidateCollection );

    std::cout << "ExampleHydraPFProducer::produce size testing: " << std::endl;
    std::cout << "  genParticleSize=" << hydraObj->genParticleSize() << std::endl;
    std::cout << "  recTrackSize=" << hydraObj->recTrackSize() << std::endl;
    std::cout << "  recHitSize=" << hydraObj->recHitSize() << std::endl;
    std::cout << "  simTrackSize=" << hydraObj->simTrackSize() << std::endl;
    std::cout << "  simVertexSize=" << hydraObj->simVertexSize() << std::endl;
    std::cout << "  simHitSize=" << hydraObj->simHitSize() << std::endl;

    for ( std::size_t i = 0 ; i < hydraObj->simHitSize() ; i++) {
        std::cout << "    SimHit " << i << " has collection index " << hydraObj->simHitCollectionIndex( i );
        edm::Ptr<PCaloHit> mySimHit = hydraObj->simHit( i );
        std::cout << " energy " << mySimHit->energy() << " detId " << mySimHit->id() << std::endl;
    }

    iEvent.put ( pfClusterCol );
    iEvent.put ( pfBlockCol );
    iEvent.put ( pfCandidateCol );
}

DEFINE_FWK_MODULE( ExampleHydraPFProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
