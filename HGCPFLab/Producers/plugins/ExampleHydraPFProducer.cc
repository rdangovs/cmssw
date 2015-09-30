// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HGCPFLab/DataFormats/interface/HighRapidityDevRecoAssociation.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace std;
using namespace edm;

class ExampleHydraPFProducer : public EDProducer
{
public:
    ExampleHydraPFProducer( const ParameterSet & );
    
private:
    void produce( Event &, const EventSetup & ) override;

    EDGetTokenT<View<Hydra> > tokenHydra_;
    unique_ptr<const Hydra> hydraObj;
};

ExampleHydraPFProducer::ExampleHydraPFProducer( const ParameterSet &iConfig ) :
    tokenHydra_( consumes<View<Hydra> >( iConfig.getParameter<InputTag> ( "HydraTag" ) ) )
{
    // produces
}

void ExampleHydraPFProducer::produce( Event &iEvent, const EventSetup & )
{
    Handle<View<Hydra> > HydraHandle;
    iEvent.getByToken(tokenHydra_, HydraHandle);
    assert ( HydraHandle->size() == 1 );
    hydraObj.reset( HydraHandle->ptrAt(0).get() );

    std::cout << "ExampleHydraPFProducer::produce size testing: " << std::endl;
    std::cout << "  genParticleSize=" << hydraObj->genParticleSize() << std::endl;
    std::cout << "  genParticleBarcodeSize=" << hydraObj->genParticleBarcodeSize() << std::endl;
    std::cout << "  genParticleMapSize=" << hydraObj->genParticleMapSize() << std::endl;
}

DEFINE_FWK_MODULE( ExampleHydraPFProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
