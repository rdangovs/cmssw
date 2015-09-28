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

class ExampleHyDRAPFProducer : public EDProducer
{
public:
    ExampleHyDRAPFProducer( const ParameterSet & );
    
private:
    void produce( Event &, const EventSetup & ) override;

    EDGetTokenT<View<HyDRA> > tokenHyDRA_;
    unique_ptr<const HyDRA> HyDRAObj;
};

ExampleHyDRAPFProducer::ExampleHyDRAPFProducer( const ParameterSet &iConfig ) :
    tokenHyDRA_( consumes<View<HyDRA> >( iConfig.getParameter<InputTag> ( "HyDRATag" ) ) )
{
    // produces
}

void ExampleHyDRAPFProducer::produce( Event &iEvent, const EventSetup & )
{
    Handle<View<HyDRA> > HyDRAHandle;
    iEvent.getByToken(tokenHyDRA_, HyDRAHandle);
    assert ( HyDRAHandle->size() == 1 );
    HyDRAObj.reset( HyDRAHandle->ptrAt(0).get() );

    std::cout << "ExampleHyDRAPFProducer::produce size testing: " << std::endl;
    std::cout << "  genParticleSize=" << HyDRAObj->genParticleSize() << std::endl;
    std::cout << "  genParticleBarcodeSize=" << HyDRAObj->genParticleBarcodeSize() << std::endl;
    std::cout << "  genParticleMapSize=" << HyDRAObj->genParticleMapSize() << std::endl;
}

DEFINE_FWK_MODULE( ExampleHyDRAPFProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
