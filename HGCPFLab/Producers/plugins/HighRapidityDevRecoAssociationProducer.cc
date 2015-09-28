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

class HighRapidityDevRecoAssociationProducer : public EDProducer
{
public:
    HighRapidityDevRecoAssociationProducer( const ParameterSet & );
    
private:
    void produce( Event &, const EventSetup & ) override;

    EDGetTokenT<View<PFRecHit> > tokenHGCrechit_;
};

HighRapidityDevRecoAssociationProducer::HighRapidityDevRecoAssociationProducer( const ParameterSet &iConfig ) :
    tokenHGCrechit_( consumes<View<PFRecHit> >( iConfig.getParameter<InputTag> ( "HGCrechitCollection" ) ) )
{
    produces<vector<HyDRA> >();
}

void HighRapidityDevRecoAssociationProducer::produce( Event &iEvent, const EventSetup & )
{
    auto_ptr<vector<HyDRA> > output( new vector<HyDRA> );
    output->emplace_back();

    Handle<View<PFRecHit> > HGCRecHitHandle;
    iEvent.getByToken(tokenHGCrechit_, HGCRecHitHandle);

    for(unsigned i=0; i<HGCRecHitHandle->size(); i++) {
        output->back().insertRecHit(HGCRecHitHandle->ptrAt(i));
    }        

    iEvent.put( output );
}

DEFINE_FWK_MODULE( HighRapidityDevRecoAssociationProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
