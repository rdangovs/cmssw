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

class HydraProducer : public EDProducer
{
public:
    HydraProducer( const ParameterSet & );
    
private:
    void produce( Event &, const EventSetup & ) override;

    EDGetTokenT<View<PFRecHit> > tokenHGCrechit_;
    EDGetTokenT<View<GenParticle> > tokenGenParticle_;
    EDGetTokenT<View<Barcode_t> > tokenGenBarcode_;
};

HydraProducer::HydraProducer( const ParameterSet &iConfig ) :
    tokenHGCrechit_( consumes<View<PFRecHit> >( iConfig.getParameter<InputTag> ( "HGCrechitCollection" ) ) ),
    tokenGenParticle_( consumes<View<GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleCollection" ) ) ),
    tokenGenBarcode_( consumes<View<Barcode_t> >( iConfig.getParameter<InputTag> ( "GenParticleCollection" ) ) )
{
    produces<vector<Hydra> >();
}

void HydraProducer::produce( Event &iEvent, const EventSetup & )
{
    auto_ptr<vector<Hydra> > output( new vector<Hydra> );
    output->emplace_back();

    Handle<View<PFRecHit> > HGCRecHitHandle;
    iEvent.getByToken(tokenHGCrechit_, HGCRecHitHandle);

    Handle<View<GenParticle> > GenParticleHandle;
    iEvent.getByToken(tokenGenParticle_, GenParticleHandle);
    Handle<View<Barcode_t> > GenBarcodeHandle;
    iEvent.getByToken(tokenGenBarcode_, GenBarcodeHandle);

    for(unsigned i=0; i<HGCRecHitHandle->size(); i++) {
        output->back().insertRecHit(HGCRecHitHandle->ptrAt(i));
    }        

    for(unsigned i=0; i<GenParticleHandle->size(); i++) {
        output->back().insertGenParticle(GenBarcodeHandle->at(i),GenParticleHandle->ptrAt(i));
    }

    iEvent.put( output );
}

DEFINE_FWK_MODULE( HydraProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
