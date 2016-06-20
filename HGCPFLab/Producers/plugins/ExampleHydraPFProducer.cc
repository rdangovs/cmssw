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

    for ( std::size_t i = 0 ; i < hydraObj->simTrackSize() ; i++) {
        edm::Ptr<SimTrack> mySimTrack = hydraObj->simTrack( i );
        auto p = mySimTrack->momentum();
        float track_e = p.E();
        float sum_e = 0.;
        std::cout << " SimTrack " << i << " has pt eta geantId " << p.Pt() << " " << p.Eta() << " " << mySimTrack->trackId() << std::endl;
        std::vector<edm::Ptr<PCaloHit> > hitList = hydraObj->simHitsFromSimTrack( i );
        std::cout << "   List of " << hitList.size() << " attached hit energies: " << std::endl;
        for ( std::size_t j = 0 ; j < hitList.size() ; j++ ) {                                                                                    
            edm::Ptr<PCaloHit> aSimHit = hitList[j];                                                                                              
            std::cout << "     * SimHit has energy " << aSimHit->energy() << " detId " << aSimHit->id() << std::endl;                          
            sum_e += aSimHit->energy();
        }
        std::cout << "   Sum of energies of hits / track energy: " << sum_e << " " << track_e << std::endl;
        sum_e = 0.;
        hitList = hydraObj->simHitsFromSimTrack( i, true );
        std::cout << "   List of " << hitList.size() << " attached hit energies including all subsequent daughters: " << std::endl;
        for ( std::size_t j = 0 ; j < hitList.size() ; j++ ) {
            edm::Ptr<PCaloHit> aSimHit = hitList[j];
            std::cout << "     * SimHit has energy " << aSimHit->energy() << " detId " << aSimHit->id() << std::endl;
            sum_e += aSimHit->energy();
        }
        std::cout << "   Sum of energies of hits / track energy: " << sum_e << " " << track_e << std::endl;
    }

    for ( std::size_t i = 0 ; i < hydraObj->simHitSize() ; i++) {
        std::cout << "    SimHit " << i << " has collection index " << hydraObj->simHitCollectionIndex( i );
        edm::Ptr<PCaloHit> mySimHit = hydraObj->simHit( i );
        std::cout << " energy " << mySimHit->energy() << " detId " << mySimHit->id();
        if ( hydraObj->hasSimTrackFromSimHit( i ) ) {
            edm::Ptr<SimTrack> mySimTrack = hydraObj->simTrackFromSimHit( i );
            std::cout << " SimTrack pt " << mySimTrack->momentum().Pt() << std::endl;
            /*
            std::cout << "      CONSISTENCY: " << mySimTrack->trackId() << " " << mySimHit->geantTrackId() << std::endl;
            Index_t mySimTrackIndex = hydraObj->simTrackFromSimHitIndex( i );
            std::vector<edm::Ptr<PCaloHit> > hitList = hydraObj->simHitsFromSimTrack( mySimTrackIndex );
            std::cout << "      List of simHits (total " << hitList.size() << ") associated to this SimTrack: " << std::endl;
            for ( std::size_t j = 0 ; j < hitList.size() ; j++ ) {
                edm::Ptr<PCaloHit> aSimHit = hitList[j];
                std::cout << "        * SimHit has energy " << aSimHit->energy() << " detId " << aSimHit->id() << std::endl;
            }
            */
        } else {
            std::cout << "   ... has no SimTrack ( listed Id " << mySimHit->geantTrackId() << " )" << std::endl;
        }
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
