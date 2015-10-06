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

class HydraFakeClusterBuilder : public EDProducer
{
public:
    HydraFakeClusterBuilder( const ParameterSet & );
    
private:
    void produce( Event &, const EventSetup & ) override;

    EDGetTokenT<View<Hydra> > tokenHydra_;
    unique_ptr<HydraWrapper> hydraObj;

    bool useGenParticles_; // if true use SimTracks corresponding to GenParticles
                           // if false use SimTracks at calorimeter face
    
    bool splitRecHits_; // if true, split the rechits in proportion to the contributing simhits
                        // if false, assign each rechit fully to the simtrack that contributed most

    bool debugPrint_;
};

HydraFakeClusterBuilder::HydraFakeClusterBuilder( const ParameterSet &iConfig ) :
    tokenHydra_( consumes<View<Hydra> >( iConfig.getParameter<InputTag> ( "HydraTag" ) ) ),
    useGenParticles_( iConfig.getParameter<bool>( "UseGenParticles" ) ),
    splitRecHits_( iConfig.getParameter<bool>( "SplitRecHits" ) ),
    debugPrint_( iConfig.getUntrackedParameter<bool>( "DebugPrint" , true ) )
{
    produces<reco::PFClusterCollection>();
}

void HydraFakeClusterBuilder::produce( Event &iEvent, const EventSetup & )
{
    std::cout << " HydraFakeClusterBuilder::produce useGenParticles=" << useGenParticles_ << std::endl;
    Handle<View<Hydra> > HydraHandle;
    iEvent.getByToken(tokenHydra_, HydraHandle);
    assert ( HydraHandle->size() == 1 );
    hydraObj.reset( new HydraWrapper( HydraHandle->ptrAt(0)) );

    auto_ptr<reco::PFClusterCollection> pfClusterCol ( new reco::PFClusterCollection );

    std::vector<Index_t> tracksToBecomeClusters;

    for ( unsigned i = 0 ; i < hydraObj->simTrackSize() ; i++ ) {
        edm::Ptr<SimTrack> currentTrack = hydraObj->simTrack( i );
        if (useGenParticles_) {
            if (!currentTrack->noGenpart() ) {
                tracksToBecomeClusters.push_back(i);
            }
        } else {
            unsigned nHitsImmediate = hydraObj->simHitsFromSimTrack( i, false ).size();
            unsigned nHitsAllDescendants = hydraObj->simHitsFromSimTrack( i, true ).size();
            if ( nHitsImmediate > 0 && nHitsAllDescendants == nHitsImmediate ) {
                tracksToBecomeClusters.push_back(i);
            }
        }
    }

    if (debugPrint_) {
        std::cout << "There are " << tracksToBecomeClusters.size() << " tracks to become clusters: " << std::endl;
        for ( auto i : tracksToBecomeClusters ) {
            edm::Ptr<SimTrack> currentTrack = hydraObj->simTrack( i );
            auto p = currentTrack->momentum();
            std::cout << "   Track for cheated cluster " << i << " has energy eta geantId pdgId " << p.E() << " " << p.Eta() 
                      << " " << currentTrack->trackId() << " " << currentTrack->type() << std::endl;
        }
    }

    for ( auto i : tracksToBecomeClusters ) {
        reco::PFCluster temp;
        edm::Ptr<SimTrack> currentTrack = hydraObj->simTrack( i );
        auto simHits = hydraObj->simHitsFromSimTrack( i, true );
        for ( unsigned j = 0 ; j < hydraObj->recHitSize(); j++ ) {
            edm::Ptr<reco::PFRecHit> currentRecHit = hydraObj->recHit( j );
            if ( hydraObj->hasSimHitFromRecHit( j ) ) {
                auto result = hydraObj->simHitAndFractionFromRecHit( j );
                auto match = std::find(simHits.begin(), simHits.end(), result.first );
                if ( match != simHits.end() ) {
                    if ( debugPrint_ ) {
                        std::cout << "   Track " << i << " matches to rechit " << j << " with fraction " << result.second;
                        std::cout << "     e(track) = " << currentTrack->momentum().E() << " e(hit)=" << currentRecHit->energy();
                        std::cout << std::endl;
                    }
                    // Need to actually add the information to the cluster!
                }
            }
        }
        pfClusterCol->push_back( temp );
    }
        

    iEvent.put ( pfClusterCol );
}

DEFINE_FWK_MODULE( HydraFakeClusterBuilder );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
