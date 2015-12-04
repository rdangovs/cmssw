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
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
//#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

//#include "RecoParticleFlow/PandoraTranslator/plugins/PandoraCMSPFCandProducer.h" // for CalibHGC

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
    double minDebugEnergy_;

    //    CalibHGC m_calibEE, m_calibHEF, m_calibHEB;
    //    bool calibInitialized;
    // https://github.com/lgray/cmssw/blob/topic_cheated_reco/RecoParticleFlow/PandoraTranslator/plugins/PandoraCMSPFCandProducer.cc#L390-L477
};

HydraFakeClusterBuilder::HydraFakeClusterBuilder( const ParameterSet &iConfig ) :
    tokenHydra_( consumes<View<Hydra> >( iConfig.getParameter<InputTag> ( "HydraTag" ) ) ),
    useGenParticles_( iConfig.getParameter<bool>( "UseGenParticles" ) ),
    splitRecHits_( iConfig.getParameter<bool>( "SplitRecHits" ) ),
    debugPrint_( iConfig.getUntrackedParameter<bool>( "DebugPrint" , true ) ),
    minDebugEnergy_( iConfig.getUntrackedParameter<double>( "MinDebugEnergy", 10. ) )
    //    m_calibEE(ForwardSubdetector::HGCEE,"EE",debugPrint), m_calibHEF(ForwardSubdetector::HGCHEF,"HEF",debugPrint), m_calibHEB(ForwardSubdetector::HGCHEB,"HEB",debugPrint), calibInitialized(false),
{
    produces<reco::PFClusterCollection>();
}

void HydraFakeClusterBuilder::produce( Event &iEvent, const EventSetup & )
{
    std::cout << " HydraFakeClusterBuilder::produce useGenParticles=" << useGenParticles_ << std::endl;
    //    if ( splitRecHits_ ) {
    //        throw cms::Exception("NotImplemented")
    //            << "splitRecHits option not available just yet";
    //    }
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
            if (p.E() > minDebugEnergy_) {
                std::cout << "   Track for cheated cluster " << i << " has energy eta geantId pdgId " << p.E() << " " << p.Eta() 
                          << " " << currentTrack->trackId() << " " << currentTrack->type() << std::endl;
            }
        }
    }

    for ( auto i : tracksToBecomeClusters ) {
        reco::PFCluster temp;
        edm::Ptr<SimTrack> currentTrack = hydraObj->simTrack( i );
        auto simHits = hydraObj->simHitsFromSimTrack( i, true );
        for ( unsigned j = 0 ; j < hydraObj->recHitSize(); j++ ) {
            edm::Ptr<reco::PFRecHit> currentRecHit = hydraObj->recHit( j );
            if ( hydraObj->hasSimHitFromRecHit( j ) ) {
                if ( !splitRecHits_ ) {
                    auto result = hydraObj->simHitAndFractionFromRecHit( j );
                    auto match = std::find(simHits.begin(), simHits.end(), result.first );
                    if ( match != simHits.end() ) {
                        if ( false && debugPrint_ ) {
                            std::cout << "   Track " << i << " matches to rechit " << j << " with fraction " << result.second;
                            std::cout << "     e(track) = " << currentTrack->momentum().E() << " e(hit)=" << currentRecHit->energy();
                            std::cout << std::endl;
                        }
                        edm::Ref<reco::PFRecHitCollection> currentRecHitRef = hydraObj->recHitRef( j );
                        if ( false && debugPrint_ ) std::cout << "    energy of currentRecHitRef: " << currentRecHitRef->energy() << std::endl;
                        std::cout << "       Unsplit: we use a fraction of 1, but it's really " << result.second << std::endl;
                        reco::PFRecHitFraction rhf( currentRecHitRef, 1. ); // unsplit case
                        temp.addRecHitFraction( rhf );
                    }
                } else {
                    auto result_vec = hydraObj->simHitsAndFractionsFromRecHit( j );
                    for ( auto result : result_vec ) {
                        auto match = std::find(simHits.begin(), simHits.end(), result.first );
                        if ( match != simHits.end() ) {
                            if ( false && debugPrint_ ) {
                                std::cout << "   Track " << i << " matches to rechit " << j << " with fraction " << result.second;
                                std::cout << "     e(track) = " << currentTrack->momentum().E() << " e(hit)=" << currentRecHit->energy();
                                std::cout << std::endl;
                            }
                            edm::Ref<reco::PFRecHitCollection> currentRecHitRef = hydraObj->recHitRef( j );
                            if ( false && debugPrint_ ) std::cout << "    energy of currentRecHitRef: " << currentRecHitRef->energy() << std::endl;
                            std::cout << "     Split: fraction=" << result.second << std::endl;
                            reco::PFRecHitFraction rhf( currentRecHitRef, result.second ); // split case
                            temp.addRecHitFraction( rhf );
                        }
                    }
                }
            }
        }

        // Cluster position calculation is not log weighted yet!
        float cl_energy = 0.;
        float cl_x =0., cl_y =0., cl_z =0.;
        int n_h = 0;
        for( const reco::PFRecHitFraction& rhf : temp.recHitFractions() ) {
            const reco::PFRecHitRef& refhit = rhf.recHitRef();
            const double rh_fraction = rhf.fraction();
            const double rh_rawenergy = refhit->energy();
            const double rh_energy = rh_rawenergy * rh_fraction;   
            cl_energy += rh_energy;
            auto cl_pos = refhit->position();
            cl_x += cl_pos.x();
            cl_y += cl_pos.y();
            cl_z += cl_pos.z();
            n_h++;
        }
        if (n_h > 0) {
            temp.setEnergy(cl_energy);
            temp.setPosition(math::XYZPoint(cl_x/n_h,cl_y/n_h,cl_z/n_h));
            if ( currentTrack->momentum().E() > minDebugEnergy_ ) {
                if (debugPrint_) std::cout << "Track " << i << " energy=" << currentTrack->momentum().E() << " fake cluster energy=" << temp.energy() << std::endl;
                if (debugPrint_) std::cout << "Track " << i << " eta=" << currentTrack->momentum().Eta() << " fake cluster eta=" << temp.eta() << std::endl;
                if (debugPrint_) std::cout << "Track " << i << " phi=" << currentTrack->momentum().Phi() << " fake cluster phi=" << temp.phi() << std::endl;
            }
            pfClusterCol->push_back( temp );
        } else {
            if ( false ) std::cout << "   cluster has no hits, skipping" << std::endl;
        }
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
