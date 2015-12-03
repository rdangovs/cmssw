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

// We probably don't need all of these
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/FlatTrd.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include <unordered_map>
#include <unordered_set>

using namespace std;
using namespace edm;
using namespace reco;

class HydraProducer : public EDProducer
{
public:
    HydraProducer( const ParameterSet & );
    
private:
    void produce( Event &, const EventSetup & ) override;
    void beginLuminosityBlock( edm::LuminosityBlock const& , const edm::EventSetup& ) override;

    unsigned int GetHGCLayer(const DetId& detid, const ForwardSubdetector& subdet) const;

    //    EDGetTokenT<View<PFRecHit> > tokenHGCRecHit_;
    std::vector<edm::InputTag> inputRecHits_;
    EDGetTokenT<View<GenParticle> > tokenGenParticle_;
    EDGetTokenT<View<Barcode_t> > tokenGenBarcode_;
    EDGetTokenT<View<PFRecTrack> > tokenPFRecTrack_;
    EDGetTokenT<View<SimTrack> > tokenSimTrack_;
    EDGetTokenT<View<SimVertex> > tokenSimVertex_;
    //    EDGetTokenT<RecoToSimCollection> tokenRecoToSim_;
    std::vector<edm::InputTag> inputSimHits_;

    edm::ESHandle<CaloGeometry> geoHandle_;
    edm::ESHandle<HGCalGeometry> hgceeGeoHandle_; 
    edm::ESHandle<HGCalGeometry> hgchefGeoHandle_; 
    edm::ESHandle<HGCalGeometry> hgchebGeoHandle_; 

    bool debug_;
   
    // TODO???  RecTrack to (simulated) TrackingParticle ???
    //    inputTagtPRecoTrackAsssociation_ = iConfig.getParameter<InputTag>("tPRecoTrackAsssociation");

};

HydraProducer::HydraProducer( const ParameterSet &iConfig ) :
    //    tokenHGCRecHit_( consumes<View<PFRecHit> >( iConfig.getParameter<InputTag> ( "HGCRecHitCollection" ) ) ),
    tokenGenParticle_( consumes<View<GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleCollection" ) ) ),
    tokenGenBarcode_( consumes<View<Barcode_t> >( iConfig.getParameter<InputTag> ( "GenParticleCollection" ) ) ),
    tokenPFRecTrack_( consumes<View<PFRecTrack> >( iConfig.getParameter<InputTag> ("RecTrackCollection") ) ),
    tokenSimTrack_( consumes<View<SimTrack> >( iConfig.getParameter<InputTag> ("SimTrackCollection") ) ),
    tokenSimVertex_( consumes<View<SimVertex> >( iConfig.getParameter<InputTag> ("SimVertexCollection") ) )
    //    tokenRecoToSim_( consumes<RecoToSimCollection>( iConfig.getParameter<InputTag> ("RecoToSimCollection") ) )
{
    inputSimHits_ = iConfig.getParameter<std::vector<InputTag> >("SimHitCollection");
    inputRecHits_ = iConfig.getParameter<std::vector<InputTag> >("HGCRecHitCollection");
    debug_ = iConfig.getUntrackedParameter<bool>("Debug",true);

    produces<std::vector<Hydra> >();
}

void HydraProducer::beginLuminosityBlock( LuminosityBlock const& iLumiBlock, const EventSetup& iSetup ) {
    iSetup.get<CaloGeometryRecord>().get(geoHandle_);
    iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",hgceeGeoHandle_) ; 
    iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",hgchefGeoHandle_) ; 
    iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive",hgchebGeoHandle_) ; 
}

void HydraProducer::produce( Event &iEvent, const EventSetup & )
{
    auto_ptr<std::vector<Hydra> > output( new std::vector<Hydra> );
    output->emplace_back(); // constructs an empty object

    //    Handle<View<PFRecHit> > HGCRecHitHandle;
    //    iEvent.getByToken(tokenHGCRecHit_, HGCRecHitHandle);
    Handle<View<GenParticle> > GenParticleHandle;
    iEvent.getByToken(tokenGenParticle_, GenParticleHandle);
    Handle<View<Barcode_t> > GenBarcodeHandle;
    iEvent.getByToken(tokenGenBarcode_, GenBarcodeHandle);
    Handle<View<PFRecTrack> > PFRecTrackHandle;
    iEvent.getByToken(tokenPFRecTrack_, PFRecTrackHandle);
    Handle<View<SimTrack> > SimTrackHandle;
    iEvent.getByToken(tokenSimTrack_, SimTrackHandle);
    Handle<View<SimVertex> > SimVertexHandle;
    iEvent.getByToken(tokenSimVertex_, SimVertexHandle);
    //    Handle<RecoToSimCollection > RecoToSimHandle;
    //    iEvent.getByToken(tokenRecoToSim_,RecoToSimHandle);
    vector<Handle<View<PCaloHit> > > simHits;  
    for( const auto& tag : inputSimHits_ ) {
        std::cout << tag << std::endl;
        simHits.emplace_back( Handle<View<PCaloHit> >() );
        iEvent.getByLabel(tag,simHits.back());
    }

    vector<Handle<View<PFRecHit> > > recHits;
    for( const auto& tag : inputRecHits_ ) {
        std::cout << tag << std::endl;
        recHits.emplace_back( Handle<View<PFRecHit> >() );
        iEvent.getByLabel(tag,recHits.back());
    }

    // setup the reco det id to sim-hit match
    unordered_multimap<uint32_t,tuple<unsigned,unsigned,float> > temp_recoDetIdToSimHit;
    unordered_set<unsigned> reco_detIds;
    for( unsigned i = 0; i < simHits.size(); ++i ) {
        for( unsigned j = 0; j < simHits[i]->size(); ++j ) {
            output->back().insertSimHit(i,simHits[i]->ptrAt(j));
            HGCalDetId simId(simHits[i]->ptrAt(j)->id());
            ForwardSubdetector mysubdet = (ForwardSubdetector)(i+3);
            const HGCalGeometry* geom = nullptr;
            switch(mysubdet) {
            case HGCEE:
                geom = hgceeGeoHandle_.product();
                break;
            case HGCHEF:
                geom = hgchefGeoHandle_.product();
                break;
            case HGCHEB:
                geom = hgchebGeoHandle_.product();
                break;
            default:
                throw cms::Exception("InvalidDetector")
                    << "Got invalid HGC subdet: " << mysubdet;
            }
            const HGCalTopology& topo = geom->topology();
            const HGCalDDDConstants& dddConst = topo.dddConstants();
      
            int layer(simId.layer()), cell(simId.cell());
            pair<int,int> recoLayerCell = dddConst.simToReco(cell,layer,topo.detectorType());
            cell  = recoLayerCell.first;
            layer = recoLayerCell.second;
            if(layer < 0) continue;
      
            uint32_t recoDetId = ( ( geom == hgceeGeoHandle_.product() ) ?
                                   (uint32_t)HGCEEDetId(ForwardSubdetector(mysubdet),simId.zside(),layer,simId.sector(),simId.subsector(),cell) :
                                   (uint32_t)HGCHEDetId(ForwardSubdetector(mysubdet),simId.zside(),layer,simId.sector(),simId.subsector(),cell)
                                   );
            reco_detIds.insert(recoDetId);
            temp_recoDetIdToSimHit.emplace(recoDetId,make_tuple(i,j,0.0f));
            std::cout << " Inserted simHit from detector " << i << " in layer " << layer << std::endl;
        }
    }

    // calculate and store the weights for particles associated to 
    // pcalohits
    for( const unsigned detid : reco_detIds ) {
        auto range = temp_recoDetIdToSimHit.equal_range(detid);
        double e_tot = 0.0;
        for( auto iter = range.first; iter != range.second; ++iter ) {
            const Ptr<PCaloHit> hit = simHits[get<0>(iter->second)]->ptrAt(get<1>(iter->second));
            if( hit->geantTrackId() > 0 ) e_tot += hit->energy();
        }
        for( auto iter = range.first; iter != range.second; ++iter ) {
            const Ptr<PCaloHit> hit = simHits[get<0>(iter->second)]->ptrAt(get<1>(iter->second));
            if( hit->geantTrackId() > 0 ) {
                float fraction = hit->energy()/e_tot;
                make_tuple(get<0>(iter->second),get<1>(iter->second),fraction).swap(iter->second);
                output->back().setRecoDetIdMatchToSimHit(get<0>(iter->second),hit,detid,fraction);
            }
        }
    }

    if (debug_) {
        for( const unsigned detid : reco_detIds ) {
            auto range = temp_recoDetIdToSimHit.equal_range(detid);
            for( auto iter = range.first; iter != range.second; ++iter ) {
                cout << "  SimHit detIndex=" << get<0>(iter->second) << " hitIndex=" << get<1>(iter->second) << " fraction=" << get<2>(iter->second) << endl;
            }
        }
    }

    for(unsigned i=0; i<SimTrackHandle->size(); i++) {
        output->back().insertSimTrack(SimTrackHandle->ptrAt(i));
    }
    for(unsigned i=0; i<SimVertexHandle->size(); i++) {
        output->back().insertSimVertex(SimVertexHandle->ptrAt(i));
    }
    for( unsigned i = 0; i < recHits.size(); ++i ) {
        for( unsigned j = 0; j < recHits[i]->size(); ++j ) {
            std::cout << " i=" << i << " j=" << j << " detId=" << recHits[i]->ptrAt(j)->detId() << " subdet=" << (ForwardSubdetector)(i+3) << std::endl;
            unsigned int layer = 999;
            try {
                layer = GetHGCLayer( recHits[i]->ptrAt(j)->detId(), (ForwardSubdetector)(i+3));
            } catch ( const cms::Exception& e ) {
                std::cout << "   caught exception " << e.what() << " but moving on" << std::endl;
            }
            std::cout << " Inserted recHit from detector " << i << " in layer " << layer << std::endl;
            output->back().insertRecHit(i,recHits[i]->ptrAt(j));
        }
    }
    for(unsigned i=0; i<GenParticleHandle->size(); i++) {
        output->back().insertGenParticle(GenBarcodeHandle->at(i),GenParticleHandle->ptrAt(i));
    }
    for(unsigned i=0; i<PFRecTrackHandle->size(); i++) {
        output->back().insertTrack(PFRecTrackHandle->ptrAt(i));

        // TODO???  RecTrack to (simulated) TrackingParticle ???
        /*
        RefToBase<Track> tr = PFRecTrackHandle->ptrAt(i)->trackRef();
        const RecoToSimCollection pRecoToSim = *(rectosimCollection.product());
        if(pRecoToSim.find(tr) != pRecoToSim.end()){
            vector<pair<TrackingParticleRef, double> > tp = pRecoToSim[tr];
            TrackingParticleRef tpr = tp.begin()->first;
        */
    }


    iEvent.put( output );
}

unsigned int HydraProducer::GetHGCLayer(const DetId& detid, const ForwardSubdetector& subdet) const {
    unsigned int layer = 0;
    if(subdet==ForwardSubdetector::HGCEE) {
        layer = (unsigned int) ((HGCEEDetId)(detid)).layer() ;
    }
    else if(subdet==ForwardSubdetector::HGCHEF){
        layer = (unsigned int) ((HGCHEDetId)(detid)).layer() ;
    }
    else if(subdet==ForwardSubdetector::HGCHEB){
        layer = (unsigned int) ((HGCHEDetId)(detid)).layer() ;
    }
    return layer;
}


DEFINE_FWK_MODULE( HydraProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
