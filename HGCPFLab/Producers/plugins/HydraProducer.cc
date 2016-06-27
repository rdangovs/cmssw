#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
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
#include "DataFormats/HGCRecHit/interface/HGCUncalibratedRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

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
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "SimDataFormats/CaloTest/interface/HGCalTestNumbering.h"

#include <unordered_map>
#include <unordered_set>
#include "TH2.h"
#include "TH3.h"

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
    //EDGetTokenT<HGCUncalibratedRecHitCollection> uncalibRHit_;
    std::vector<edm::InputTag> uncalibRHit_;
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
//     edm::ESHandle<HGCalGeometry> hgchebGeoHandle_; 

    bool debug_;
   
    // TODO???  RecTrack to (simulated) TrackingParticle ???
    //    inputTagtPRecoTrackAsssociation_ = iConfig.getParameter<InputTag>("tPRecoTrackAsssociation");

    std::array<TH2D*,3> h_recHit_E_vs_simHit_E;
    std::array<TH2D*,3> h_recHit_E_vs_simHit_E_FH;
    
    std::array<TH2D*,3> h_recHit_uncalibE_vs_simHit_E;
    std::array<TH2D*,3> h_recHit_uncalibE_vs_simHit_E_FH;
    
    std::array<TH1D*,3> h_recHit_uncalibE;
    std::array<TH1D*,3> h_recHit_uncalibE_FH;

    std::array<TH1D*,3> h_simHit_E;
    std::array<TH1D*,3> h_simHit_E_FH;

    TH2D * h_recLayer_vs_simLayer;
    TH2D * h_recDet_vs_simDet;
    TH2D * h_recWafer_vs_simWafer;   
    TH2D * h_recCell_vs_simCell;
    TH1D * h_frac;
    TH2D * h_cType;
    TH2D * h_recZ_vs_simZ;
};

HydraProducer::HydraProducer( const ParameterSet &iConfig ) :
    tokenGenParticle_( consumes<View<GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleCollection" ) ) ),
    tokenGenBarcode_( consumes<View<Barcode_t> >( iConfig.getParameter<InputTag> ( "GenParticleCollection" ) ) ),
    tokenPFRecTrack_( consumes<View<PFRecTrack> >( iConfig.getParameter<InputTag> ("RecTrackCollection") ) ),
    tokenSimTrack_( consumes<View<SimTrack> >( iConfig.getParameter<InputTag> ("SimTrackCollection") ) ),
    tokenSimVertex_( consumes<View<SimVertex> >( iConfig.getParameter<InputTag> ("SimVertexCollection") ) )
{
	edm::Service<TFileService> fs;     
	h_recHit_E_vs_simHit_E[2] = fs->make<TH2D>("h_recHit_E_vs_simHit_E_300","h_recHit_E_vs_simHit_E",2000,0,0.010,2000,0,0.010);
    h_recHit_E_vs_simHit_E[1] = fs->make<TH2D>("h_recHit_E_vs_simHit_E_200","h_recHit_E_vs_simHit_E",2000,0,0.010,2000,0,0.010);
    h_recHit_E_vs_simHit_E[0] = fs->make<TH2D>("h_recHit_E_vs_simHit_E_100","h_recHit_E_vs_simHit_E",2000,0,0.010,2000,0,0.010);
	h_recHit_E_vs_simHit_E_FH[2] = fs->make<TH2D>("h_recHit_E_vs_simHit_E_FH_300","h_recHit_E_vs_simHit_E_FH",2000,0,0.010,2000,0,0.010);
    h_recHit_E_vs_simHit_E_FH[1] = fs->make<TH2D>("h_recHit_E_vs_simHit_E_FH_200","h_recHit_E_vs_simHit_E_FH",2000,0,0.010,2000,0,0.010);
    h_recHit_E_vs_simHit_E_FH[0] = fs->make<TH2D>("h_recHit_E_vs_simHit_E_FH_100","h_recHit_E_vs_simHit_E_FH",2000,0,0.010,2000,0,0.010);
	h_recHit_uncalibE_vs_simHit_E[2] = fs->make<TH2D>("h_recHit_uncalibE_vs_simHit_E_300","h_recHit_uncalibE_vs_simHit_E",2000,0,0.010,300,0,300);
    h_recHit_uncalibE_vs_simHit_E[1] = fs->make<TH2D>("h_recHit_uncalibE_vs_simHit_E_200","h_recHit_uncalibE_vs_simHit_E",2000,0,0.010,300,0,300);
    h_recHit_uncalibE_vs_simHit_E[0] = fs->make<TH2D>("h_recHit_uncalibE_vs_simHit_E_100","h_recHit_uncalibE_vs_simHit_E",2000,0,0.010,300,0,300);
	h_recHit_uncalibE_vs_simHit_E_FH[2] = fs->make<TH2D>("h_recHit_uncalibE_vs_simHit_E_FH_300","h_recHit_uncalibE_vs_simHit_E_FH",2000,0,0.010,1000,0,300);
    h_recHit_uncalibE_vs_simHit_E_FH[1] = fs->make<TH2D>("h_recHit_uncalibE_vs_simHit_E_FH_200","h_recHit_uncalibE_vs_simHit_E_FH",2000,0,0.010,1000,0,300);
    h_recHit_uncalibE_vs_simHit_E_FH[0] = fs->make<TH2D>("h_recHit_uncalibE_vs_simHit_E_FH_100","h_recHit_uncalibE_vs_simHit_E_FH",2000,0,0.010,1000,0,300);

    h_recHit_uncalibE[2] = fs->make<TH1D>("h_recHit_uncalibE_300","h_recHit_uncalibE",100,0,10);
    h_recHit_uncalibE[1] = fs->make<TH1D>("h_recHit_uncalibE_200","h_recHit_uncalibE",100,0,10);
    h_recHit_uncalibE[0] = fs->make<TH1D>("h_recHit_uncalibE_100","h_recHit_uncalibE",100,0,10);
	h_recHit_uncalibE_FH[2] = fs->make<TH1D>("h_recHit_uncalibE_FH_300","h_recHit_uncalibE_FH",100,0,10);
    h_recHit_uncalibE_FH[1] = fs->make<TH1D>("h_recHit_uncalibE_FH_200","h_recHit_uncalibE_FH",100,0,10);
    h_recHit_uncalibE_FH[0] = fs->make<TH1D>("h_recHit_uncalibE_FH_100","h_recHit_uncalibE_FH",100,0,10);

    h_simHit_E[2] = fs->make<TH1D>("h_simHit_E_300","h_simHit_E",100,0,0.0005);
    h_simHit_E[1] = fs->make<TH1D>("h_simHit_E_200","h_simHit_E",100,0,0.0005);
    h_simHit_E[0] = fs->make<TH1D>("h_simHit_E_100","h_simHit_E",100,0,0.0005);
	h_simHit_E_FH[2] = fs->make<TH1D>("h_simHit_E_FH_300","h_simHit_E_FH",100,0,0.0005);
    h_simHit_E_FH[1] = fs->make<TH1D>("h_simHit_E_FH_200","h_simHit_E_FH",100,0,0.0005);
    h_simHit_E_FH[0] = fs->make<TH1D>("h_simHit_E_FH_100","h_simHit_E_FH",100,0,0.0005);

    h_recLayer_vs_simLayer = fs->make<TH2D>("h_recLayer_vs_simLayer","h_recLayer_vs_simLayer",100,0,100,50,0,50);
    h_recDet_vs_simDet = fs->make<TH2D>("h_recDet_vs_simDet","h_recDet_vs_simDet",10,0,10,10,0,10);
	h_recWafer_vs_simWafer = fs->make<TH2D>("h_recWafer_vs_simWafer","h_recWafer_vs_simWafer",1000,0,1000,1000,0,1000);
	h_recCell_vs_simCell = fs->make<TH2D>("h_recCell_vs_simCell","h_recCell_vs_simCell",260,0,260,260,0,260);
    h_recZ_vs_simZ = fs->make<TH2D>("h_recZ_vs_simZ","h_recZ_vs_simZ",5,-2,3,5,-2,3);
	h_frac = fs->make<TH1D>("h_frac","h_frac",50,0,5);
	h_cType = fs->make<TH2D>("h_cType","h_cType",10,-5,5, 6, -3, 3 );

    inputSimHits_ = iConfig.getParameter<std::vector<InputTag> >("SimHitCollection");
    for( const auto& tag : inputSimHits_ ) {
        consumes<View<PCaloHit> >(tag);
    }
    inputRecHits_ = iConfig.getParameter<std::vector<InputTag> >("HGCRecHitCollection");
    for( const auto& tag : inputRecHits_ ) {
        consumes<View<PFRecHit> >(tag);
    }
    uncalibRHit_ = iConfig.getParameter<std::vector<InputTag> >("HGCalUncalibRecHitCollection");
    for( const auto& tag : uncalibRHit_ ) {
        consumes<HGCUncalibratedRecHitCollection>(tag);
    }
    
    debug_ = iConfig.getUntrackedParameter<bool>("Debug",false);

    produces<std::vector<Hydra> >();
}

void HydraProducer::beginLuminosityBlock( LuminosityBlock const& iLumiBlock, const EventSetup& iSetup ) {
    iSetup.get<CaloGeometryRecord>().get(geoHandle_);
    iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",hgceeGeoHandle_) ; 
    iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",hgchefGeoHandle_) ; 
    //iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive",hgchebGeoHandle_) ; 
}

void HydraProducer::produce( Event &iEvent, const EventSetup & )
{
    auto_ptr<std::vector<Hydra> > output( new std::vector<Hydra> );
    output->emplace_back(); // constructs an empty object


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


    vector<Handle<View<PCaloHit> > > simHits;  
    for( const auto& tag : inputSimHits_ ) {
        if(debug_)std::cout << tag << std::endl;
        simHits.emplace_back( Handle<View<PCaloHit> >() );
        iEvent.getByLabel(tag,simHits.back());
    }


    vector<Handle<HGCUncalibratedRecHitCollection>> uncalibRecHits;
    for( const auto& tag : uncalibRHit_ ) {
        if(debug_)std::cout << tag << std::endl;
		uncalibRecHits.emplace_back( Handle<HGCUncalibratedRecHitCollection>() );
		iEvent.getByLabel(tag, uncalibRecHits.back());
    }

    vector<Handle<View<PFRecHit> > > recHits;
    for( const auto& tag : inputRecHits_ ) {
        if(debug_)std::cout << tag << std::endl;
        recHits.emplace_back( Handle<View<PFRecHit> >() );
        iEvent.getByLabel(tag,recHits.back());
    }

    // setup the reco det id to sim-hit match
    unordered_multimap<uint32_t,tuple<unsigned,unsigned,float> > temp_recoDetIdToSimHit;
    unordered_set<unsigned> reco_detIds;
    int numSimHit=0;
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
//             case HGCHEB:
//                 geom = hgchebGeoHandle_.product();
//                 break;
            default:
                throw cms::Exception("InvalidDetector")
                    << "Got invalid HGC subdet: " << mysubdet;
            }
            const HGCalTopology& topo = geom->topology();
            const HGCalDDDConstants& dddConst = topo.dddConstants();
      
            int layer, cell, sec, subsec, zp;
            int subdet;
            uint32_t intSimId = simId.rawId();
            uint32_t recoDetId = 0;
            if (dddConst.geomMode() == HGCalGeometryMode::Square) {
				HGCalTestNumbering::unpackSquareIndex(intSimId, zp, layer, sec, subsec, cell);
				if(!HGCalTestNumbering::isValidSquare(zp, layer, sec, subsec, layer)) {
					throw cms::Exception("BadGeometry") << "Sim DetID -> Reco DetID transform was invalid!";
				}
            } else {
             
              HGCalTestNumbering::unpackHexagonIndex(simId, subdet, zp, layer, sec, subsec, cell); 
              mysubdet = (ForwardSubdetector)(subdet);
              //sec is wafer and subsec is celltyp
            }
            //skip this hit if after ganging it is not valid
            //std::cout << simHits[i]->ptrAt(j)->id() << "   " << simHits[i]->ptrAt(j)->energy() << "\n";
            numSimHit++;
            //std::cout << "output Vito Cell, Layer from SimHit= " << cell<< " , " << layer << std::endl; 
            std::pair<int,int> recoLayerCell=dddConst.simToReco(cell,layer,sec,topo.detectorType());
            cell  = recoLayerCell.first;
            layer = recoLayerCell.second;
			//std::cout << "output Vito Cell, Layer from recHit= " << cell<< " , " << layer << std::endl; 

            if (layer<0 || cell<0) {
              //hitRefs[i]=std::make_tuple( i, 0, 0.);
              continue;
            }

            //assign the RECO DetId
            DetId id;
            if (dddConst.geomMode() == HGCalGeometryMode::Square) {
				recoDetId = ( ( geom == hgceeGeoHandle_.product() ) ?
					(uint32_t)HGCEEDetId(ForwardSubdetector(mysubdet),zp,layer,sec,subsec,cell) :
					(uint32_t)HGCHEDetId(ForwardSubdetector(mysubdet),zp,layer,sec,subsec,cell)
				);
            } else {
              recoDetId = HGCalDetId(mysubdet,zp,layer,subsec,sec,cell);

            }

            reco_detIds.insert(recoDetId); //unordered set with the detId values or the recHit
            //int waferNumber = ((HGCalDetId)(simId)).wafer();
            float cellTypeL=-1;
            cellTypeL = dddConst.waferTypeL(sec);
            h_cType->Fill(cellTypeL, subsec);
            //std::cout << " cel longitudinal size =  "<< celTypeL  << std::endl;
            temp_recoDetIdToSimHit.emplace(recoDetId,make_tuple(i,j,cellTypeL));
            unordered_multimap<uint32_t,tuple<unsigned,unsigned,float> >::const_iterator it;
            if (debug_) std::cout << " Inserted simHit from detector " << i << " cobined with recHit "<< j <<" recoDetId = " << recoDetId << " in layer " << layer << "  cellTypeL"<< cellTypeL << std::endl;
			//std::cout << cellTypeL << "\n";
        }
    }
    //2147483647
    //std::cout  << "999999999"<<numSimHit << "\n";
	//loop over the multimap
    unordered_multimap<uint32_t,tuple<unsigned,unsigned,float> >::const_iterator it;
    for (it = temp_recoDetIdToSimHit.begin();it != temp_recoDetIdToSimHit.end();++it){
    	//std::cout << "recoDetIdToSimId recHit Id = " << it->first << " i=" << get<0>(it->second) << " j=" << get<1>(it->second) << " frac=" << get<2>(it->second) << std::endl;
    }
    // calculate and store the weights for particles associated to 
    // pcalohits
    for( const unsigned detid : reco_detIds ) {
        auto range = temp_recoDetIdToSimHit.equal_range(detid);
        double e_tot = 0.0;
        for( auto iter = range.first; iter != range.second; ++iter ) {
            const Ptr<PCaloHit> hit = simHits[get<0>(iter->second)]->ptrAt(get<1>(iter->second));
            if( hit->geantTrackId() > 0 ) e_tot += hit->energy();
            //std::cout << " good fraction " << get<2>(iter->second) << std::endl;
        }
        for( auto iter = range.first; iter != range.second; ++iter ) {
            const Ptr<PCaloHit> hit = simHits[get<0>(iter->second)]->ptrAt(get<1>(iter->second));
            if( hit->geantTrackId() > 0 ) {
            	HGCalDetId simId(hit->id());
            	//std::cout << " geant4 track id = " << hit->geantTrackId()  << " sim id " << simId << std::endl;
                float fraction = hit->energy()/e_tot;
				float cellSize = (get<2>(iter->second));
                //std::cout << "final fraction " <<  cellSize <<" , " << fraction << " , "<< (cellSize + fraction) << std::endl;

                make_tuple(get<0>(iter->second),get<1>(iter->second),(cellSize + fraction)).swap(iter->second);

                output->back().setRecoDetIdMatchToSimHit(get<0>(iter->second),hit,detid,(cellSize + fraction) );

                
                ////// IMPORTANT MODIFIED LUCA                
                //output->back().setRecoDetIdMatchToSimHit(get<0>(iter->second),get<1>(iter->second),hit,detid,fraction);
                
                //std::cout << "setRecoDetIdMatchToSimHit i = "<<  get<0>(iter->second) << " j "<< get<1>(iter->second) << std::endl;
                //std::cout << "setRecoDetIdMatchToSimHit associated SimTrack " << hit->geantTrackId() << std::endl; 
            }
        }
    }

    if (false) {
    	//std::cout << "Number of recHit " << recHits[0].size() << "\n";
        int numbHit=0;
        for( const unsigned detid : reco_detIds ) {
            auto range = temp_recoDetIdToSimHit.equal_range(detid);
            DetId recId(detid);
            HGCalDetId recId2(detid);    
            //std::cout <<"--->reco det Id = " << detid << "   from HGCalDetId ==>> "<< recId2 << "\n";// << recId <<  "\n"; 
           	for( unsigned i = 0; i < recHits.size(); ++i ) {
        		for( unsigned j = 0; j < recHits[i]->size(); ++j ) {
					if(recHits[i]->ptrAt(j)->detId() == recId){

                        numbHit++;
                    }
                }
            }
            for( auto iter = range.first; iter != range.second; ++iter ) {
            	const Ptr<PCaloHit> hit = simHits[get<0>(iter->second)]->ptrAt(get<1>(iter->second));
                DetId hitid(hit->id());
            }
        }
    	//std::cout << "Number of recHit " << numbHit << "\n";
    }

    for(unsigned i=0; i<SimTrackHandle->size(); i++) {
        output->back().insertSimTrack(SimTrackHandle->ptrAt(i));
    }
    for(unsigned i=0; i<SimVertexHandle->size(); i++) {
        output->back().insertSimVertex(SimVertexHandle->ptrAt(i));
    }

    for( unsigned i = 0; i < recHits.size(); ++i ) {
        for( unsigned j = 0; j < recHits[i]->size(); ++j ) {
            if (debug_)std::cout << " i=" << i << " j=" << j << " detId=" << recHits[i]->ptrAt(j)->detId() << " subdet=" << (ForwardSubdetector)(i+3) << std::endl;

        	auto range = temp_recoDetIdToSimHit.equal_range(recHits[i]->ptrAt(j)->detId());
        	double e_tot = 0.0;
        	for( auto iter = range.first; iter != range.second; ++iter ) {
            	const Ptr<PCaloHit> hit = simHits[get<0>(iter->second)]->ptrAt(get<1>(iter->second));
            	if( hit->id() > 0 ) e_tot += hit->energy();
        	}
        	for( auto iter = range.first; iter != range.second; ++iter ) {
            	const Ptr<PCaloHit> hit = simHits[get<0>(iter->second)]->ptrAt(get<1>(iter->second));
            	if( hit->id() > 0 ) {
            		
                    HGCalDetId simId(hit->id());            		
      
            		int Sim_subdet, Sim_z, Sim_lay, Sim_wafer, Sim_celltyp, Sim_cell;
                	int mySimsubdet=0;
                      
                    HGCalTestNumbering::unpackHexagonIndex(simId, Sim_subdet, Sim_z, Sim_lay, Sim_wafer, Sim_celltyp, Sim_cell); 
                    mySimsubdet = (ForwardSubdetector)(Sim_subdet);
                    
            		//std::cout << " geant4 track id = " << hit->geantTrackId()  << " sim id " << simId << std::endl;
					//float fraction = hit->energy()/e_tot;
					float fraction = (get<2>(iter->second));
                	//std::cout << "final fraction 2 "  << fraction << std::endl;
            		make_tuple(get<0>(iter->second),get<1>(iter->second),fraction).swap(iter->second);
                    
                    DetId rechitid(recHits[i]->ptrAt(j)->detId());
                    
                    HGCalDetId recId(recHits[i]->ptrAt(j)->detId());
					int layer = recId.layer();
                    //int wafTp = recId.waferType();
                    
					//std::cout << " new fraction = " << fraction <<  "waferType " << wafTp <<  std::endl;
                    if((ForwardSubdetector)rechitid.subdetId()>0)
                    {
                    	h_recLayer_vs_simLayer->Fill(Sim_lay,layer);
                    	h_recDet_vs_simDet->Fill(mySimsubdet,(ForwardSubdetector)rechitid.subdetId() );
                    	h_recWafer_vs_simWafer->Fill(Sim_wafer, recId.wafer());
                    	h_recCell_vs_simCell->Fill(Sim_cell, recId.cell());
                        h_recZ_vs_simZ->Fill(Sim_z, recId.zside());
                        
                    	h_frac->Fill((fraction));
                    }					
            		//output->back().setRecoDetIdMatchToSimHit(get<0>(iter->second),hit,detid,fraction);                                       
            	}
        	}
            
            bool amplitudeSet = false;
            double uncalibAmplitude = -10.;
            
            //match recHit with uncalib
            DetId key(recHits[i]->ptrAt(j)->detId());
            auto iter = uncalibRecHits[key.subdetId()-3]->find(key);
            if( iter != uncalibRecHits[key.subdetId()-3]->end() ) {
                uncalibAmplitude = double(iter->amplitude());
                amplitudeSet = true;
                //std::cout << "YEP YEP YEP YEP" << std::endl;
            } else {
                //std::cout << "NOPE NOPE NOPE NOPE" << std::endl;
            }
            
            if(amplitudeSet){
                if(debug_){
                    //std::cout << " >>> PF = " << recHits[i]->ptrAt(j)->energy() << " detId = " << recHits[i]->ptrAt(j)->detId() << std::endl;
                    //std::cout << " >>> Uncalib = " << uncalibAmplitude << " detId = " << std::endl;
                    
                    //std::cout << " px = " << recHits[i]->ptrAt(j)->position().x() 
                    //        << " py = " << recHits[i]->ptrAt(j)->position().y()
                    //        << " pz = " << recHits[i]->ptrAt(j)->position().z()
                        //<< " Ax = " << recHits[i]->ptrAt(j)->getAxisXYZ().X() 
                        //<< " Ay = " << recHits[i]->ptrAt(j)->getAxisXYZ().Y() 
                        //<< " Az = " << recHits[i]->ptrAt(j)->getAxisXYZ().Z() 
                    //<< std::endl;
                    
                }
            }            
            const auto& hit = *(recHits[i]->ptrAt(j));
            DetId rechitid(hit.detId());

            ForwardSubdetector mysubdet = (ForwardSubdetector)rechitid.subdetId();
            const HGCalGeometry* geom = nullptr;
            switch(mysubdet) {
            case HGCEE:
                geom = hgceeGeoHandle_.product();
                break;
            case HGCHEF:
                geom = hgchefGeoHandle_.product();
                break;            
            default:
                throw cms::Exception("InvalidDetector")
                    << "Got invalid HGC subdet: " << mysubdet;
            }
            const HGCalTopology& topo = geom->topology();
            const HGCalDDDConstants& dddConst = topo.dddConstants();

            int wafTpL = dddConst.waferTypeL(HGCalDetId(rechitid).wafer());
            const double cos_theta = std::abs(std::cos(hit.position().theta()));
            const double corr_e_tot = e_tot*cos_theta;
            const double corr_uncAmpl = uncalibAmplitude*cos_theta;
            const double corr_ehit = recHits[i]->ptrAt(j)->energy()*cos_theta;
            if((ForwardSubdetector)rechitid.subdetId()==3){                
                if( e_tot > 0 ) {
                    h_recHit_E_vs_simHit_E[wafTpL-1]->Fill(corr_e_tot, corr_ehit);
                    h_recHit_uncalibE_vs_simHit_E[wafTpL-1]->Fill(corr_e_tot, corr_uncAmpl);                
                    h_recHit_uncalibE[wafTpL-1]->Fill(corr_uncAmpl);
                    h_simHit_E[wafTpL-1]->Fill(corr_e_tot);
                }
            }else if((ForwardSubdetector)rechitid.subdetId()==4){               
                if( e_tot > 0 ) {
                    h_simHit_E_FH[wafTpL-1]->Fill(corr_e_tot);
                    h_recHit_E_vs_simHit_E_FH[wafTpL-1]->Fill(corr_e_tot, corr_ehit);
                    h_recHit_uncalibE_vs_simHit_E_FH[wafTpL-1]->Fill(corr_e_tot, corr_uncAmpl);
                    h_recHit_uncalibE_FH[wafTpL-1]->Fill(corr_uncAmpl);
                    h_simHit_E_FH[wafTpL-1]->Fill(corr_e_tot);
                }
            }        
            
            unsigned int layer = 999;
            int det=-1;
            try {
                DetId hitid(recHits[i]->ptrAt(j)->detId());
                det = hitid.subdetId();
                layer = GetHGCLayer( hitid, (ForwardSubdetector)hitid.subdetId() );
            } catch ( const cms::Exception& e ) {
                if (debug_)std::cout << "   caught exception " << e.what() << " but moving on" << std::endl;
            }
            if (debug_)std::cout << " Inserted recHit from detector " << det << " in layer " << layer << std::endl;
            
			output->back().insertRecHit(i,recHits[i]->ptrAt(j));
            //std::cout << "	insertRecHit "<< i << "    id " <<  recHits[i]->ptrAt(j)->detId() <<  std::endl;
       		//auto range = temp_recoDetIdToSimHit.equal_range(recHits[i]->ptrAt(j)->detId());
        	// for( auto iter = range.first; iter != range.second; ++iter ) {
//             	const Ptr<PCaloHit> hit = simHits[get<0>(iter->second)]->ptrAt(get<1>(iter->second));
//         		std::cout << "i = "<<  get<0>(iter->second) << " j "<< get<1>(iter->second) << std::endl;
//                 std::cout << "    associated SimTrack " << hit->geantTrackId() << std::endl; 
//             }

        }
    }
    


    for(unsigned i=0; i<GenParticleHandle->size(); i++) {
        output->back().insertGenParticle(GenBarcodeHandle->at(i),GenParticleHandle->ptrAt(i));
    }
    
    if( PFRecTrackHandle.isValid() ) {
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
    }

    iEvent.put( output );
}


unsigned int HydraProducer::GetHGCLayer(const DetId& detid, const ForwardSubdetector& subdet) const {
    unsigned int layer = 0;
	
    layer = (unsigned int) ((HGCalDetId)(detid)).layer() ;
	
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
