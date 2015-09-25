// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HGCPFLab/DataFormats/interface/HighRapidityDevRecoAssociation.h"

using namespace std;
using namespace edm;

typedef HighRapidityDevRecoAssociation HyDRA;

class HighRapidityDevRecoAssociationProducer : public edm::EDProducer
{
public:
    HighRapidityDevRecoAssociationProducer( const edm::ParameterSet & );
    
private:
    void produce( edm::Event &, const edm::EventSetup & ) override;
    
  };

HighRapidityDevRecoAssociationProducer::HighRapidityDevRecoAssociationProducer( const ParameterSet &iConfig )
{
    produces<vector<HyDRA> >();
}

void HighRapidityDevRecoAssociationProducer::produce( Event &evt, const EventSetup & )
{
    std::auto_ptr<vector<HyDRA> > output( new vector<HyDRA> );
    output->emplace_back();
    evt.put( output );
}

DEFINE_FWK_MODULE( HighRapidityDevRecoAssociationProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
