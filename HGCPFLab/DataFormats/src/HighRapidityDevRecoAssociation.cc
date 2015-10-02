#include "HGCPFLab/DataFormats/interface/HighRapidityDevRecoAssociation.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Utilities/interface/Exception.h"

HighRapidityDevRecoAssociation::HighRapidityDevRecoAssociation()
{}

HighRapidityDevRecoAssociation::~HighRapidityDevRecoAssociation()
{}

void HighRapidityDevRecoAssociation::insertTrack(const edm::Ptr<reco::PFRecTrack> & trk) {
    m_trackPtrs.push_back(trk);
}

void HighRapidityDevRecoAssociation::insertRecHit(const edm::Ptr<reco::PFRecHit> & hit) {
    m_recHitPtrs.push_back(hit);
}

void HighRapidityDevRecoAssociation::insertGenParticle(Barcode_t barcode, const edm::Ptr<reco::GenParticle> & gp) 
{
    m_genParticlePtrs.push_back(gp);
    m_genParticleBarcodes.push_back(barcode);
}

void HighRapidityDevRecoAssociation::insertSimTrack(const edm::Ptr<SimTrack> & st)
{
    Index_t i = m_simTrackPtrs.size();
    /*
    std::cout << "HighRapidityDevRecoAssociation::insertSimTrack" << std::endl;
    std::cout << "  Index: " << i << std::endl;
    std::cout << "  TrackId: " << st->trackId() << std::endl;
    std::cout << "  vertIndex: " << st->vertIndex() << std::endl;
    std::cout << "  genpartIndex: " << st->genpartIndex() << std::endl;
    */
    /*
    if ( (int)st->trackId() == st->genpartIndex() ) {
        throw cms::Exception( "DuplicateBarcode" ) << "trackId()  == genpartIndex() == " << st->genpartIndex() << std::endl;
    }
    for (auto it = m_simTrackBarcodes.begin() ; it != m_simTrackBarcodes.end() ; it++) {
        if (it->first == (int)st->trackId()) {
            throw cms::Exception( "DuplicateBarcode" ) << "trackId() is giving a barcode that is already in m_simTrackBarcodes: " << st->trackId();
        }
        if (it->first == st->genpartIndex()) {
            throw cms::Exception( "DuplicateBarcode" ) << "genpartIndex() is giving a barcode that is already in m_simTrackBarcodes: " << st->genpartIndex();
        }
    }
    */
    m_simTrackPtrs.push_back(st);
    m_simTrackBarcodes.push_back(st->trackId());
    if( !st->noVertex() ) {
        m_simVertexAndSimTracks.push_back(std::make_pair(st->vertIndex(), i));
    }
    if( !st->noGenpart() ) {
        m_simTrackGenBarcodes.push_back(st->genpartIndex());
    }
}

void HighRapidityDevRecoAssociation::insertSimVertex(const edm::Ptr<SimVertex> & sv)
{
    Index_t i = m_simVertexPtrs.size();
    m_simVertexPtrs.push_back(sv);
    m_simVertexBarcodes.push_back(sv->vertexId());
    if (!sv->noParent()) {
        m_simTrackAndSimVertex.push_back( std::make_pair(sv->parentIndex(), i ) );
        m_simVertexAndSimTrackParent.push_back( std::make_pair( sv->vertexId(), sv->parentIndex() ) );
    }
}

void HighRapidityDevRecoAssociation::insertSimHit(Index_t iCol,const edm::Ptr<PCaloHit> & sh)
{
    m_simHitPtrs[iCol].push_back(sh);
    m_simHitBarcodes[iCol].push_back(sh->geantTrackId());
}

void HighRapidityDevRecoAssociation::setRecoDetIdMatchToSimHit(Index_t iCol, const edm::Ptr<PCaloHit> & sh, RecoDetId_t detid, float fraction)
{
    auto iter = find(m_simHitBarcodes[iCol].begin(),m_simHitBarcodes[iCol].end(),sh->geantTrackId());
    if ( iter == m_simHitBarcodes[iCol].end() ) {
        throw cms::Exception( "Missing SimHit" ) << " We don't have the requested barcode for this simhit, should be in the list already " << sh->geantTrackId();
    }
    Index_t iSimHit = std::distance(m_simHitBarcodes[iCol].begin(), iter);
    m_recoDetIdAndSimHits.push_back(std::make_pair(detid,std::make_pair(std::make_pair(iCol,iSimHit),fraction)));
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

