import FWCore.ParameterSet.Config as cms

#### PF CLUSTER ECAL ####

#energy corrector for corrected cluster producer
imagingClusterHGCal =  cms.EDProducer(
    "HGCalClusterTestProducer",
    detector = cms.string("EE"),
    doSharing = cms.bool(False),
    deltac = cms.double(2.),
    ecut = cms.double(0.060),
    kappa = cms.double(12.0),
    verbosity = cms.untracked.uint32(3)
    )

