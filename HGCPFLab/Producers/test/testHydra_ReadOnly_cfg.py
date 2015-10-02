import FWCore.ParameterSet.Config as cms

process = cms.Process('TESTHYDRA')

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:step3_HydraOnly.root'))

process.ExampleReader = cms.EDProducer("ExampleHydraPFProducer",HydraTag=cms.InputTag("Hydra"))

process.p = cms.Path(process.ExampleReader)

