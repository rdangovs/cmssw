import FWCore.ParameterSet.Config as cms

process = cms.Process('TESTHYDRA')

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:step3_HydraOnly.root'))

process.ExampleReader = cms.EDProducer("ExampleHydraPFProducer",HydraTag=cms.InputTag("Hydra"))

process.FakeClusterGen = cms.EDProducer("HydraFakeClusterBuilder",HydraTag=cms.InputTag("Hydra"),
                                        SplitRecHits=cms.bool(False),
                                        UseGenParticles=cms.bool(True),
                                        MinDebugEnergy=cms.untracked.double(30.)
                                       )

process.FakeClusterCaloFace = cms.EDProducer("HydraFakeClusterBuilder",HydraTag=cms.InputTag("Hydra"),
                                             SplitRecHits=cms.bool(False),
                                             UseGenParticles=cms.bool(False),
                                             MinDebugEnergy=cms.untracked.double(30.)
                                        )

#process.testSequence = cms.Sequence(process.ExampleReader+process.FakeClusterGen+process.FakeClusterCaloFace)
process.testSequence = cms.Sequence(process.FakeClusterGen+process.FakeClusterCaloFace)
process.p = cms.Path(process.testSequence)

