import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

def addL1(process, path=None):

   process.l1bits = cms.EDProducer("L1TriggerResultsConverter",
         src = cms.InputTag("gtStage2Digis"),
         legacyL1 = cms.bool(False),
         storeUnprefireableBit = cms.bool(True),
         src_ext = cms.InputTag("simGtExtUnprefireable")
   )
  
   process.l1Task = cms.Task(process.l1bits)

   if path is None:
      process.schedule.associate(process.l1Task)
   else:
      getattr(process, path).associate(process.l1Task) 
