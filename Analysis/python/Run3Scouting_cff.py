import FWCore.ParameterSet.Config as cms
from Run3ScoutingAnalysisTools.Analysis.Run3ScoutingPF_cff import *

def customise(process):

   scoutingToReco(process)
   addParticles(process)
   addJets(process)

   process.particleTask = cms.Task(process.pfcands,
                                   process.recojets,
                                   process.particleTable,
                                   process.jetTable)

   process.schedule.associate(process.particleTask)

   #process.load('PhysicsTools.NanoAOD.triggerObjects_cff')
   #process.schedule.associate(process.triggerObjectTablesTask)

   return process


