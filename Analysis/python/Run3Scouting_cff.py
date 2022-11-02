import FWCore.ParameterSet.Config as cms
from Run3ScoutingAnalysisTools.Analysis.Run3ScoutingPF_cff import *

def customise(process):

   scoutingToReco(process)
   addParticles(process)
   addAK4Jets(process)
   addAK8Jets(process)

   process.particleTask = cms.Task(process.pfcands,
                                   process.ak4Jets,
                                   process.particleTable,
                                   process.ak4MatchGen,
                                   process.ak4JetTable,
                                   process.ak8Jets,
                                   process.ak8MatchGen,
                                   process.ak8JetTable,
  )

   process.schedule.associate(process.particleTask)

   #process.load('PhysicsTools.NanoAOD.triggerObjects_cff')
   #process.schedule.associate(process.triggerObjectTablesTask)

   return process


