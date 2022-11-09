import FWCore.ParameterSet.Config as cms
from Run3ScoutingAnalysisTools.Analysis.Run3ScoutingPF_cff import *
from Run3ScoutingAnalysisTools.Analysis.Run3ScoutingOther_cff import *

def customise(process):

   scoutingToReco(process)
   addParticles(process)
   addAK4Jets(process)
   addAK8Jets(process)
   addScouting(process)

   process.particleTask = cms.Task(process.pfcands,
                                   process.ak4Jets,
                                   process.particleTable,
                                   process.ak4JetTable,
                                   #process.ak4MatchGenTable,
                                   process.ak8Jets,
                                   process.ak8JetTable,
                                   #process.ak8MatchGenTable,
                                   process.pfParticleNetMassRegressionJetTagInfos,
                                   process.pfParticleNetJetTags,
                                   process.pfParticleNetDiscriminatorsJetTags,
                                   process.pfParticleNetMassRegressionJetTags,
                                   process.run3ScoutingTable,
  )

   process.schedule.associate(process.particleTask)

   #process.load('PhysicsTools.NanoAOD.triggerObjects_cff')
   #process.schedule.associate(process.triggerObjectTablesTask)

   return process


