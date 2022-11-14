import FWCore.ParameterSet.Config as cms
from Run3ScoutingAnalysisTools.Analysis.Run3ScoutingPF_cff import *
from Run3ScoutingAnalysisTools.Analysis.Run3ScoutingOther_cff import *

def run3Scouting_customiseMC(process):

   scoutingToReco(process)
   addParticles(process)
   addAK4Jets(process, isMC=True)
   addAK8Jets(process, isMC=True)
   addScouting(process)

   return process


def run3Scouting_customiseData(process):

   scoutingToReco(process)
   addParticles(process)
   addAK4Jets(process, isMC=False)
   addAK8Jets(process, isMC=False)
   addScouting(process)

   return process
