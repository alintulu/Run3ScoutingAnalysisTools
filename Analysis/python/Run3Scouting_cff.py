import FWCore.ParameterSet.Config as cms
from Run3ScoutingAnalysisTools.Analysis.Run3ScoutingPF_cff import *
from Run3ScoutingAnalysisTools.Analysis.Run3ScoutingOther_cff import *

def customise(process):

   scoutingToReco(process)
   addParticles(process)
   addAK4Jets(process)
   addAK8Jets(process)
   addScouting(process)

   return process


