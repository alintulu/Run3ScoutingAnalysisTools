import FWCore.ParameterSet.Config as cms
from Run3ScoutingAnalysisTools.Analysis.Run3ScoutingToReco import *
from Run3ScoutingAnalysisTools.Analysis.Run3ScoutingPFCand import *
from Run3ScoutingAnalysisTools.Analysis.Run3ScoutingAK4PFJet import *
from Run3ScoutingAnalysisTools.Analysis.Run3ScoutingAK8PFJet import *
from Run3ScoutingAnalysisTools.Analysis.Run3Scouting import *

def run3Scouting_customiseMC(process):

   scoutingToReco(process)
   addParticles(process)
   addAK4Jets(process, isMC=True)
   addAK4SoftKillerJets(process, isMC=True)
   addAK4CHSJets(process, isMC=True)
   addAK8Jets(process, isMC=True)
   addAK8SoftKillerJets(process, isMC=True)
   addAK8CHSJets(process, isMC=True)
   addScouting(process)

   return process


def run3Scouting_customiseData(process):

   scoutingToReco(process)
   addParticles(process)
   addAK4Jets(process, isMC=False)
   addAK4SoftKillerJets(process, isMC=False)
   addAK4CHSJets(process, isMC=False)
   addAK8Jets(process, isMC=False)
   addAK8SoftKillerJets(process, isMC=False)
   addAK8CHSJets(process, isMC=False)
   addScouting(process)

   return process
