#!/bin/bash

cmsDriver.py \
    scouting \
    --filein root://cmseos.fnal.gov//store/user/dmahon/condor/RPVSingleStopRun3MC/MINIAOD/MINIAOD-2000_100-1.root \
    --eventcontent NANOAODGEN \
    --datatier NANOAOD \
    --fileout file:nanogen.root \
    --conditions auto:phase1_2022_realistic \
    --step NANOGEN \
    --geometry DB:Extended \
    --era Run3 \
    --no_exec \
    --mc \
    --customise Run3ScoutingAnalysisTools/Analysis/Run3Scouting_cff.customise \
    --customise_commands "process.NANOAODGENoutput.outputCommands.extend(['keep edmTriggerResults_*_*_DQM'])" \
    -n 1
