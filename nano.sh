#!/bin/bash

cmsDriver.py \
    scouting \
    --filein root://cmseos.fnal.gov//store/user/dmahon/condor/RPVSingleStopRun3MC/MINIAOD/MINIAOD-2000_100-1.root \
    --eventcontent NANOAODSIM \
    --datatier NANOAOD \
    --fileout file:nano.root \
    --conditions auto:phase1_2022_realistic \
    --step NANO \
    --geometry DB:Extended \
    --era Run3 \
    --no_exec \
    --mc \
    --customise Run3ScoutingAnalysisTools/Analysis/Run3Scouting_cff.customise \
    --customise_commands "process.NANOAODSIMoutput.outputCommands.extend(['drop edmTriggerResults_*_*_HLT'])" \
    -n 1
