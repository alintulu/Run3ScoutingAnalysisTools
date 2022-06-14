_This repository creates jet-level ntuples for training mass regression with ParticleNet. It was originally created for the "Jet tagging using the Run3 Scouting data" hackathon on November 8-11 2021. The code and instructions were updated on June 14 2022 for more completeness._

## Download, compile and run code

The code was develop in `CMSSW_12_3_0`.

1. Log into lxplus
2. Prepare the CMSSW release

```
cmsrel CMSSW_12_3_0
cd CMSSW_12_3_0/src
cmsenv
```

3. Clone this repository and compile

```
git clone git@github.com:alintulu/Run3ScoutingAnalysisTools.git -b ak8jet-massregression
scram b -j 8
```
4. Run locally over one file

```
cmsRun Run3ScoutingAnalysisTools/Analysis/test/ScoutingNanoAOD_cfg.py inputFiles=file:/eos/cms/store/group/ml/Tagging4ScoutingHackathon/Adelina/DeepNtuples/12_3_0/ScoutingAK4-v00/BulkGraviton_hh_GF_HH_14TeV_TuneCP5_pythia8/DeepNtuplesAK4-v00/220502_112159/0000/mini_1.root outputFile=test.root maxEvents=1
```

Change the argument `maxEvents` to -1 to run over all events.

Done!
