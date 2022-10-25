## Setup

```
cmsrel CMSSW_12_4_3
cd CMSSW_12_4_3/src
cmsenv
git clone git@github.com:alintulu/Run3ScoutingAnalysisTools.git -b event-vars
scram b -j 8
```

## Run

### Data

```
cmsRun Run3ScoutingAnalysisTools/Analysis/test/ScoutingNanoAOD_cfg.py \
       isMC=False \
       GlobalTagData=124X_dataRun3_HLT_v4 \
       inputFiles=/store/data/Run2022F/ScoutingPFRun3/RAW/v1/000/360/919/00000/c650b47a-4964-4fc2-9e6d-e7e3073f51f9.root \
       maxEvents=10
```
