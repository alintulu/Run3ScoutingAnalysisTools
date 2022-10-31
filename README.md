## Download, compile and run code

The code was develop in `CMSSW_12_4_5`.

1. Log into lxplus
2. Prepare the CMSSW release

```
cmsrel CMSSW_12_4_5
cd CMSSW_12_4_5/src
cmsenv
```

3. Clone this repository and compile

```
git clone git@github.com:alintulu/Run3ScoutingAnalysisTools.git -b nanoaod
scram b -j 8
```
4. Run locally over one file

```
cmsRun Run3ScoutingAnalysisTools/Analysis/test/ScoutingNanoAOD_cfg.py inputFiles=/store/data/Run2022F/ScoutingPFRun3/RAW/v1/000/360/919/00000/c650b47a-4964-4fc2-9e6d-e7e3073f51f9.root maxEvents=10
```

Done!
