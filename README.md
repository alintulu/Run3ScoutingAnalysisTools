This repository creates jet-level ntuples for training ParticleNet

# How to run

## Download and compile code

The code runs on `CMSSW_11_2_1_Patatrack`

```
1. prepare CMSSW release

cmsrel CMSSW_11_2_1_Patatrack
cd CMSSW_11_2_1_Patatrack/src
cmsenv
git cms-init
git cms-addpkg HLTrigger/JetMET DataFormats/Scouting
```

Changes were made to the HLT Scouting producer, hence we need to re-run the Scouting reconstruction with the updated code

```
2. download the updated CMSSW code

git remote add alintulu git@github.com:alintulu/cmssw.git
git checkout alintulu CMSSW_11_2_1_Patatrack-ScoutPNet

3. clone this repository and compile

git clone git@github.com:alintulu/Run3ScoutingAnalysisTools.git -b particlenet-pat
scram b
```

## Re-run HLT step

With the code compiled we are ready to re-run the Scouting reconstruction

```bash
4. Create the following file at $CMSSW_BASE/src

$ cat reHLT.sh
#!/bin/bash

INPUT_FILE=$1
OUTPUT_FILE=$2
PYTHON_CFG=$3

cmsDriver.py \
    reHLT \
    --python_filename PYTHON_CFG \
    --eventcontent RAWMINIAODSIM \
    --customise HLTrigger/Configuration/customizeHLTforPatatrack.customizeHLTforPatatrackTriplets \
    --filein file:$INPUT_FILE \
    --fileout file:$OUTPUT_FILE \
    --conditions 112X_mcRun3_2021_realistic_v16 \
    --step HLT:GRun \
    --geometry DB:Extended \
    --era Run3 \
    --no_exec \
    --mc \
    -n 10

5. run the file

source reHLT.py file:/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/Scouting/Run3/ML_210512/GluGluHToBB_M125_masseffects_14TeV_TuneCP5_powheg_pythia8/ML_210512/210602_090726/0000/scouting_75.root file:/eos/user/a/adlintul/scouting/particlenet/particle_features/reHLT/edm/scouting_75.root scoutingPF_default.py
```

```diff
6. change the following line in $PYHTON_CFG to differentiate the new objects from the old ones

- process = cms.Process('HLT',Run3)
+ process = cms.Process('reHLT',Run3)
```

## Create ntuples

```
cmsRun Run3ScoutingAnalysisTools/Analysis/test/ScoutingNanoAOD_cfg.py inputFile=/eos/user/a/adlintul/scouting/particlenet/particle_features/reHLT/edm/scouting_75.root outputFile=/eos/user/a/adlintul/scouting/particlenet/particle_features/reHLT/nano/scouting_75.root
```

Done!
