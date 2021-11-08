This repository creates jet-level ntuples for training ParticleNet used at the "Jet tagging using the Run3 Scouting data" hackathon (8th-11th Nov 2021).

## Download and compile code

The code was develop in `CMSSW_11_2_1_Patatrack`

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
git fetch alintulu
git checkout alintulu/CMSSW_11_2_1_Patatrack-ScoutPNet

3. clone this repository and compile

git clone git@github.com:alintulu/Run3ScoutingAnalysisTools.git -b particlenet-reco
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
    --python_filename $PYTHON_CFG \
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

source reHLT.sh /eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/Scouting/Run3/ML_210512/GluGluHToBB_M125_masseffects_14TeV_TuneCP5_powheg_pythia8/ML_210512/210602_090726/0000/scouting_75.root /eos/user/a/adlintul/scouting/particlenet/particle_features/reHLT/edm/scouting_75.root scouting_GluGluHToBB.py
```

```diff
6. change the following line in $PYHTON_CFG to differentiate the new objects from the old ones

- process = cms.Process('HLT',Run3)
+ process = cms.Process('reHLT',Run3)
```

```bash
7. run script

cmsRun scouting_GluGluHToBB.py
```

## Create ntuples

```
cmsRun Run3ScoutingAnalysisTools/Analysis/test/ScoutingNanoAOD_cfg.py inputFiles=/eos/user/a/adlintul/scouting/particlenet/particle_features/reHLT/edm/scouting_75.root outputFile=/eos/user/a/adlintul/scouting/particlenet/particle_features/reHLT/nano/scouting_75.root isQCD=False isMC=True useWeights=False GlobalTagMC=112X_mcRun3_2021_realistic_v16
```

Done!
