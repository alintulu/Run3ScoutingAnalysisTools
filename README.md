This repository creates jet-level ntuples for training flavour jet tagging with ParticleNet (originally created for the "Jet tagging using the Run3 Scouting data" hackathon 8th-11th Nov 2021).

# Start here

Every step described in this README has be tested on lxplus only. There are three main steps to the full worflow: 

1. [Creating the ntuples](#creating-the-ntuples)
2. [Training the network](#training-the-network)
3. [Testing the network performance](#testing-the-network-performance)

Next I'll try to describe each step as carefully as possible.

## Creating the ntuples

The first step is to create the jet-level ntuples which are used for training and testing the network. If your Scouting samples are produced in the 12_2_0 or a newer CMSSW release you can skip to [Running the ntupiliser](#running-the-ntupiliser). If this is not the case please continue reading here.

These are the samples we used:
- BulkGraviton to HH to 4Q (`/BulkGraviton_hh_GF_HH_14TeV_TuneCP5_pythia8/mkomm-ML_210512-d1606bbcd4ad268d18a3191de15a9732/USER`)
- QCD (`/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/mkomm-ML_210512-d1606bbcd4ad268d18a3191de15a9732/USER`)
- TTbar (`/TTbar_TuneCP5_14TeV-pythia8/mkomm-ML_210512-d1606bbcd4ad268d18a3191de15a9732/USER`)

These samples are created in the 11_2_1_Patatrack CMSSW release and are in the format of RAWMINIAODSIM. Since they are RAW we are able to re-run the HLT step. This is necessary since the 11_2_1_Patatrack Scouting format does not include all necessary information. The correct Scouting format is available in the 12_2_0 release and newer releases. Next I will describe how we reconstruct the samples above in the 12_3_0 release.

### Re-running the HLT step

1. Prepare the CMSSW release on lxplus

```
cmsrel CMSSW_12_3_0
cd CMSSW_12_3_0/src
cmsenv
```

2. Create a file called `reHLT.sh` at `$CMSSW_BASE/src`

```
lxplus> cd $CMSSW_BASE/src
lxplus> cat reHLT.sh
#!/bin/bash

INPUT_FILE=$1
OUTPUT_FILE=$2
PYTHON_CFG=$3

cmsDriver.py \
    reHLT \
    --processName reHLT \
    --python_filename $PYTHON_CFG \
    --eventcontent RAWMINIAODSIM \
    --customise HLTrigger/Configuration/customizeHLTforPatatrack.customizeHLTforPatatrackTriplets \
    --filein file:$INPUT_FILE \
    --fileout file:$OUTPUT_FILE \
    --conditions auto:phase1_2021_realistic \
    --step HLT:GRun \
    --geometry DB:Extended \
    --era Run3 \
    --no_exec \
    --mc \
    -n 10
 ```

3. Create the HLT config file by running `reHLT.sh`

```
INPUT=/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/Scouting/Run3/ML_210512/GluGluHToBB_M125_masseffects_14TeV_TuneCP5_powheg_pythia8/ML_210512/210602_090726/0000/scouting_75.root
OUTPUT=<INSERT YOUR OWN PATH>
source reHLT.sh $INPUT $OUTPUT reHLT.py
```

4. Re-run the HLT step by running `reHLT.py`. Here you are two options

    a. Run the file locally over one file

    ```
    cmsRun reHLT.py
    ```

    b. Run over all samples using CRAB

    Here is an example of a CRAB config file you can use. You need to edit `outLFNDirBase` and `storageSite` to match the site you have access to write to. You can find more information about CRAB [here](https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile). The other fields can stay as they are.

    ```bash
    from WMCore.Configuration import Configuration
    config = Configuration()
    config.section_('General')
    config.General.transferLogs = False
    config.General.transferOutputs = True
    config.General.workArea = 'crab_projects'
    config.General.requestName = 'BulkGraviton_hh_GF_HH_14TeV_TuneCP5_pythia8'
    config.section_('JobType')
    config.JobType.numCores = 1
    config.JobType.sendExternalFolder = True
    config.JobType.pyCfgParams = ['outputFile=mini.root']
    config.JobType.pluginName = 'Analysis'
    config.JobType.allowUndistributedCMSSW = True
    config.JobType.psetName = 'reHLT.py'
    config.JobType.maxMemoryMB = 2000
    config.section_('Data')
    config.Data.inputDataset = '/BulkGraviton_hh_GF_HH_14TeV_TuneCP5_pythia8/mkomm-ML_210512-d1606bbcd4ad268d18a3191de15a9732/USER'
    config.Data.outputDatasetTag = 'DeepNtuplesAK4-v00'
    config.Data.publication = False
    config.Data.unitsPerJob = 10
    config.Data.inputDBS = 'phys03'
    config.Data.splitting = 'FileBased'
    config.Data.allowNonValidInputDataset = True
    config.Data.outLFNDirBase = '/store/group/ml/Tagging4ScoutingHackathon/Adelina/DeepNtuples/12_3_0/ScoutingAK4-v00'
    config.section_('Site')
    config.Site.storageSite = 'T2_CH_CERN'
    config.section_('User')
    config.section_('Debug')
    ```

### Running the ntupiliser

Either your samples are in a CMSSW release >= 12_2_0 or you just succesfully re-ran the HLT step on your samples using the instructions above. We are now ready to create the ntuples. Make sure you have prepared your CMSSW release, then clone this repository and compile.

```
cd CMSSW_12_3_0/src
cmsenv

git clone git@github.com:alintulu/Run3ScoutingAnalysisTools.git -b ak4jet_gennano
scram b
```

1. Next you can either test the ntupiliser locally by running over one file 

```
INPUT=<INSERT YOUR SCOUTING SAMPLE>
OUTPUT=<INSERT YOUR OWN PATH>
cmsRun Run3ScoutingAnalysisTools/Analysis/test/ScoutingNanoAOD_NanoGEN_cfg.py inputFiles=$INPUT outputFile=$OUTPUT
```
2. Or run over all files with CRAB

Here is an example of a CRAB config file. You need to edit `outLFNDirBase` and `storageSite` to match the site you have access to write to. You can find more information about CRAB [here](https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile). `psetName` is the relative path to the ntupiliser config file. `userInputFiles` is a text file with your Scouting files listen one file per line.

```
from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferLogs = False
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.General.requestName = 'BulkGraviton_hh_GF_HH_14TeV_TuneCP5_pythia8'
config.section_('JobType')
config.JobType.numCores = 1
config.JobType.sendExternalFolder = True
config.JobType.pyCfgParams = ['outputFile=nano.root']
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.psetName = '../test/ScoutingNanoAOD_NanoGEN_cfg.py'
config.JobType.maxMemoryMB = 2000
config.section_('Data')
config.Data.userInputFiles = open('files_BulkGraviton.txt').readlines()
config.Data.outputPrimaryDataset = 'BulkGraviton_hh_GF_HH_14TeV_TuneCP5_pythia8'
config.Data.outputDatasetTag = 'BulkGraviton_hh_GF_HH_14TeV_TuneCP5_pythia8'
config.Data.publication = False
config.Data.unitsPerJob = 10
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.allowNonValidInputDataset = True
config.Data.outLFNDirBase = '/store/group/ml/Tagging4ScoutingHackathon/Adelina/DeepNtuples/12_3_0/ScoutingAK4-v00/Ntuples'
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.section_('User')
config.section_('Debug')
```

Ntuples successfully created!

## Training the network

We are now ready to start training the network. But first we need to download the code and install all dependencies on lxplus.
