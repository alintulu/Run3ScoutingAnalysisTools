This repository creates jet-level ntuples for training flavour jet tagging with ParticleNet (originally created for the "Jet tagging using the Run3 Scouting data" hackathon 8th-11th Nov 2021).

# Start here

Every step described in this README has be tested on lxplus only. There are three main steps to the full worflow: 

1. [Creating the ntuples](#creating-the-ntuples)
2. [Training the network](#training-the-network)
3. [Checking the network performance](#checking-the-network-performance)

Next I'll try to describe each step as carefully as possible.

## Creating the ntuples

The first step is to create the jet-level ntuples which are used for training and testing the network. If your Scouting samples are produced in the 12_2_0 or a newer CMSSW release you can skip to [Running the ntupiliser](#running-the-ntupiliser). If this is not the case please continue reading here.

These are the samples we used:
- BulkGraviton to HH to 4Q (`/BulkGraviton_hh_GF_HH_14TeV_TuneCP5_pythia8/mkomm-ML_210512-d1606bbcd4ad268d18a3191de15a9732/USER`)
- QCD (`/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/mkomm-ML_210512-d1606bbcd4ad268d18a3191de15a9732/USER`)
- TTbar (`/TTbar_TuneCP5_14TeV-pythia8/mkomm-ML_210512-d1606bbcd4ad268d18a3191de15a9732/USER`)

These samples are created in the 11_2_1_Patatrack CMSSW release and are in the format of RAWMINIAODSIM. Since they are RAW we are able to re-run the HLT step. This is necessary since the 11_2_1_Patatrack Scouting format does not include all necessary information. The correct Scouting format is available in the 12_2_0 release and newer releases. Next I will describe how we reconstruct the samples above in the 12_3_0 release.

### Re-running the HLT step

1. Prepare the CMSSW release on lxplus.

```
cmsrel CMSSW_12_3_0
cd CMSSW_12_3_0/src
cmsenv
```

2. Create a file called `reHLT.sh` at `$CMSSW_BASE/src`.

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

3. Create the HLT config file by running `reHLT.sh`.

```
INPUT=/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/Scouting/Run3/ML_210512/GluGluHToBB_M125_masseffects_14TeV_TuneCP5_powheg_pythia8/ML_210512/210602_090726/0000/scouting_75.root
OUTPUT=<INSERT YOUR OUTPUT FILE PATH>
source reHLT.sh $INPUT $OUTPUT reHLT.py
```

4. Re-run the HLT step by running `reHLT.py`. Here you are two options.

    a. Run the file locally over one file.

    ```
    cmsRun reHLT.py
    ```

    b. Run over all samples using CRAB.

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

Either your samples are already in a CMSSW release >= 12_2_0 or you just succesfully re-ran the HLT step on your samples using the instructions above. We are now ready to create the ntuples. Make sure you have prepared your CMSSW release, then clone this repository and compile.

```
cd CMSSW_12_3_0/src
cmsenv

git clone git@github.com:alintulu/Run3ScoutingAnalysisTools.git -b ak4jet_gennano
scram b
```

1. Next you can either test the ntupiliser locally by running over one file.

```
INPUT=<INSERT YOUR SCOUTING FILE PATH>
OUTPUT=<INSERT YOUR OUTPUT FILE PATH>
cmsRun Run3ScoutingAnalysisTools/Analysis/test/ScoutingNanoAOD_NanoGEN_cfg.py inputFiles=$INPUT outputFile=$OUTPUT
```
2. Or run over all files with CRAB.

Here is an example of a CRAB config file. You need to edit `outLFNDirBase` and `storageSite` to match the site you have access to write to. You can find more information about CRAB [here](https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile). `psetName` is the relative path to the ntupiliser config file. `userInputFiles` is a text file with your Scouting files listed one file per line.

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

## Training the network

Hopefully you have now succesfully created your own training ntuples. We are then ready to start training the network. But first we need to install Miniconda, Pytorch (cuda10.2) and Weaver on a normal lxplus node. 

### Installing the network

Log into lxplus and follow these [instructions](https://github.com/hqucms/weaver#set-up-your-environment) to set up the necessary environment for training.

With everything succesfully installed, we can start training.

### Accessing GPU

If you do not have your own GPU resources, you can access GPU using the CERN HTCondor batch system. To start training on HTCondor we need (1) the shell script and (2) the submission file.

Here is an example of the shell script. I name mine `run.sh`.
```
#!/bin/bash

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/afs/cern.ch/work/a/adlintul/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/afs/cern.ch/work/a/adlintul/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/afs/cern.ch/work/a/adlintul/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/afs/cern.ch/work/a/adlintul/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
conda activate py3

echo "Setting up CUDA..."
export PATH="/usr/local/cuda/bin:$PATH"
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

nvcc --version
nvidia-smi
export TMPDIR=`pwd`

FILES="/eos/cms/store/group/ml/Tagging4ScoutingHackathon/Adelina/DeepNtuples/12_3_0/ScoutingAK4-v00/Ntuples_matching_fixed/TTbar_TuneCP5_14TeV-pythia8/TTbar_TuneCP5_14TeV-pythia8/220525_142744/0000/nano_*.root /eos/cms/store/group/ml/Tagging4ScoutingHackathon/Adelina/DeepNtuples/12_3_0/ScoutingAK4-v00/Ntuples_matching_fixed/BulkGraviton_hh_GF_HH_14TeV_TuneCP5_pythia8/BulkGraviton_hh_GF_HH_14TeV_TuneCP5_pythia8/220525_142803/0000/nano_*.root /eos/cms/store/group/ml/Tagging4ScoutingHackathon/Adelina/DeepNtuples/12_3_0/ScoutingAK4-v00/Ntuples_matching_fixed/QCD_Pt15to7000_TuneCP5_14TeV-pythia8/DeepNtuplesAK4-v00/220525_142827/0000/nano_*.root"

MODEL_PREFIX="/eos/user/a/adlintul/scouting/particlenet/kubeflow/outputs/train_pixelOnly_vs_Full-{auto}/net"
LOG="/eos/user/a/adlintul/scouting/particlenet/kubeflow/logs/train_pixelOnly_vs_Full-{auto}.log"
CONFIG="/afs/cern.ch/work/a/adlintul/public/condor-gpu/config.yaml"

WEAVER_PATH="/afs/cern.ch/work/a/adlintul/public/scouting/particlenet/weaver"
cd $WEAVER_PATH

wget https://raw.githubusercontent.com/hqucms/weaver-benchmark/main/networks/top/particlenet_pf.py

cp $CONFIG data_config.yaml

python train.py --data-train $FILES --data-test $FILES --data-config "data_config.yaml" --network-config "particlenet_pf.py" --model-prefix $MODEL_PREFIX --predict-output "pred.root" --num-workers "0" --fetch-step 0.01 --data-fraction "1" --batch-size "512" --num-epochs "5" --start-lr "5e-3" --optimizer "ranger" --log $LOG --gpus "0"
```

It looks chaotic but its use is very simple. First it initializes the conda environment. It sets a few necessary environment variables and then we `cd` into our weaver installation directory. When we are in the weaver installation directory we are able to run the `train.py` script and initiate the training. There are a couple of things you need to change in order to run the shell script succesfully in your environment.

1. Change the conda path in the "conda initialize" part of the script to match your path. Your conda path should be known from installing Miniconda in [Installing the network](#installing-the-network).
2. Change the `$FILES` variable to point to your training ntuples.
3. Change `$MODEL_PREFIX` and `$LOG` to point to a path you have write access to (note: by including `{auto}` here, weaver will create a unique output name for each training making sure that no output is overwritten by mistake).
4. Change the `$WEAVER_PATH` variable to point to your weaver installation directory.
5. Change the `$CONFIG` variable to point to your data config. You can find an example of such a config in [Data config](#data-config).


And here is an example of the submission file. I name mine `submit.jdl`. You do not need to change anything in this file as long as your shell script is called `run.sh`.

```
request_cpus          = 4
request_gpus          = 1
requirements          = regexp("V100", TARGET.CUDADeviceName)
executable            = run.sh
should_transfer_files = YES
output                = train.$(ClusterId).$(ProcId).out
error                 = train.$(ClusterId).$(ProcId).err
log                   = train.$(ClusterId).log
WhenToTransferOutput  = ON_EXIT
want_graceful_removal = true

+JobFlavour = "testmatch"

queue
```

### Data config

With the shell script and HTCondor submission file we almost have everything needed to train our network. The only thing missing is our data config.

Here is an example of a data config. You do not need to change anything in this file.

```
selection:
   ### use `&`, `|`, `~` for logical operations on numpy arrays
   ### can use functions from `math`, `np` (numpy), and `awkward` in the expression
   >-
   (jet_hadronFlavour>-90) & (np.abs(jet_eta)<2.5) & (jet_pt>15) & (jet_pt<600) & (event_no%7!=0)
   & ( (jet_hadronFlavour>0) | ((jet_hadronFlavour==0) & ( (np.abs(jet_partonFlavour)>0) & (np.abs(jet_partonFlavour)<4) | (jet_partonFlavour==21) ) ) )

test_time_selection:
   (jet_hadronFlavour>-90) & (np.abs(jet_eta)<2.5) & (jet_pt>15) & (jet_pt<600) & (event_no%7==0)
   & ( (jet_hadronFlavour>0) | ((jet_hadronFlavour==0) & ( (np.abs(jet_partonFlavour)>0) & (np.abs(jet_partonFlavour)<4) | (jet_partonFlavour==21) ) ) )

treename: full/tree

new_variables:
   ### [format] name: formula
   ### can use functions from `math`, `np` (numpy), and `awkward` in the expression
   pfcand_mask: awkward.JaggedArray.ones_like(pfcand_etarel)
   label_b:     (jet_nBHadrons==1)
   label_bb:    (jet_nBHadrons>1)
   label_c:     (jet_nBHadrons==0) & (jet_nCHadrons==1)
   label_cc:    (jet_nBHadrons==0) & (jet_nCHadrons>1)
   label_uds:   (jet_hadronFlavour==0) & (np.abs(jet_partonFlavour)>0) & (np.abs(jet_partonFlavour)<4)
   label_g:     (jet_hadronFlavour==0) & (jet_partonFlavour==21)
   label_catB:  (jet_hadronFlavour==5)
   label_catC:  (jet_hadronFlavour==4)
   label_undef: (jet_hadronFlavour==0) & (jet_partonFlavour!=21) & ((np.abs(jet_partonFlavour)>=4) | (jet_partonFlavour==0))


preprocess:
  ### method: [manual, auto] - whether to use manually specified parameters for variable standardization
  method: manual
  ### data_fraction: fraction of events to use when calculating the mean/scale for the standardization
  data_fraction:

inputs:
   pf_points:
      length: 50
      vars:
         - [pfcand_etarel, null]
         - [pfcand_phirel, null]
   pf_features:
      length: 50
      vars:
      ### [format 1]: var_name (no transformation)
      ### [format 2]: [var_name, 
      ###              subtract_by(optional, default=None, no transf. if preprocess.method=manual, auto transf. if preprocess.method=auto), 
      ###              multiply_by(optional, default=1), 
      ###              clip_min(optional, default=-5), 
      ###              clip_max(optional, default=5), 
      ###              pad_value(optional, default=0)]
         - [pfcand_pt_log_nopuppi, 0.8, 0.7]
         - [pfcand_e_log_nopuppi, 1.2, 0.6]
         - [pfcand_etarel, null]
         - [pfcand_phirel, null]
         - [pfcand_abseta, 0.8, 1.2]
         - [pfcand_charge, null]
         - pfcand_isEl
         - pfcand_isMu
         - pfcand_isGamma
         - pfcand_isChargedHad
         - pfcand_isNeutralHad
         - [pfcand_lostInnerHits, null]
         - [pfcand_normchi2, 5, 0.2]
         - [pfcand_quality, 0, 0.2]
         - [pfcand_dz, 0, 150]
         - [pfcand_dzsig, 0, 0.7]
         - [pfcand_dxy, 0, 270]
         - [pfcand_dxysig, 0, 1.1]
         - [pfcand_btagEtaRel, 2.1, 0.5]
         - [pfcand_btagPtRatio, 0, 1]
         - [pfcand_btagPParRatio, 0, 1]
   pf_mask:
      length: 50
      vars:
         - [pfcand_mask, null]

labels:
   ### type can be `simple`, `custom`
   ### [option 1] use `simple` for binary/multi-class classification, then `value` is a list of 0-1 labels
   type: simple
   value: [label_b, label_bb, label_c, label_cc, label_uds, label_g, label_undef]
   ### [option 2] otherwise use `custom` to define the label, then `value` is a map
   # type: custom
   # value:
      # target_mass: np.where(fj_isQCD, fj_genjet_sdmass, fj_gen_mass) 

observers:
   - event_no
   - jet_pt
   - jet_eta
   - jet_phi
   - jet_mass
   - jet_hadronFlavour
   - jet_partonFlavour
   - jet_nBHadrons
   - jet_nCHadrons

monitor_variables:
   # - jet_pt
   # - jet_eta
   # - jet_hadronFlavour
   # - jet_partonFlavour
   # - jet_nBHadrons
   # - jet_nCHadrons
   # - jet_genjet_pt

weights:
   ### [option 1] use precomputed weights stored in the input files
   # use_precomputed_weights: true
   # weight_branches: [weight, class_weight]
   ### [option 2] compute weights on-the-fly using reweighting histograms
   use_precomputed_weights: false
   reweight_method: flat
   reweight_vars:
      # np.round(np.exp(np.linspace(np.log(15), np.log(1000), 9))).astype('int').tolist()
      # jet_pt: [15, 25, 43, 72, 122, 207, 350, 592, 1000]
      # np.round(np.exp(np.linspace(np.log(15), np.log(1000), 16))).astype('int').tolist()
      jet_pt: [15, 20, 26, 35, 46, 61, 80, 106, 141, 186, 247, 326, 432, 600]
      jet_eta: [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
   reweight_classes: [label_catB, label_catC, label_uds, label_g]
   class_weights: [1, 1, 2.5, 5]
   reweight_threshold: 15
   reweight_hists:
```

### Submitting the training

We are now ready to submit the training to HTCondor.

```
condor_submit submit.jdl
```

In my experience it takes about 2 days to complete the training. When the training is complete you will find a `log`, `out` and `error` file in your directory. Look through them carefully to make sure everything ran successfully. At the bottom of the `out` file you'll find a path to a file containing the prediction. We will need this file in the next step.

## Checking the network performance

It is now time to check the network performance. We have already automatically performed a test after training the network and the a predicrion file has been created. We can use this jupyer notebook to create ROC curves.
