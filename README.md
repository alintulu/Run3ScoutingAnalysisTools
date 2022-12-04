## Download, compile and run code

The code was develop in `CMSSW_12_4_5`. It needs a custom branch of this release.

### Instructions

1. Log into lxplus
2. Prepare the CMSSW release

```
cmsrel CMSSW_12_4_5
cd CMSSW_12_4_5/src
cmsenv
git cms-init
git cms-addpkg CommonTools PhysicsTools RecoBTag
git remote add alintulu git@github.com:alintulu/cmssw.git
git fetch alintulu scoutingNanoAOD
git checkout --track alintulu/scoutingNanoAOD
```

3. Clone this repository and compile

```
git clone git@github.com:alintulu/Run3ScoutingAnalysisTools.git -b nanoaod
scram b -j8
```
4. Create a custom scouting NanoAOD

Data
    
```
cmsRun scouting_data.py
```

MC
    
```
cmsRun scouting_mc.py
```

5. Inspect your file

```
$ root -b -l scoutingnano_data.root
root> Events->Print()
*............................................................................*
*Br  180 :ScoutingMuon_trk_dxyError : Float_t trk_dxyError                   *
*Entries :        1 : Total  Size=        768 bytes  File Size  =        110 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br  181 :ScoutingMuon_trk_dxy_dsz_cov : Float_t trk_dxy_dsz_cov             *
*Entries :        1 : Total  Size=        783 bytes  File Size  =        113 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br  182 :ScoutingMuon_trk_dz : Float_t trk_dz                               *
*Entries :        1 : Total  Size=        738 bytes  File Size  =        104 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br  183 :ScoutingMuon_trk_dzError : Float_t trk_dzError                     *
*Entries :        1 : Total  Size=        763 bytes  File Size  =        109 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br  184 :ScoutingMuon_trk_eta : Float_t trk_eta                             *
*Entries :        1 : Total  Size=        743 bytes  File Size  =        105 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
.
.
.
```

### AOD to MINIAOD with scouting

```
cmsDriver.py \
  --python_filename aod2miniaod.py \
  --eventcontent MINIAODSIM \
  --customise Configuration/DataProcessing/Utils.addMonitoring \
  --datatier MINIAODSIM \
  --fileout file:miniaod.root \
  --conditions auto:phase1_2022_realistic \
  --step PAT \
  --geometry DB:Extended \
  --filein "/store/mc/Run3Summer22DRPremix/MinBias_TuneCP5_14TeV-pythia8/AODSIM/Pilot_124X_mcRun3_2022_realistic_v11-v2/50000/b255101d-13ab-4068-a46c-a5ea9e1a8a6c.root" \
  --era Run3 \
  --no_exec \
  --mc \
  -n 1 \
 --customise_commands "process.MINIAODSIMoutput.outputCommands.extend(['drop *_hlt*_*_*','keep *_hltScoutingEgammaPacker_*_*','keep *_hltScoutingMuonPacker_*_*','keep *_hltScoutingPFPacker_*_*','keep *_hltScoutingPrimaryVertexPacker_*_*','keep *_hltScoutingTrackPacker_*_*','drop *_*_*_PAT','keep *_*GenParticles_*_*','keep *_slimmedGenJets*_*_*', 'keep *_slimmedMET*_*_*', 'keep edmTriggerResults_*_*_PAT'])"
```

### Left to figure out

- How to get number of daughters for the jet
