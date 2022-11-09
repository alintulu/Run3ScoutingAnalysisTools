## Download, compile and run code

The code was develop in `CMSSW_12_4_5`.

### Instructions

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
4. Create a NanoAOD file that also contains some scouting collections

```
source nano.sh
```

You can find an example of such a file in lxplus at `~adlintul/public/scouting/nanoaod`

5. Create a NanoAODGEN file that also contains some scouting collections

```
source nanogen.sh
```

You can find an example of such a file in lxplus at `~adlintul/public/scouting/nanoaod`

E.g.

```
$ root -b -l ~adlintul/public/scouting/nanoaod/nanogen.root
root> Events->Print()
*............................................................................*
*Br   68 :Run3ScoutingJet_phi : Float_t phi                                  *
*Entries :        1 : Total  Size=        812 bytes  File Size  =        172 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br   69 :Run3ScoutingJet_pt : Float_t pt                                    *
*Entries :        1 : Total  Size=        807 bytes  File Size  =        171 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br   70 :Run3ScoutingJet_nCh : Int_t number of charged hadrons in the jet   *
*Entries :        1 : Total  Size=        838 bytes  File Size  =        172 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br   71 :Run3ScoutingJet_nElectrons : Int_t number of electrons in the jet  *
*Entries :        1 : Total  Size=        860 bytes  File Size  =        172 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.04     *
*............................................................................*
*Br   72 :Run3ScoutingJet_nMuons : Int_t number of muons in the jet          *
*Entries :        1 : Total  Size=        840 bytes  File Size  =        168 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.04     *
*............................................................................*
*
*
* etc
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
 --customise_commands "process.MINIAODSIMoutput.outputCommands.extend(['drop *_hlt*_*_*','keep *_hltScoutingEgammaPacker_*_*','keep *_hltScoutingMuonPacker_*_*','keep *_hltScoutingPFPacker_*_*','keep *_hltScoutingPrimaryVertexPacker_*_*','keep *_hltScoutingTrackPacker_*_*','drop *_*_*_PAT','keep *_*GenParticles_*_*','keep *_slimmedGenJets*_*_*', 'keep *_slimmedMET*_*_*'])"
```

### Left to figure out

- How to get softdrop mass from AK8 jets
- How to use bDiscriminator with PF (not PAT) jets
- Why PF jets does not take mass and area
- How to add n2b1 to ak8 jets
