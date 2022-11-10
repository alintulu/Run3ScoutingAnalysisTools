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
git clone git@github.com:alintulu/Run3ScoutingAnalysisTools.git -b nanoaod_simple
scram b -j 8
```
4. Create a NanoAOD file

```
cmsRun scouting_mc.py
```

5. Check the contents of the file you just produced

```
$ root -b -l scoutingnano_mc.root
root> Events->Print()
*............................................................................*
*Br   68 :ScoutingJet_phi : Float_t phi                                  *
*Entries :        1 : Total  Size=        812 bytes  File Size  =        172 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br   69 :ScoutingJet_pt : Float_t pt                                    *
*Entries :        1 : Total  Size=        807 bytes  File Size  =        171 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br   70 :ScoutingJet_nCh : Int_t number of charged hadrons in the jet   *
*Entries :        1 : Total  Size=        838 bytes  File Size  =        172 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
*............................................................................*
*Br   71 :ScoutingJet_nElectrons : Int_t number of electrons in the jet  *
*Entries :        1 : Total  Size=        860 bytes  File Size  =        172 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.04     *
*............................................................................*
*Br   72 :ScoutingJet_nMuons : Int_t number of muons in the jet          *
*Entries :        1 : Total  Size=        840 bytes  File Size  =        168 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.04     *
*............................................................................*
*
*
* etc
```
