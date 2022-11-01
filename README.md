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

Done!
