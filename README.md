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

5. Create a NanoAODGEN file that also contains some scouting collections

```
source nanogen.sh
```

Done!
