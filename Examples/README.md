# How to work with the customised scouting NanoAOD format

1. [Coffea](#Coffea)
2. RDataFrame (TODO)
3. NanoAODTools (TODO)

## Coffea

The plan is to add the ScoutingNanoAOD schema (read more about schemas [here](https://coffeateam.github.io/coffea/notebooks/nanoevents.html)) into the official coffea release. However until then, one must do the following to use this schema


1. Log into lxplus or lpc, then

```
git clone git@github.com:alintulu/coffea.git -b scouting
```

2. Start the coffea singularity image

```
singularity shell -B ${PWD}:/work /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-dask:latest
cd /work/coffea
pip install .
```
Done!

If you want to work in a jupyter notebook, remeber to log into lxplus/lpc with

```
ssh -L 8766:localhost:8766 -Y <username>@lxplus.cern.ch
```

And after doing all the steps above, do

```
jupyter notebook --no-browser --port=8766
```
