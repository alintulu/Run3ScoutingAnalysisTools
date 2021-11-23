#!/usr/bin/env python

import os, sys
from optparse import OptionParser

#get inputs
from optparse import OptionParser
parser=OptionParser()
parser.add_option("-q","--flavour",dest="jobFlavour",type="str",default="longlunch",help="job FLAVOUR",metavar="FLAVOUR")
parser.add_option("-p","--proxy",dest="proxyPath",type="str",default="noproxy",help="Proxy path")
parser.add_option("-f","--inputFile",dest="inputFile",type="str",help="Input file")
parser.add_option("-c","--cmssw",dest="cmssw",type="str",help="CMSSW path")
parser.add_option("-n","--maxEvents",dest="maxEvents",type="int",default=100,help="Max number of events per file")
opts, args = parser.parse_args()

job_dir = os.getcwd()+"/Jobs/"

with open(opts.inputFile) as file:
    for i, line in enumerate(file):
        input_file = line.rstrip()
        input_folder, input_name = input_file.split('/')[-2:]
        if "QCD_Pt15to7000" in input_file:
            sample = "QCD_Pt15to7000"
            isQCD = True
        elif "QCD_Pt300to7000_" in input_file:
            sample = "QCD_Pt300to7000"
            isQCD = True
        elif "TTbar" in input_file:
            sample= "TTbar"
            isQCD = False
        elif "BulkGraviton" in input_file:
            sample = "BulkGraviton"
            isQCD = False
        else:
            sys.exit("Does not recognise input file")
        output_edm = "/eos/user/a/adlintul/scouting/particlenet/particle_features/reHLT/prod/edm/" + sample + "/" + input_folder
        os.system('mkdir -p %s'%output_edm)
        os.system('mkdir -p %s'%output_edm.replace("edm", "nano_msd"))
        output_edm += "/" + input_name
        if (opts.maxEvents != -1):
           output_edm_num = output_edm[:-5]+"_numEvent%d"%opts.maxEvents+".root"
        else:
           output_edm_num = output_edm
        if (opts.maxEvents != -1):
           output_nano = output_edm.replace("edm", "nano_msd")[:-5]+"_numEvent%d"%opts.maxEvents+".root"
        else:
           output_nano = output_edm.replace("edm", "nano_msd")
        os.system('mkdir -p %s'%job_dir+sample+"/"+input_folder+"/Job_%d"%i)

        job_name=sample+"/"+input_folder+"/Job_%d"%i+"/sub_%s.sh"%(str(i))
        job=open(job_dir+"/"+job_name,'w')
        job.write("#!/bin/sh\n\n")
        job.write("INPUT={0}\n".format(input_file))
        #job.write("OUTPUT_EDM={0}\n".format(output_edm))
        #job.write("OUTPUT_EDM_NUM={0}\n".format(output_edm_num))
        job.write("OUTPUT_NANO={0}\n\n".format(output_nano))
        if opts.proxyPath != "noproxy":
            job.write("export X509_USER_PROXY=$1\n")
            job.write("voms-proxy-info -all\n")
            job.write("voms-proxy-info -all -file $1\n")
        job.write("ulimit -v 5000000\n")
        job.write("cd $TMPDIR\n")
        job.write("mkdir Job_%s\n"%str(i))
        job.write("cd Job_%s\n"%str(i))
        job.write("cd %s\n"%(opts.cmssw))
        job.write("eval `scramv1 runtime -sh`\n")
        job.write("cd -\n")
        job.write("cmsRun {reHLT_file} inputFiles=file:$INPUT outputFile=$OUTPUT_NANO maxEvents={maxEvents} isQCD={isQCD}".format(reHLT_file=opts.cmssw+"/Run3ScoutingAnalysisTools/Analysis/test/reHLT_2in1.py", maxEvents=opts.maxEvents, isQCD=isQCD))
        #job.write("cmsRun {reHLT_file} inputFiles=file:$INPUT outputFile=$OUTPUT_EDM maxEvents={maxEvents}\n".format(reHLT_file=opts.cmssw+"/reHLT.py", maxEvents=opts.maxEvents)) #reHLT
        #job.write("cmsRun {scoutingNano} inputFiles=file:$OUTPUT_EDM_NUM outputFile=$OUTPUT_NANO isQCD={isQCD} isMC=True useWeights=False GlobalTagMC=112X_mcRun3_2021_realistic_v16\n".format(scoutingNano=opts.cmssw+"/Run3ScoutingAnalysisTools/Analysis/test/ScoutingNanoAOD_cfg.py", isQCD=isQCD)) #ntuple creation
        #job.write("rm $OUTPUT_EDM_NUM\n")
        job.close()
        os.system("chmod +x %s"%(job_dir+job_name))

condor_str = "executable = $(filename)\n"
if opts.proxyPath != "noproxy":
    condor_str += "Proxy_path = %s\n"%opts.proxyPath
    condor_str += "arguments = $(Proxy_path) $Fp(filename) $(ClusterID) $(ProcId)\n"
else:
    condor_str += "arguments = $Fp(filename) $(ClusterID) $(ProcId)\n"
condor_str += "output = $Fp(filename)hlt.stdout\n"
condor_str += "error = $Fp(filename)hlt.stderr\n"
condor_str += "log = $Fp(filename)hlt.log\n"
condor_str += '+JobFlavour = "%s"\n'%opts.jobFlavour
condor_str += '+AccountingGroup = "group_u_CMST3.all"\n'
condor_str += "queue filename matching ("+job_dir+"*/*/Job_*/*.sh)"
condor_name = os.getcwd()+"/condor.sub"
condor_file = open(condor_name, "w")
condor_file.write(condor_str)
