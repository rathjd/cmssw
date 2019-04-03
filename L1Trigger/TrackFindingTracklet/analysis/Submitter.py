import os,sys

directory = 'Batch3/'
os.system("mkdir "+directory)
with open('files3.txt') as f:
	lines = f.readlines()
	lines = [x.strip() for x in lines]
	jobCounter = 0
	for line in lines:
		print(line+" for job "+str(jobCounter))
		filename=directory+"job"+str(jobCounter)+".slrm"
		file = open(filename,"w")
		file.write("#!/bin/bash\n")
		file.write("#SBATCH -J TraMuNt\n")
		file.write("#SBATCH -p stakeholder-4g\n")
		file.write("#SBATCH --time=4:00:00\n")
		file.write("#SBATCH -n1\n")
		file.write("#SBATCH --mem-per-cpu=4000\n")
		file.write("#SBATCH -o log_TraMuNt_-%j.out\n")
		file.write("#SBATCH -e log_TraMuNt_-%j.err\n")
		file.write("export X509_USER_PROXY=$HOME/x509up_u20083\n");
		file.write("echo \"starting at `date` on `hostname`\"\n")
		file.write("echo \"SLURM_JOBID=$SLURM_JOBID\"\n")
		command="cmsRun /fdata/hepx/store/user/rathjd/IntegratedStubTruthTests/CMSSW_10_5_0/src/L1Trigger/TrackFindingTracklet/test/L1TrackNtupleMaker_cfg.py inputFile="+line+" outputFile="+directory+"output_"+str(jobCounter)+".root GEOMETRY=D35"
		file.write(command+"\n")
		file.write("echo \"ended at `date` on `hostname`\"\n")
		file.write("exit 0")
		file.close()
		os.system("sbatch "+directory+"job"+str(jobCounter)+".slrm")
		jobCounter+=1
