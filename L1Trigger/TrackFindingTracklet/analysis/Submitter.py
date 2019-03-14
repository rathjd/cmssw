import os,sys

with open('filelist.txt') as f:
	lines = f.readlines()
	lines = [x.strip() for x in lines]
	jobCounter = 0
	for line in lines:
		print(line+" for job "+str(jobCounter))
		filename="job"+str(jobCounter)+".slrm"
		file = open(filename,"w")
		file.write("#!/bin/bash\n")
		file.write("#SBATCH -J TraMuNt\n")
		file.write("#SBATCH -p stakeholder-4g\n")
		file.write("#SBATCH --time=4:00:00\n")
		file.write("#SBATCH -n1\n")
		file.write("#SBATCH --mem-per-cpu=4000\n")
		file.write("#SBATCH -o log_TraMuNt_-%j.out\n")
		file.write("#SBATCH -e log_TraMuNt_-%j.err\n")
		file.write("echo \"starting at `date` on `hostname`\"\n")
		file.write("echo \"SLURM_JOBID=$SLURM_JOBID\"\n")
		command="cmsRun /fdata/hepx/store/user/rathjd/TrackletMuonFusion/CMSSW_9_3_7/src/L1Trigger/TrackFindingTracklet/test/L1TrackNtupleMaker_cfg.py inputFile=file:/fdata/hepx/store/mc/PhaseIIFall17D/SingleMu_FlatPt-2to100/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/00000/"+line+" outputFile=output_"+line
		file.write(command+"\n")
		file.write("echo \"ended at `date` on `hostname`\"\n")
		file.write("exit 0")
		file.close()
		os.system("sbatch job"+str(jobCounter)+".slrm")
		jobCounter+=1
