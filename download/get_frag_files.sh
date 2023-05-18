#!/bin/bash
#SBATCH --job-name=get_encode_adrenal_frags  ## Name of the job
#SBATCH --partition carter-compute              ## partition/queue name
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=atacdata-%J.out ## output log file
#SBATCH --error=atacdata-%J.err ## error log file
#SBATCH --mem=128G

mkdir fragments/
cd fragments/

xargs -L 1 curl -O -J -L < ../adrenal_10x_fragment_files.txt

for f in *.tar.gz
do
	mkdir "${f%.tar.gz}"
	tar xf "$f" -C "${f%.tar.gz}"
	rm "$f"
done
