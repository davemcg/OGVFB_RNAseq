#!/bin/bash

# run in the data folder for this project
# on biowulf2:
# /data/mcgaugheyd/projects/nei/brooks/oca_rna-seq

source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh
mamba activate base

mkdir -p 00log

#module load snakemake || exit 1

sbcmd="sbatch --cpus-per-task={threads} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"


snakemake -s /home/mcgaugheyd/git/OGVFB_RNAseq/Snakefile_EiaD \
  --rerun-triggers mtime \
  -pr --local-cores 2 --jobs 500 \
  --cluster-config /home/mcgaugheyd/git/OGVFB_RNAseq/cluster.json \
  --cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
  --configfile $1 --use-conda \
  -k --restart-times 0 \
  --resources parallel=4 
