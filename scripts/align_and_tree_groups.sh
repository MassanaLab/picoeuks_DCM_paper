#!/bin/bash

#SBATCH --job-name=align_and_tree
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=data/logs/align_and_tree_%A_%a.out
#SBATCH --error=data/logs/align_and_tree_%A_%a.err
#SBATCH --array=1-23%16

module load mafft
module load iqtree3

GROUP=$(awk "NR == ${SLURM_ARRAY_TASK_ID}" data/trees/groups.txt)
QUERY="data/trees/fasta/${GROUP}.fasta"
ALIGNMENT=${QUERY/.fasta/.pir}
OUT_PREFIX=data/trees/${GROUP}

## mafft

mafft --maxiterate 1000 --localpair --thread ${SLURM_CPUS_PER_TASK} ${QUERY} > ${ALIGNMENT}

## iqtree (model in auto mode)

iqtree3 \
 -bb 1000 \
 -s ${ALIGNMENT} \
 -alrt 1000 \
 -T ${SLURM_CPUS_PER_TASK} \
 -pre ${OUT_PREFIX}
