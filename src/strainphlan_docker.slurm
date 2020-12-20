#!/bin/sh
#SBATCH --ntasks=36
#SBATCH --partition=batch
#SBATCH --time=7-0:0
#SBATCH --mem=500G
#SBATCH --job-name=tree
#SBATCH --error=job.tree.%J.err
#SBATCH --output=job.tree.%J.out

module load singularity
unset REGISTRY

singularity exec -B metaphlan_databases:/opt/conda/bin/metaphlan_databases,zempleni:/tmp docker://biobakery/strainphlan:2.7.7 sh $1
