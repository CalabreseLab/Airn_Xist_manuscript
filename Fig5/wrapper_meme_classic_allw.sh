#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBATCH --mem=32g

module load meme

for file in *.fasta
do
  base=${file%.fasta}
  sbatch --export=ALL,base=${base} meme_classic_allw.sh ${base}
done
