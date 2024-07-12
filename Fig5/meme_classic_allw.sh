#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBATCH --mem=64g


meme -oc $1-meme-classic-allw -objfun classic -bfile unique_lncRNA_500_withfullairn_markov -mod anr -maxw 8 -minw 4 -minsites 2 -maxsites 1000 -nmotifs 50 -rna -allw $1.fasta
