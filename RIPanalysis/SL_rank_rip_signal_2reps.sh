#!/bin/bash

#SBATCH -J ripRank_%j
#SBATCH -p general
#SBATCH --mem=100G
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=3:00:00
#SBATCH -o ripRank_%j.out
#SBATCH -e ripRank_%j.err

#### Goal: rank transcript based on binding signal to a RBP

# Example command:
    # sbatch 22sp_rank_rip_signal_proj.sh $1
    # $1 = rbp name (all lower cases) e.g., alyref

module load kallisto
module load samtools

# input files needed
    # rip files:
        # /proj/calabrlb/users/Zhiyue/22_sp/rip/*$1*
    # alignment of rip:
        # /proj/calabrlb/users/Zhiyue/22_sp/rip_sam/$1_Aligned_filteredsq30.out.sam
    # bam alignment of igg (sorted and indexed):
        # /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/igg_sorted_neg.bam
        # /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/igg_sorted_pos.bam
        # /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/igg_sorted_neg.bam.bai
        # /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/igg_sorted_pos.bam.bai
    # rip peaks (enriched over 2igg in => 2 replicates from vet_peaks pipeline):
        # /proj/calabrlb/users/Zhiyue/22_sp/rip_peaks_2igg/$1_peaks_2igg_2reps.bed
    # kallisto reference index:
        # /work/users/s/h/shuang9/rip/mm10_kallisto_transcriptome.index
    
    
# do following work in $1_kallisto
mkdir $1_kallisto
cd $1_kallisto

# move input here
mkdir rip_fastq
cp /work/users/s/h/shuang9/rip/rip_all_data/*$1* rip_fastq
gunzip rip_fastq/*.fastq.gz
cp /work/users/s/h/shuang9/rip/rip_all_stats/$1_classify/$1_Aligned_filteredsq30.out.sam .

cp /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/igg_sorted_neg.bam .
cp /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/igg_sorted_pos.bam .
cp /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/igg_sorted_neg.bam.bai .
cp /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/igg_sorted_pos.bam.bai .

################## need to change this ##############
cp /work/users/s/h/shuang9/rip/rip_all_stats/$1_classify/$1_peaks_2igg_2reps.bed .

### 1. find aligned reads that overlap with enriched peaks
    # Note: peak file is $1_peaks_2igg.bed

# split peak strands for strand-specific filtering
cat $1_peaks_2igg_2reps.bed | awk '{if ($6=="+") print $0}' > $1_peaks_2igg_pos.bed
cat $1_peaks_2igg_2reps.bed | awk '{if ($6=="-") print $0}' > $1_peaks_2igg_neg.bed

# split sam by strand
samtools view -h -F 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_neg.out.sam
samtools view -h -f 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_pos.out.sam

## convert SAM to BAM
    # -S: take in SAM (by default it expects BAM)
    # -b: output is BAM (by default it produces BAM) 
samtools view -S -b $1_Aligned_neg.out.sam > $1_neg.bam
samtools view -S -b $1_Aligned_pos.out.sam > $1_pos.bam

## sort and index BAM
samtools sort $1_neg.bam -o $1_sorted_neg.bam
samtools sort $1_pos.bam -o $1_sorted_pos.bam
samtools index $1_sorted_neg.bam
samtools index $1_sorted_pos.bam

# save files to proj: bam
# cp $1_sorted_neg.bam /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/
# cp $1_sorted_pos.bam /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/
# cp $1_sorted_neg.bam.bai /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/
# cp $1_sorted_pos.bam.bai /proj/calabrlb/users/Zhiyue/22_sp/rip_bam/

## filtering: extract SAM entries overlapping peaks.bed from both rip and igg
    # -L: only output alignments overlapping the BED FILE.
    # -h: keep header
samtools view $1_sorted_neg.bam -h -L $1_peaks_2igg_neg.bed > $1_peaks_neg.out.sam
samtools view $1_sorted_pos.bam -h -L $1_peaks_2igg_pos.bed > $1_peaks_pos.out.sam

samtools view igg_sorted_neg.bam -h -L $1_peaks_2igg_neg.bed > igg_$1_peaks_neg.out.sam
samtools view igg_sorted_pos.bam -h -L $1_peaks_2igg_pos.bed > igg_$1_peaks_pos.out.sam

## To convert peaks.sam and igg_peaks.sam back to FASTQ, must do bam processing first
# convert SAM to BAM
samtools view -S -b $1_peaks_neg.out.sam > $1_peaks_neg.bam
samtools view -S -b $1_peaks_pos.out.sam > $1_peaks_pos.bam

samtools view -S -b igg_$1_peaks_neg.out.sam > igg_$1_peaks_neg.bam
samtools view -S -b igg_$1_peaks_pos.out.sam > igg_$1_peaks_pos.bam

# convert BAM to FASTQ
samtools fastq $1_peaks_neg.bam > $1_peaks.fastq
samtools fastq $1_peaks_pos.bam >> $1_peaks.fastq

samtools fastq igg_$1_peaks_neg.bam > igg_$1_peaks.fastq
samtools fastq igg_$1_peaks_pos.bam >> igg_$1_peaks.fastq

### 2. Kallisto: align both rip and igg peak.fastq to transcriptome

# var: kallisto index (input file)
    # derived from gencode.vM25.basic.annotation.complete.ERCC.fa:
    # sbatch --mem 250g -N 1 -n 16 -t 3:00:00 -J kallisto_idx -o kallisto_idx.out -e kallisto_idx.err --wrap='kallisto index -i gencode.vM25.basic.annotation.complete.ERCC.idx gencode.vM25.basic.annotation.complete.ERCC.fa'
kallisto_idx='/work/users/s/h/shuang9/rip/mm10_kallisto_transcriptome.index'

## kallisto alignment
    # -i: location of the transcriptome index
    # -o: kallisto output folder for each dataset
    # --single -l 200 -s 50: single-end data
    # -l and -s are required for SE reads: estimated average and standard deviation of fragment length
    # rf-stranded: reverse stranded
kallisto quant -i $kallisto_idx -o $1_kallisto --single -l 200 -s 50 --rf-stranded $1_peaks.fastq
kallisto quant -i $kallisto_idx -o igg_$1_kallisto --single -l 200 -s 50 --rf-stranded igg_$1_peaks.fastq

## Make table summarizing estimated counts
cat $1_kallisto/abundance.tsv | cut -f1,2,3,4 | sed 1d > $1_kallisto_id_est_cts.txt
cat igg_$1_kallisto/abundance.tsv | cut -f4 | sed 1d > igg_$1_kallisto_est_cts.txt
paste $1_kallisto_id_est_cts.txt igg_$1_kallisto_est_cts.txt > $1_kallisto_igg.txt

## Normalization: convert estimated counts to rpm
# calculate rip reads
fastq_lin_num=$(wc -l rip_fastq/*$1*.fastq | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
(( reads = fastq_lin_num / 4 ))

# col6 = rip rpm; col7 = igg rpm; col8 = rip_rpm - igg_rpm = 6 - 7
cat $1_kallisto_igg.txt | awk -v OFS="\t" -v reads=$reads -v reads_igg=124344953 '{print $0, (($4*1000000)/reads), (($5*1000000)/reads_igg)}' | awk -v OFS="\t" '{print $0, $6-$7}' > $1_kallisto_igg_rpm.txt
sed -i "1i gene_ID\tlength\teff_length\t$1_est_cts\tigg_est_cts\t$1_rpm\tigg_rpm\t$1_rpm_over_igg" $1_kallisto_igg_rpm.txt
cat $1_kallisto_igg_rpm.txt | tr '\t' ',' > $1_kallisto_igg_rpm.csv

# save files to proj: kallisto tables
cp $1_kallisto_igg_rpm.txt /work/users/s/h/shuang9/rip/rip_all_stats/combkallisto/
cp $1_kallisto_igg_rpm.csv /work/users/s/h/shuang9/rip/rip_all_stats/combkallisto/