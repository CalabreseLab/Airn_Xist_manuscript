#!/bin/bash

#SBATCH -J peak_%j
#SBATCH -p general
#SBATCH --mem=50G
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=12:00:00
#SBATCH -o peak_%j.out
#SBATCH -e peak_%j.err

### Goals: find true peaks of rip enriched over 2*igg;
    # count reads within each peak, including reads from all rips, igg, and individual rips (if rip is more than one);
    # classify peaks based on exonic/intronic/utr5/utr3, and coding/lncRNA/other

## input files needed: 
    # rips in $2
        # If you supply >1 rip files, they will be indexed in alphabetical order (the order of "ls" command).
        # Ex. For alyref.fastq and ciz1.fastq, alyref will be rip1 and ciz1 will be rip2.
    # igg alignment:
        # /proj/calabrlb/users/Zhiyue/22_sp/rip_sam/igg_Aligned_filteredsq30.out.sam
    # randomization code:
        # /proj/calabrlb/users/Zhiyue/22_sp/code/macs_strand_rand_sam.pl
    # gencode annotation (with utf info):
        # /proj/calabrlb/users/Zhiyue/22_sp/gtf/gencode.vM25.basic.annotation.utr.gtf

## Example command:
# sbatch 230301_3pr_classify.sh ezh2 "/proj/calabrlb/users/Zhiyue/23_sp/0221_3pr_rip/*ezh2*"
    # $1 = rbp name (all lower cases) e.g., alyref
    # $2 = all rip files; you can put them in someDirectory and use "someDirectory/*"

## load modules
module load star/2.5.4b
module load samtools
module load macs/2.2.7.1
module load subread/2.0.6
module load bedtools/2.30
#module load python/2.7.12

#### Part 1: from rip, STAR genome alignment and MACS peak-calling 

## do following work in directory named $1_classify
mkdir $1_classify
cd $1_classify

### 1A. STAR alignment

## STAR alignment of all rips
# move rip files here
mkdir $1
cp ${2} ${1}
gunzip $1/*.fastq.gz

# use STAR aligner to align reads to mm10 mouse genome.
rip_filenames=$(ls -m $1/*.fastq | sed -r 's/\s+//g' | tr -d '\n')
STAR --genomeDir /proj/seq/data/STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix $1_ --readFilesIn $rip_filenames

# Filter and output to sam file
samtools view -h -Sq 30 $1_Aligned.out.sam > $1_Aligned_filteredsq30.out.sam

## STAR alignment of individual rips
# decide if there are more than 1 rip
rip_file_num=$(ls $1 | wc -l)
if [ $rip_file_num -gt 1 ]
then

    IFS=',' read -r -a rip_arr <<< "$rip_filenames" # make array with rip file names

    # align each rip file
    for i in ${!rip_arr[@]}; do # loop on indices of array
        rip_name="${rip_arr[$i]}" # for each file name
        echo The rip${i} file is $rip_name. >> $1_rip_info.txt # NOTE: check in $1_rip_info.txt to see which rip gets assigned to which number!
        j=$(($i + 1)) # index+1 so that it's 1 based for the file name

        # align to mm10 genome
        STAR --genomeDir /proj/seq/data/STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix $1_rip${j}_ --readFilesIn $rip_name
        
        # filter to sam file
        samtools view -h -Sq 30 $1_rip${j}_Aligned.out.sam > $1_rip${j}_Aligned_filteredsq30.out.sam
    done

fi


### 1B. MACS peak-calling

# split sam by strand
samtools view -h -F 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_neg.out.sam
samtools view -h -f 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_pos.out.sam

# load randomization code
cp /proj/calabrlb/users/Zhiyue/22_sp/code/macs_strand_rand_sam.pl .

# randomize strand
perl macs_strand_rand_sam.pl $1_Aligned_pos.out.sam $1_Aligned_pos.out.rand.sam
perl macs_strand_rand_sam.pl $1_Aligned_neg.out.sam $1_Aligned_neg.out.rand.sam

# convert sam to bam to avoid error in macs callpeak
samtools view -S -b $1_Aligned_pos.out.rand.sam > $1_Aligned_pos.out.rand.bam
samtools view -S -b $1_Aligned_neg.out.rand.sam > $1_Aligned_neg.out.rand.bam

# macs using bam file
    # merge size = 304 = 4 * 76
macs2 callpeak -t $1_Aligned_pos.out.rand.bam --keep-dup all --broad --broad-cutoff 0.3 --max-gap 76 --outdir $1_pos_peaks -n $1_pos
macs2 callpeak -t $1_Aligned_neg.out.rand.bam --keep-dup all --broad --broad-cutoff 0.3 --max-gap 76 --outdir $1_neg_peaks -n $1_neg 
    # add trackline for UCSC
cat $1_neg_peaks/$1_neg_peaks.broadPeak | sed "1 i\track type=broadPeak name="$1_peaks_neg" description="$1_peaks_neg" nextItemButton=on" > $1_neg_peaks_ucsc.broadPeak
cat $1_pos_peaks/$1_pos_peaks.broadPeak | sed "1 i\track type=broadPeak name="$1_peaks_pos" description="$1_peaks_pos" nextItemButton=on" > $1_pos_peaks_ucsc.broadPeak



#### Part 2: featureCounts count number of reads in each peak

## convert broadPeak to SAF file in preparation for featureCounts
cat $1_neg_peaks_ucsc.broadPeak | sed 1d | awk -v OFS="\t" '{print $1,$2,$3,"-"}' > $1_saf.txt
cat $1_pos_peaks_ucsc.broadPeak | sed 1d | awk -v OFS="\t" '{print $1,$2,$3,"+"}' >> $1_saf.txt
cat $1_saf.txt | awk -v OFS="\t" '{print NR,$0}' | sed '1 i\GeneID\tChr\tStart\tEnd\tStrand' > $1_peaks.saf

## featureCounts with SAF
# count reads in each peak for all_rips and igg
featureCounts -s 2 -F SAF -a $1_peaks.saf -o $1_fc $1_Aligned_filteredsq30.out.sam
featureCounts -s 2 -F SAF -a $1_peaks.saf -o $1_igg_fc /proj/calabrlb/users/Zhiyue/22_sp/rip_sam/igg_Aligned_filteredsq30.out.sam

# extract raw reads
cut -f7 $1_igg_fc > $1_igg_counts.txt

# count reads in each peak for individual rips
if [ $rip_file_num -gt 1 ]
then

    for i in ${!rip_arr[@]}; do # loop on indices
        j=$(($i + 1)) # index+1 so that it's 1 based
        featureCounts -s 2 -F SAF -a $1_peaks.saf -o $1_rip${j}_fc $1_rip${j}_Aligned_filteredsq30.out.sam
        
        # extract raw reads
        cut -f7 $1_rip${j}_fc > $1_rip${j}_counts.txt
        
        # prep for concatenate count files
        if [ $j == 1 ]
        then
        rip_count_filename="$1_rip1_counts.txt"
        else
        rip_count_filename="${rip_count_filename} $1_rip${j}_counts.txt"
        fi
    done

    # make summary file: all_rip, igg, individual rips
    paste $1_fc $1_igg_counts.txt $rip_count_filename | sed 1d > $1_fc.txt

else

    # make summary file: rip, igg (when only one rip file)
    paste $1_fc $1_igg_counts.txt | sed 1d > $1_fc.txt

fi



#### Part 3: find true peaks enriched over ctrl (rpm > 2 * ctrl_rpm)

## calculate rpm = (count * 1000,000 / reads)

# rpm conversion for all rips:
# var: reads = number of lines in fastq files / 4
fastq_lin_num=$(wc -l $1/*$1*.fastq | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
(( reads = fastq_lin_num / 4 ))
cat $1_fc.txt | sed 1d | awk -v OFS="\t" -v reads=$reads '{print $0, (($7*1000000)/reads)}' > $1_fc_allRpm.txt

# rpm conversion for igg:
    # reads_igg 
        # only calculate once:
            # cp /proj/calabrlb/users/Zhiyue/22_sp/rip/*igg* .
            # gunzip *igg*.fastq.gz
            # fastq_lin_num_igg=$(wc -l *igg*.fastq | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
            # (( reads_igg = fastq_lin_num_igg / 4 ))
            # echo $reads_igg
cat $1_igg_counts.txt  | sed 1,2d | awk -v reads_igg=124344953 '{print (($1*1000000)/reads_igg)}' > $1_igg_rpm.txt

# rpm conversion for individual rips:

if [ $rip_file_num -gt 1 ]
then

    for i in ${!rip_arr[@]}; do # loop on indices
        j=$(($i + 1)) # index+1 so that it's 1 based
        # count reads
        fastq_lin_num=$(wc -l ${rip_arr[$i]} | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
        (( reads = fastq_lin_num / 4 ))
        # rpm conversion
        cat $1_rip${j}_counts.txt | sed 1,2d | awk -v reads=$reads '{print (($1*1000000)/reads)}' > $1_rip${j}_rpm.txt

        # prep for concatenate rpm files and file header
        if [ $j == 1 ]
        then
        rip_count_header="$1_rip1_counts"
        rip_rpm_header="$1_rip1_rpm"
        rip_rpm_filename="$1_rip1_rpm.txt"
        else
        rip_count_header="${rip_count_header}\t$1_rip${j}_counts"
        rip_rpm_header="${rip_rpm_header}\t$1_rip${j}_rpm"
        rip_rpm_filename="${rip_rpm_filename} $1_rip${j}_rpm.txt"
        fi
    done

    # make summary of rpm: all rips, igg, individual rips
    paste $1_fc_allRpm.txt $1_igg_rpm.txt $rip_rpm_filename | sed "1 i\GeneID\tChr\tStart\tEnd\tStrand\tLength\t$1_counts\tigg_counts\t${rip_count_header}\t$1_rpm\tigg_rpm\t${rip_rpm_header}" > $1_fc_rpm.txt

else

    # make summary of rpm: rip, igg
    paste $1_fc_allRpm.txt $1_igg_rpm.txt | sed "1 i\GeneID\tChr\tStart\tEnd\tStrand\tLength\t$1_counts\tigg_counts\t$1_rpm\tigg_rpm" > $1_fc_rpm.txt

fi

## filter peaks with rpm 2 folds over igg

if [ $rip_file_num -gt 1 ]
then
    ((rip_all_col = rip_file_num + 9)) # calculate column number of rip_rpm
    ((igg_col = rip_all_col + 1)) # calculate column number of igg_rpm
    cat $1_fc_rpm.txt | awk -v c1=$rip_all_col -v c2=$igg_col  '{if(NR == 1) print $0; else if($c1 > 2 * $c2) print $0}' > $1_fc_rpm_2igg.txt
else
    cat $1_fc_rpm.txt | awk '{if(NR == 1) print $0; else if($9 > 2 * $10) print $0}' > $1_fc_rpm_2igg.txt
fi

# Note: all following commands are about peaks and do not take into account how many rip files you have.

# convert enricheded peaks to BED format
cat $1_fc_rpm_2igg.txt | sed 1d | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > $1_peaks_2igg.bed

# add trackline for UCSC
cat $1_peaks_2igg.bed | sed "1 i\track name=$1_peaks_2igg description=$1_peaks_2igg colorByStrand='<0,0,255 255,0,0>'" > $1_peaks_2igg_ucsc.bed


