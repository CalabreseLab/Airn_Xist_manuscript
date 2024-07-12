#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBATCH --mem=64g

featureCounts -s 2 -F SAF -a xak_chunks_new.saf -o xak_chunks_prc_new.txt bmi1_rip1_Aligned_filteredsq30.out.sam bmi1_rip2_Aligned_filteredsq30.out.sam epop_rip1_Aligned_filteredsq30.out.sam epop_rip2_Aligned_filteredsq30.out.sam ezh2_rip1_Aligned_filteredsq30.out.sam ezh2_rip2_Aligned_filteredsq30.out.sam hnrnpk_rip1_Aligned_filteredsq30.out.sam hnrnpk_rip2_Aligned_filteredsq30.out.sam igg_rip1_Aligned_filteredsq30.out.sam igg_rip2_Aligned_filteredsq30.out.sam igg_rip3_Aligned_filteredsq30.out.sam igg_rip4_Aligned_filteredsq30.out.sam igg_rip5_Aligned_filteredsq30.out.sam igg_rip6_Aligned_filteredsq30.out.sam jarid2_rip1_Aligned_filteredsq30.out.sam jarid2_rip2_Aligned_filteredsq30.out.sam mtf2_rip1_Aligned_filteredsq30.out.sam mtf2_rip2_Aligned_filteredsq30.out.sam ring1b_rip1_Aligned_filteredsq30.out.sam ring1b_rip2_Aligned_filteredsq30.out.sam rybp_rip1_Aligned_filteredsq30.out.sam rybp_rip2_Aligned_filteredsq30.out.sam suz12_rip1_Aligned_filteredsq30.out.sam suz12_rip2_Aligned_filteredsq30.out.sam