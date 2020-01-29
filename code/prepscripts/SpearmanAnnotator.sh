#!/bin/bash

echo -e "Target_Gene\tHit_Gene\tSpearman_Correlation\tSpearman_pVal\tFDR_CorPval\tHit_Tissue\tHit_Data_Set" > $tissue"_SpearmanCor.txt"
awk -v tissue=$tissue -vcohort=$cohort '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"tissue"\t"cohort}' Attemp.txt >> $tissue"_SpearmanCor.txt"
awk -v tissue=$tissue -vcohort=$cohort '{print $2"\t"$1"\t"$3"\t"$4"\t"$5"\t"tissue"\t"cohort}' Attemp.txt >> $tissue"_SpearmanCor.txt"

