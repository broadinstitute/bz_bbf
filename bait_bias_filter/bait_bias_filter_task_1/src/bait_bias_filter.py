# -*- coding: utf-8 -*-
import sys
import pandas as pd
from scipy import stats

"""
Created on Thu Feb 21 12:23:37 2019

@author: bzhitomi

Resolve bait bias by removing mutations from the biased type from 
the lowest alternative allele counts until a desired P-value is reached for a
one sided binomial test.
clear_bait_bias does the above for a single sampla MAF
multi_sample_clear_bait_bias uses clear_bait_bias to itterate over multiple 
samples in a single MAF.

clear_bait_bias(MAF,ref_base,alt_base,min_pval) accepts:
    - MAF: MAF file in the form of a dataframe with data for one sample.
    - ref_base: 'A','T','C' or 'G' - the referance base in the suspect mutation.
    - alt_base: 'A','T','C' or 'G' - the alternative base in the suspect mutation.
    - min_pval: the P-value to reach by removing mutations
clear_bait_bias returns:  
    - MAF: MAF file with the lowest count mutations from the suspect mutation 
    type removed until the point where a one sided binomial test p_value is 
    greater than min_pval.
    - Counter-1: the highest count of mutations removed.
    
multi_sample_clear_bait_bias(MAF,ref_base,alt_base,min_pval) accepts:
    - MAF: MAF file in the form of a dataframe with data for multiple.
    - ref_base: 'A','T','C' or 'G' - the referance base in the suspect mutation.
    - alt_base: 'A','T','C' or 'G' - the alternative base in the suspect mutation.
    - min_pval: the P-value to reach by removing mutations
multi_sample_clear_bait_bias returns:
    - result_MAF: MAF file with the lowest count mutations from the suspect mutation 
    type removed until the point where a one sided binomial test p_value is 
    greater than min_pval.
    - removed: a dictionary with the highest count of mutations removed for each sample.
"""

compdict={'A':'T','T':'A','C':'G','G':'C'}

def clear_bait_bias (MAF,ref_base,alt_base,min_pval):
    mut_counts=MAF[(MAF["Reference_Allele"]==ref_base) & (MAF["Tumor_Seq_Allele2"]==alt_base)].count()["Reference_Allele"]
    rev_mut_count=MAF[(MAF["Reference_Allele"]==compdict[ref_base]) & (MAF["Tumor_Seq_Allele2"]==compdict[alt_base])].count()["Reference_Allele"]
    counter=min(MAF["t_alt_count"])
    p0=stats.binom_test([mut_counts,rev_mut_count],p=0.5,alternative='greater')
    while stats.binom_test([mut_counts,rev_mut_count],p=0.5,alternative='greater')<min_pval:
        MAF=MAF[((MAF["t_alt_count"]>counter) & (MAF["Reference_Allele"]==ref_base) & (MAF["Tumor_Seq_Allele2"]==alt_base)) | (MAF["Reference_Allele"]!=ref_base) | (MAF["Tumor_Seq_Allele2"]!=alt_base)]
        mut_counts=MAF[(MAF["Reference_Allele"]==ref_base) & (MAF["Tumor_Seq_Allele2"]==alt_base)].count()["Reference_Allele"]
        rev_mut_count=MAF[(MAF["Reference_Allele"]==compdict[ref_base]) & (MAF["Tumor_Seq_Allele2"]==compdict[alt_base])].count()["Reference_Allele"]
        counter+=1
    
    p1=stats.binom_test([mut_counts,rev_mut_count],p=0.5,alternative='greater')
    return (MAF,counter-1,p0,p1)


def multi_sample_clear_bait_bias(MAF_file,ref_base,alt_base,name,min_pval):

    MAF=pd.read_csv(MAF_file,sep='\t',index_col=None,low_memory=False)
    samples=set(MAF["Tumor_Sample_Barcode"])
    result_MAF=MAF.iloc[0:0].copy()
    removed={}
    for sample in samples:
        samp_MAF,cnt,p0,p1=clear_bait_bias(MAF[(MAF["Tumor_Sample_Barcode"]==sample)],ref_base,alt_base,min_pval)
        result_MAF=result_MAF.append(samp_MAF)
        removed[sample]=cnt
        
    ncut = len(MAF)-len(result_MAF)
    f1=open(name+".BaitBiasfilt_"+ref_base+"to"+alt_base+"_ncut.txt",'w')
    f1.write(str(ncut))
    f1.close()
    result_MAF.to_csv(path_or_buf=name+".BaitBiasfilt_"+ref_base+"to"+alt_base+".maf",sep='\t', index=None)
    
    return(0)


def main():
    # print command line arguments
    for arg in sys.argv[1:]:
        print(arg)

    MAF_file=sys.argv[1]
    ref_base=sys.argv[2]
    alt_base=sys.argv[3]
    name=sys.argv[4]

    if len(sys.argv)>5:
        min_pval=float(sys.argv[5])
    else:
        min_pval=0.1
    
    
    ret = multi_sample_clear_bait_bias(MAF_file,ref_base,alt_base,name,min_pval)

    return(ret)

if __name__ == "__main__":
    main()
