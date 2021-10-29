#### Extract hits of putative receptors from blastn result sample by sample; and then output in one file;
##### Author yyzou ######

import pandas as pd
import re
import os
from sys import argv

def process(sample,receptor,path):
    df = pd.DataFrame(columns=receptor[0],index=sample[0])
    for m in range(sample.shape[0]):
        i = sample[0].loc[m]
        f = path + i + "_merge.blast"
        if os.path.getsize(f):
            r = pd.read_csv( f ,sep="\t",header=None)
            r[12]=r.apply(lambda x:x[1].split("_")[0] ,axis=1)
       ## 1C = 978*10^6 bp 

       #### normalization ######
            genome_length = sample[2].loc[m]*978*1000000
            
            sample_length = float(sample[3].loc[m]) ## total length of one sample

        
            for n in range(receptor.shape[0]):
                j = receptor[0].loc[n]  ## receptor 
            
                reads_num = float(len(r[r[12]==j][0].unique()))
                gene_length = float(receptor[1].loc[n])
            
                df[j].loc[i] = (reads_num/sample_length)*(genome_length/gene_length)
       #         df[j].loc[i] = reads_num
    return df



receptor = "/home/yyzou/COVID19/genome/seq/receptor.length"
rec = pd.read_csv(receptor,header=None,sep="\t")
sample = pd.read_csv(argv[1],header=None,sep="\t")  ## sample : dataframe ，col=["sra_id","species_name","genome_length",sequencing_length"]

rec_result = process(sample,rec,argv[2])   ### argv[2] : path where store the blast result
rec_result.to_csv(argv[3],sep="\t")       ### argv[3] : output 
