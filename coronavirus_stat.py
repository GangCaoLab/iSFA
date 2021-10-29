#### Extract hits of 38 coronaviruses from blastn result sample by sample; and then output in one file;

##### Author yyzou ######

import pandas as pd
import re
import os
from sys import argv

def process(sample,virus_id,path):
    df = pd.DataFrame(columns=virus_id[1],index=sample[0])
    for i in sample[0]:
        f = path + i + "_merge.blast"
        
        ## check empty file
        if os.path.getsize(f):
            r = pd.read_csv( f ,sep="\t",header=None)
            for j in range(virus_id.shape[0]):
                df[virus_id[1].loc[j]].loc[i] = len(r[r[1]==virus_id[0].loc[j]][0].unique())
    return df

#virus_id = pd.read_csv("/home/yyzou/COVID19/genome/id_virus",header=None,sep="\t")

virus_id = pd.read_csv("/home/yyzou/COVID19/genome/added_coronavirus/add_id_virus",header=None,sep="\t")
print(virus_id)
sample = pd.read_csv(argv[1],header=None)  ## argv[1]: sample, list of sra id

rec_result = process(sample,virus_id,argv[2])   ### argv[2] : path where store the blast result
rec_result.to_csv(argv[3],sep="\t")   ### output file
