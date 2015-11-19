import pandas as pd
import numpy as np



def ready():
    filename = 'fitnorm.gLoess2Dt25s04.txt'
    data = pd.read_csv(filename,sep='\t')

    data = data.iloc[:,2:]
    data = data.iloc[1:15,1:20]

    gene_num = len(data)
    time_num = len(data.transpose())

    A = []

    for t in range(gene_num):
        A.append(pd.DataFrame(np.random.rand(gene_num,gene_num)))

    return (gene_num,time_num,data,A)
