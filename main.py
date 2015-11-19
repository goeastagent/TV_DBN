import pandas as pd
import numpy as np
import math 

import os
import sys
# Author : Garam Lee , Ajou University.
# Nov 7th. 2015
# 
# This program builds Time-Varing Dynamic Bayesian Network from time-course data
# Algorithm used for construction of network is refered to "Time-Varying Dynamic Bayesian Networks"
# from Le Song, Mladen Kolar and Eric P.Xing.

# sample file name is "fitnorm.gLoess2Dt25s04.txt"
#
#
# Usage :
# console$ python main.py [name of sample file] [parameter lambda] [kernel parameter h]
# paremeter lambda referes to optimizaing sparse model as shown in equation (5) in paper.
# kernel parameter is used for reweighting the observations as shown in paper describing w.
#
# This program creates new file in form of "[name of sample file]_output.csv",
# which is coefficient vectors to build time-varying network at individual time points.  





# This function is the main function to produce coefficients to build time-varying network.
# Returns Dataframe that contains coefficient vectors at individual time points.

def estimate_TV_DBN(data,pLambda,pBandwidth):
    gene_num = len(data)   # the number of genes
    time_num = len(data.transpose()) # the number of time points

    pBandwidth = float(time_num)/7
    A = []                           # coefficient vectors
    # initialize the A vector
    A.append(pd.DataFrame(np.random.rand(gene_num,gene_num)))
    for t in range(time_num):
        A.append(pd.DataFrame(np.zeros(shape=(gene_num,gene_num))))

    for p in range(gene_num):
        for t in range(1,time_num+1):
            A[t].iloc[p] = A[t-1].iloc[p].copy()
            w_t = []

            # Scale time series
            # for i in range(time_num):
            #     w_t.append(get_weight(t,i,pBandwidth,time_num))
        
            while True:
                prev = A[t].iloc[p].copy()

                # compute and update
                prewt = []
                for i in range(time_num):
                    prewt.append(get_weight(t,i,pBandwidth,time_num))
                    
                
                for j in range(gene_num):
                    # b = (2.0/time_num)*(data.iloc[j,:]*data.iloc[j,:]*prewt).sum()
                    S = 0.0
                    b = 0
                    for i in range(1,time_num):
                        small_S = 0.0
                        w_t = get_weight(t,i,pBandwidth,time_num)
                        b = b + data.iloc[j,i-1]*data.iloc[j,i-1]*w_t
                        for k in range(gene_num):
                            if k==j:
                                continue
                            small_S = small_S+(A[t].iloc[p,k]*data.iloc[k,i-1])
                        S = S +(((small_S - data.iloc[p,i])*data.iloc[j,i-1])*w_t)

                    S = S/time_num*2
                    if abs(S) > pLambda:
                        sign = 0
                        if S-pLambda > 0:
                            sign = +1
                        else:
                            sign = -1
                        A[t].iloc[p,j] = (sign*pLambda - S)/b
                    else:
                        A[t].iloc[p,j] = 0

                #print A[t].iloc[p,:]
                if not(False in list(prev == A[t].iloc[p])):
                    break
            
            print A[t]

    return A[1:]                # Return A except A0


def get_weight(t_star,t,pBandwidth,time_num):
    return calculate_kernel_function(t_star,t,pBandwidth)/calculate_kernel_sum(t_star,pBandwidth,time_num)

def calculate_kernel_function(t_star,t,h):
    temp = math.pow(t_star-t,2)    
    return math.exp(temp*(-1)/float(h))

      
def calculate_kernel_sum(t_star,h,time_num):
    s = 0.0
    for t in range(time_num):
        s = s + calculate_kernel_function(t_star,t,h)
    return s


# main function
if __name__ == "__main__":
    filename = sys.argv[1]
    pLambda = sys.argv[2]
    pBandwidth = sys.argv[3]

    filename = 'fitnorm.gLoess2Dt25s04.txt'
    pLambda = 0.003
    pBandwidth = 40

    data = pd.read_csv(filename,sep='\t')
    data = data.iloc[:,2:]
    data = data.iloc[1:20,1:20]

    TV_network = estimate_TV_DBN(data,pLambda,pBandwidth)

    for i,p in enumerate(TV_network):
        p.to_csv('output/'+filename+'_output' + str(i) + '.csv',sep='\t')
    # TV_network.to_csv(filename+'_'+'output.csv',sep='\t')
