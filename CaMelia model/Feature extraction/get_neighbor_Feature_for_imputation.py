# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 09:16:08 2020

@author: tjx
"""

from __future__ import division
from sys import argv
import math
import pandas as pd
import numpy as np
import os
import time

from numba import jit
     
from concurrent.futures import ProcessPoolExecutor,as_completed
import warnings
warnings.filterwarnings('ignore')
##############################################################
def reduce_mem(df):
    #starttime = time.time()
    numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
    #start_mem = df.memory_usage().sum() / 1024**2
    for col in df.columns:
        col_type = df[col].dtypes
        if col_type in numerics:
            c_min = df[col].min()
            c_max = df[col].max()
            if pd.isnull(c_min) or pd.isnull(c_max):
                continue
            if str(col_type)[:3] == 'int':
                if c_min > np.iinfo(np.int8).min and c_max < np.iinfo(np.int8).max:
                    df[col] = df[col].astype(np.int8)
                elif c_min > np.iinfo(np.int16).min and c_max < np.iinfo(np.int16).max:
                    df[col] = df[col].astype(np.int16)
                elif c_min > np.iinfo(np.int32).min and c_max < np.iinfo(np.int32).max:
                    df[col] = df[col].astype(np.int32)
                elif c_min > np.iinfo(np.int64).min and c_max < np.iinfo(np.int64).max:
                    df[col] = df[col].astype(np.int64)
            else:
                if c_min > np.finfo(np.float16).min and c_max < np.finfo(np.float16).max:
                    df[col] = df[col].astype(np.float16)
                elif c_min > np.finfo(np.float32).min and c_max < np.finfo(np.float32).max:
                    df[col] = df[col].astype(np.float32)
                else:
                    df[col] = df[col].astype(np.float64)
    #end_mem = df.memory_usage().sum() / 1024**2
    #print('-- Mem. usage decreased to {:5.2f} Mb ({:.1f}% reduction),time spend:{:2.2f} min'.format(end_mem,                                                                                                          100*(start_mem-end_mem)/start_mem,                                                                                                           (time.time()-starttime)/60))
    return df

@jit
def chu(a,n):
    return round(a/(2*n),4)
 

@jit
def corre(a,b):
    return round(sum(a==b)/len(a),4)   

@jit
def getlog2(a):    
    return round(np.log2(a+1.01),4)

@jit
def hah(a,b,c,d):
    return round(np.log2(a+1.01),4)* (1-(abs(b - c)/d)) 

##############################################################
#union
def unionfile(file_dir,meragefiledir,filenames,filenames1,i1):
    if i1+2 >= len(filenames1):
        i2 = len(filenames1)
    else:
        i2 = i1+2   
    for j in range(i1,i2):
        df = pd.DataFrame(np.random.randn(0, 2), columns=['chrom','location'])       
        for i in range(len(filenames)):
            path = r'%s/%s/%s' %(meragefiledir,filenames[i],filenames1[j])
            data = pd.read_csv(path,header=0,sep='\t')

            if i==0:
                df = pd.merge(df,data,how='outer')
            else:
                df = pd.concat([df,data])
            #print (j,i,len(df))
        df = df.drop_duplicates(['chrom','location'])
        df[list(df)[2:]] = df[list(df)[2:]].astype('float16')
        df[list(df)[2:]] = df[list(df)[2:]].round(4)       
        df.to_csv(r'%s/%s' % (file_dir,filenames1[j]),sep='\t',header=True,index=False)     
    return ('jincheng%d:Done!' % i1)

def test(data1,value,neighbor_region,file_dir,colname,i1):
    if i1+2 >= ( len(list(data1))-1 ):
        i2 = len(list(data1))-1 
    else:
        i2 = i1+2
    for i in range(i1,i2): 
        print ("%s: start!" % list(data1)[i] )
        data3 = data1[['chrom','location','%s' % list(data1)[i]]]               
        d_c = pd.DataFrame(np.random.randn(0, 2+2*neighbor_region), columns=['chrom','location']+colname)  
        
        for j in range(len(value)): 
            data2 = data3[data3['chrom'] ==  value[j]]  
            data2 = data2.sort_values(by='location') 
            
            serchtable = list(data2[list(data2)[2]])
            serchtable_l = list(data2[list(data2)[1]])
            
            data2 = data2.dropna(subset=['%s' % list(data3)[2]]) 
            
            data2_value =data2.values   

            #
            if len(data2) > 2*neighbor_region:
                loca = list(data2['location'])
                data2_value = data2.values
                
                data_r = pd.DataFrame(data2_value[neighbor_region:len(loca)-neighbor_region,0:2])
                data_r = data_r.rename(columns={0:'chrom',
                                1:'location'})
                
                loca_for_imputation = []        

                
                data_neigh = pd.DataFrame(np.random.randn(0, 2+2*neighbor_region), columns=['chrom','location']+colname)
                data_neigh['location'] = serchtable_l
                data_neigh['chrom'] = value[j]
                data_neigh = data_neigh.values
                
                s = -1
                ss = []
                for l in range(len(serchtable)):                    
                    if not(np.isnan(serchtable[l])):
                        s += 1
                    else:
                        s = s
                    ss.append(s)
                                        
                    if len(ss)>= 2 :
                        if ss[-1]>=len(data2_value)-neighbor_region:
                            break
                        if (ss[-1] == ss[-2]) and (ss[-1] != -1):
                            if (ss[-1] >= neighbor_region) and (ss[-1]<len(data2_value)-neighbor_region):                                
                                location = serchtable_l[l]
                                location_1 = list( data2_value[ (ss[-1]-9):(ss[-1]+1),1] )
                                location_2 = list( data2_value[ (ss[-1]+1):(ss[-1]+11),1] )
                                if ( location_2[-1]-location )>= ( location- location_1[0] ):
                                    maxdis = location_2[-1]-location 
                                else:
                                    maxdis = location- location_1[0] 
                                            
                                loca = location_1+location_2
                                ind = list( range((ss[-1]-9),(ss[-1]+1)) ) + list( range((ss[-1]+1),(ss[-1]+11)) )
                                mmm = list(data2_value[ind,2])
                                          
                                for p in range(len(mmm)):
                                   m4 = hah(mmm[p],loca[p],location,maxdis)                      
                                   data_neigh[l,p+2] = round(m4,4)
                                         
                                loca_for_imputation.append(location)
                
                data_neigh = pd.DataFrame(data_neigh,columns=['chrom','location']+colname)
                data_neigh = data_neigh.dropna(axis=0,how='any')
                                                                
            if len(data_neigh) != 0:               
                d_c = pd.concat([d_c,data_neigh])  
          
        d_c.to_csv(r'%s/%s_neighbor_methFeature.txt' %  (file_dir,list(data1)[i]) ,sep='\t',header=True,index=False)   
        print ("%s: Done!" % list(data1)[i] )
    return ('jincheng%d:Done!' % i1)


    
##############################################################
#get neighbor_methFeature
if __name__ == '__main__':

    #All data located in the same directory
    DataPath = r'%s' % argv[1]
    #Input data     
    InputDataName = '%s' % argv[2]  
    #neighboring range
    neighbor_region = int(argv[3])

    gse = DataPath
    ff  = InputDataName	
    
    start = time.clock()   

    
    file_dir = r'%s/Forimputation/neighbor_methFeature_%d' % (gse,neighbor_region)
    if not os.path.exists(file_dir):
        os.makedirs(file_dir) 

    path = r'%s/%s' % (gse,ff)
    data = pd.read_csv(path,header=0,sep='\t')
    data = reduce_mem(data)
    if list(data)[0] != 'chrom':
        del data['%s' % list(data)[0]]
    
    #human or mouse
    if 'chr22' in list(data['chrom']):
        #human chrome
        ch = []
        for i in range(1,23):
            ch.append('chr%d' % i)
        ch.append('chrX')
        ch.append('chrY')
        
        #1
        ch1 = []
        for i in range(1,7):
            ch1.append('chr%d' % i)
        #2
        ch2 = []
        for i in range(7,14):
            ch2.append('chr%d' % i)  
        #3
        ch3 = []
        for i in range(14,21):
            ch3.append('chr%d' % i)    
        #4
        ch4 = []
        for i in range(21,23):
            ch4.append('chr%d' % i) 
        ch4.append('chrX')
        ch4.append('chrY')
        
        chh = []
        chh.append(ch1)
        chh.append(ch2)
        chh.append(ch3)
        chh.append(ch4)
    else:
        #mouse chrome
        ch=[]
        for i in range(1,20):
            ch.append('chr%d' % i)
        ch.append('chrX')     
        ch.append('chrY')	
    		
        #1
        ch1 = []
        for i in range(1,7):
            ch1.append('chr%d' % i)
        #2
        ch2 = []
        for i in range(7,14):
            ch2.append('chr%d' % i)    
        #3
        ch3 = []
        for i in range(14,20):
            ch3.append('chr%d' % i) 
        ch3.append('chrX')     
        ch3.append('chrY')     
    
        chh = []
        chh.append(ch1)
        chh.append(ch2)
        chh.append(ch3)
                        

    data = data[data.chrom.isin(ch)]         
    data = data.drop_duplicates(['chrom','location'])
    #colnumes
    colname = []
    for i in range(2*neighbor_region):
        colname.append('neighbor_%d' % i)    
    #####################

    a = len(list(data))-3
    ll = []
    k = 0
    for i in range(int(math.floor(a/2))):
        k += 2
        ll.append(k)
    if a%2 != 0:
        ll.append(k+1) 
    for k in range(len(chh)):
       print ('Block%d : start!' % k)
       file_dir_Block = r'%s/Block%d' %  (file_dir,k)
       if not os.path.exists(file_dir_Block):
           os.makedirs(file_dir_Block)    
       data1 = data[data.chrom.isin(chh[k])] 
       value = chh[k]        
       with ProcessPoolExecutor(max_workers=len(chh[k])) as pool:
            futures = [pool.submit(test,data1,value,neighbor_region,file_dir_Block,colname,i) for i in ll]                               
            for j in as_completed(futures):
                print(j.result())    

    elapsed = (time.clock() - start)
    print("Time used: %d s" % round(elapsed,4))
    
    del data
    del data1

    ####################################    
    #
    print("union neighbor methFeature: Start!" )

    meragefiledir = r'%s' % file_dir
    filenames=os.listdir(meragefiledir)
    
    meragefiledir1 = r'%s/%s' % (file_dir,filenames[0])
    filenames1=os.listdir(meragefiledir1)

    file_dir = r'%s/Forimputation/neighbor_methFeature_%d/localRegion_%d' % (gse,neighbor_region,neighbor_region)
    if not os.path.exists(file_dir):
        os.makedirs(file_dir)  
  
    a = len(filenames1)
    ll = [0]
    k = 0
    for i in range(int(math.floor(a/2))):
        k += 2
        ll.append(k)
    if a%2 != 0:
        ll.append(k+1)
    
    with ProcessPoolExecutor() as pool:
        futures = [pool.submit(unionfile,file_dir,meragefiledir,filenames,filenames1,i) for i in ll]
   
         
    print("union neighbor methFeature: Finished!" )     


 

