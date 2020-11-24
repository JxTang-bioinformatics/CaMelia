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
            df = pd.merge(df,data,how='outer')
        df = df.drop_duplicates(['chrom','location'])
        df.to_csv(r'%s/%s' % (file_dir,filenames1[j]),sep='\t',header=True,index=False)     
    return ('jincheng%d:Done!' % i1)


def test(data,value,neighbor_region,file_dir,colname,i1):
    if i1+2 >= ( len(list(data))-1 ):
        i2 = len(list(data))-1 
    else:
        i2 = i1+2
    for i in range(i1,i2): 
        print ("%s: start!" % list(data)[i] )
        data3 = data[['chrom','location','%s' % list(data)[i]]]               
        
        d_c = pd.DataFrame(np.random.randn(0, 2+2*neighbor_region), columns=['chrom','location']+colname) 
        
        for j in range(len(value)): 
            data2 = data3[data3['chrom'] ==  value[j]]  
            data2 = data2.sort_values(by='location')  
            data2 = data2.dropna(subset=['%s' % list(data3)[2]]) 
            data2 = data2.reindex()    
            
            #Èç¹ûÎ»µã²»Âú×ãlocalregionµÄ·¶Î§ÔòÌáÇ°ÍË³ö
            if len(data2) > 2*neighbor_region:
                loca = list(data2['location'])
                data2_value = data2.values
                
                data_r = pd.DataFrame(np.random.randn(0, 2+2*neighbor_region), columns=['chrom','location']+colname)
                data_r['chrom'] = data2_value[neighbor_region:len(loca)-neighbor_region,0]
                data_r['location'] = data2_value[neighbor_region:len(loca)-neighbor_region,1]
                
                data_r_value = data_r.values
                       
                for l in range(neighbor_region,len(loca)-neighbor_region):
                    ind = list(range(l-neighbor_region,l+neighbor_region+1))
                    ind.remove(l)
                    location_1 = list(data2_value[l-neighbor_region:l,1])
                    location = data2_value[l,1]
                    location_2 = list(data2_value[l+1:l+neighbor_region+1,1])
                    
                    if ( location_2[-1]-location )>= ( location- location_1[0] ):
                        maxdis = location_2[-1]-location 
                    else:
                        maxdis = location- location_1[0] 
                                
                    loca = location_1+location_2
                    mmm = list(data2_value[ind,2])
                                
                    for p in range(len(mmm)):
                       m4 = hah(mmm[p],loca[p],location,maxdis)
                       data_r_value[l-neighbor_region,p+2] = round(m4,4)
                                                                                                                    
                data_r = pd.DataFrame(data_r_value,columns=['chrom','location']+colname)
                
                d_c = pd.merge(d_c,data_r,how='outer')  
          
        d_c.to_csv(r'%s/%s_neighbor_methFeature.txt' %  (file_dir,list(data)[i]) ,sep='\t',header=True,index=False)   
        print ("%s: Done!" % list(data)[i] )
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
    #gse = r'E:\singlecellimputation\test_HCC_0622'
    #ff = r'GSE65364_Ca.txt'   
    #neighbor_region = 10
    
    file_dir = r'%s/neighbor_methFeature_%d' % (gse,neighbor_region)
    if not os.path.exists(file_dir):
        os.makedirs(file_dir) 

    path = r'%s/%s' % (gse,ff)
    data = pd.read_csv(path,header=0,sep='\t')
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
       #°´È¾É«Ìå·Ö¸îdata
       data1 = data[data.chrom.isin(chh[k])] 
       value = chh[k]        
       with ProcessPoolExecutor(max_workers=8) as pool:
            #´«µÝdata-·Ö¿é
            futures = [pool.submit(test,data1,value,neighbor_region,file_dir_Block,colname,i) for i in ll]
                           
    
            for j in as_completed(futures):
                print(j.result()) # ¶ÔÓ¦½ø³ÌÍê³ÉË³ÐòÊä³ö£¬Êä³ö2£¬3£¬5£¬6£¬10    

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

    file_dir = r'%s/neighbor_methFeature_10/localRegion_%d' % (gse,neighbor_region)
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
        for j in as_completed(futures):
            print(j.result())   
         
    print("union neighbor methFeature: Finished!" )     
    
    


 

