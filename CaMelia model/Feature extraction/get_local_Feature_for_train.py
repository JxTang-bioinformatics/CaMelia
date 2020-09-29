# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 09:45:08 2020

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

####################################

@jit
def chu(a,n):
    return round(a/(2*n),4)
 

@jit
def corre(a,b):
    return round(sum(a==b)/len(a),4)   

@jit
def getlog2(a):    
    return round(np.log2(a+1.01),4)


####################################
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
#get local_methFeature
def test(data,value,neighbor_region,file_dir,file_dir_for_impu,i1):    
    if i1+2 >= ( len(list(data))-1 ):
        i2 = len(list(data))-1 
    else:
        i2 = i1+2
    for i in range(i1,i2):
        print ('%s: start!' % list(data)[i])       
        data_all = pd.DataFrame(np.random.randn(0, 2), columns=['chrom','location'])  
        data_all_r_for_impu = pd.DataFrame(np.random.randn(0, 2), columns=['chrom','location'])      
        for k in range(i+1,len(list(data))):
            data1 = data[['%s' % list(data)[0],'%s' % list(data)[1],'%s' % list(data)[i],'%s' % list(data)[k]]]
            data_chr = pd.DataFrame(np.random.randn(0, 4), columns=['chrom','location','%s_%s_r' % (list(data1)[2],list(data1)[3]),'%s_%s_methy' % (list(data1)[2],list(data1)[3])])
            data_r_for_impu = pd.DataFrame(np.random.randn(0, 3), columns=['chrom','location','%s_%s_r' % (list(data1)[2],list(data1)[3])])
            for j in range(len(value)): 
                data2 = data1[data1['chrom'] ==  value[j]] 
                data2 = data2.dropna(subset=list(data2)[2:])
                data2 = data2.sort_values(by='location')   

                if len(data2) > 2*neighbor_region:
                    loca = list(data2['location'])
                    data2_value = data2.values
                    
                    data_r = pd.DataFrame(data2_value[neighbor_region:len(loca)-neighbor_region,0:2])
                    data_r = data_r.rename(columns={0:'chrom',1:'location'})
                    
                    data_rr = pd.DataFrame(data2_value[neighbor_region:len(loca)-neighbor_region,0:2])
                    data_rr = data_rr.rename(columns={0:'chrom',1:'location'})
                    
                    r = []
                    rrr = []
                    meth = []
                                           
                    ind_l = list(range(0,neighbor_region))
                    ind_r = list(range(neighbor_region+1,2*neighbor_region+1))                    
                    ind = ind_l + ind_r                    
                    cell_a = pd.Series(data2_value[ind,2])
                    cell_b = pd.Series(data2_value[ind,3])                                 
                    rr = sum(cell_a==cell_b)
                    
                    r.append( chu(rr,neighbor_region)  )
                    meth.append( getlog2(data2_value[neighbor_region,3]) )
                   
                    rrr.append( r[0] )
                                        
                    cell_a_l = list(data2_value[ind_l,2])
                    cell_a_r = list(data2_value[ind_r,2]) 
                                        
                    cell_b_l = list(data2_value[ind_l,3])
                    cell_b_r = list(data2_value[ind_r,3] )     
                    
                    for l in range(neighbor_region+1,len(loca)-neighbor_region):
                        if cell_a_l[0] == cell_b_l[0]:
                            rr = rr -1
                        
                        if data2_value[l-1,2] == data2_value[l-1,3]:
                            rr = rr+1    
                            
                        if cell_a_r[0] == cell_b_r[0]:
                            rr = rr-1   
                            
                        if data2_value[l+neighbor_region,2] == data2_value[l+neighbor_region,3]:
                            rr = rr+1
                       
                        ind_l.pop(0)
                        ind_l.append(l-1)
                        
                        cell_a_l.pop(0)
                        cell_a_l.append( data2_value[l-1,2] )
                        cell_a_r.pop(0)
                        cell_a_r.append( data2_value[l+neighbor_region,2] ) 
                        
                        
                        ind_r.pop(0)
                        ind_r.append( l+neighbor_region )
                        
                        cell_b_l.pop(0)
                        cell_b_l.append( data2_value[l-1,3] )
                        cell_b_r.pop(0)
                        cell_b_r.append( data2_value[l+neighbor_region,3] )  
                        
                        r_for2 = chu(rr,neighbor_region)
                        r.append( r_for2  )
                        rrr.append( r_for2  )
                        meth.append( getlog2(data2_value[l,3]) )
                                                                                   
                    data_r['%s_%s_r' % (list(data2)[2],list(data2)[3]) ]  = r
                    data_r['%s_%s_methy' % (list(data2)[2],list(data2)[3]) ]  = meth      
                    
                    data_rr['%s_%s_r' % (list(data2)[2],list(data2)[3]) ]  = rrr
                    
                    data_chr = pd.merge(data_chr,data_r,how='outer')  
                    data_r_for_impu = pd.merge(data_r_for_impu,data_rr,how='outer')  
            if len(data_chr) != 0:
                data_all = pd.merge(data_all,data_chr,how='outer',on=['chrom','location'])
            if len(data_r_for_impu) != 0:
                data_all_r_for_impu = pd.merge(data_all_r_for_impu,data_r_for_impu,how='outer',on=['chrom','location'])
                                
            print ('%s-%s: Done!' % (list(data1)[2],list(data1)[3]))
            
        data_all.to_csv(r'%s/%s_local_methFeature.txt' %  (file_dir,list(data)[i]) ,sep='\t',header=True,index=False) 
                
        data_all_r_for_impu.to_csv(r'%s/%s_local_r.txt' %  (file_dir_for_impu,list(data)[i]) ,sep='\t',header=True,index=False) 
        print ('%s: Done!' % list(data1)[2])
    return ('jincheng%d:Done!' % i1)

if __name__ == '__main__': 
	
    #All data located in the same directory
    DataPath = r'%s' % argv[1]
    #Input data     
    InputDataName = '%s' % argv[2]  
    #neighboring range
    neighbor_region = int(argv[3])
    #correlation threshold
    threshold = float(argv[4])
   
   
    gse = DataPath
    ff  = InputDataName	
	   
    start = time.clock()
    #gse = r'E:\singlecellimputation\test_HCC_0622'
    #ff = r'GSE65364_Ca.txt'   
    #neighbor_region = 10
    #threshold = 0.8
    
    file_dir = r'%s/local_methFeature_spares' %  gse
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
    
              
    cell_num = []
    for i in range(2,len(list(data))-1):
        cell_num.append(list(data)[i])
    cell_num=list(set(cell_num))  
                            
    data = data.drop_duplicates(['chrom','location'])
    data = data[data.chrom.isin(ch)] 
              
    print('get local methFeature: Start!')

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
        file_dir_for_impu = r'%s/for_impu_r/Block%d' %  (gse,k)
        if not os.path.exists(file_dir_for_impu):
            os.makedirs(file_dir_for_impu)                    
        #split by chrom
        data1 = data[data.chrom.isin(chh[k])] 
        value = chh[k]       
        with ProcessPoolExecutor() as pool:
            #run by bock
            futures = [pool.submit(test,data1,value,neighbor_region,file_dir_Block,file_dir_for_impu,i) for i in ll]
            for j in as_completed(futures):
                print(j.result()) 
    
    elapsed = (time.clock() - start)
    print("Time used: %d s" % round(elapsed,4))

    print("get local methFeature: Finished!" )

    del data
    del data1

    ############################################################################
    #
    print("union local methFeature: Start!" )

    meragefiledir = r'%s' % file_dir
    filenames=os.listdir(meragefiledir)
    
    meragefiledir1 = r'%s/%s' % (file_dir,filenames[0])
    filenames1=os.listdir(meragefiledir1)

    file_dir = r'%s/local_methFeature/localRegion_%d' % (gse,neighbor_region)
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
      
   
    print("union local methFeature: Finished!" )
    ############################################################################
    
    region = neighbor_region
    meragefiledir = r'%s/local_methFeature/localRegion_%d' % (gse,region)
    filenames=os.listdir(meragefiledir)

    path = r'%s/local_methFeature/localRegion_%d/%s' % (gse,region,filenames[0])
    data = pd.read_csv(path,header=0,sep='\t')

    for i in range(1,len(filenames)):
        path = r'%s/%s' % (meragefiledir,filenames[i])
        df = pd.read_csv(path,header=0,sep='\t')
        data = pd.merge(data,df,how='outer',on=['chrom','location'])
    print("union" )        
    ############################################################################    
    
    file_dir_1 = r'%s/local_methFeature_cellbycell/region%d/corr' % (gse,region)
    if not os.path.exists(file_dir_1):
        os.makedirs(file_dir_1)

    file_dir_2 = r'%s/local_methFeature_cellbycell/region%d/methy' % (gse,region)
    if not os.path.exists(file_dir_2):
        os.makedirs(file_dir_2)
		    
    for i in range(len(cell_num)): 
        cell_r = []
        cell_meth = []
        for j in range(2,len(list(data))):
            if (list(data)[j].split('_')[0] == cell_num[i]) or (list(data)[j].split('_')[1] == cell_num[i]):
                if list(data)[j].split('_')[-1] == 'r':
                    cell_r.append(list(data)[j])
                if list(data)[j].split('_')[-1] == 'methy':
                    cell_meth.append(list(data)[j])            
        
        df_r =data[ ['chrom','location'] + cell_r ] 
        df_m = data[ ['chrom','location'] + cell_meth ] 
        
        df_r = df_r.dropna(axis=0,subset=list(df_r)[2:],how='all') 
        df_m = df_m.dropna(axis=0,subset=list(df_m)[2:],how='all') 
        
        df_r.to_csv(r'%s/%s_r.txt' % (file_dir_1,cell_num[i])  ,sep='\t',header=True,index=False)    
        df_m.to_csv(r'%s/%s_m.txt' % (file_dir_2,cell_num[i])  ,sep='\t',header=True,index=False)    
    print("split" )     
    ############################################################################
    for i in range(len(cell_num)):
        #corr
        path = r'%s/%s_r.txt' % (file_dir_1,cell_num[i])
        data_r = pd.read_csv(path,header=0,sep='\t')
        #methy
        path = r'%s/%s_m.txt' % (file_dir_2,cell_num[i])
        data_m = pd.read_csv(path,header=0,sep='\t')
        
        df_r = data_r[['chrom','location']]
    
        name = list(data_r)[2:]
               
        local_matched_r = []
        local_matched_methy = []
        
        data_r.fillna(0, inplace=True)
        data_m.fillna(0, inplace=True)
    
        r_value = data_r.values
        m_value = data_m.values    
        
        col_ind = []
        col_ind_m = []
        for j in range(len(r_value)):       
            a = np.array(r_value[j,2:])
            b = np.array(m_value[j,2:])
            #c = np.where((a>=threshold)&(a==max(a)) )
            c = np.where( a>=threshold )
            
            if len(c[0]) != 0:
                ss = 0
                for k in range(len(c[0])):           
                    ss += a[ c[0][k] ]*b[ c[0][k] ]
                col_ind_m.append(round(ss/len(c[0]),4))    
                if len(c[0]) == 1 :
                    if name[c[0][0]].split('_')[-2] ==  cell_num[i]:
                        col_ind.append(name[c[0][0]].split('_')[0])
                    else:
                        col_ind.append(name[c[0][0]].split('_')[-2])
                else:
                    col_ind.append(len(c[0]))
            else:
                col_ind.append(np.nan)
                col_ind_m.append(np.nan)

        df_r['aver_meth'] = col_ind_m
        df_r['matched_cell'] = col_ind
                
        file_dir = r'%s/region%d_localmatched_morethan08' % (gse,region)
        if not os.path.exists(file_dir):
            os.makedirs(file_dir)
        
        path = r'%s/%s.txt' % (file_dir,cell_num[i])
        
        df_r.to_csv(path, sep='\t', header=True, index=False)   
    print("split: Done!" )