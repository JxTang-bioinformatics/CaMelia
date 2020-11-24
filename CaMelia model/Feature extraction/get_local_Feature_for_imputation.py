# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 16:44:39 2020

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
def unionlargefile(i1,gse,region,chh,cell_num): 
    print ("Block%d:start!" % i1 )
    file_dir_Block = r'%s/Forimputation/local_methFeature_cellbycell/Block%d' %  (gse,i1)
    if not os.path.exists(file_dir_Block):
        os.makedirs(file_dir_Block) 
    print('here1')    
    
    meragefiledir = r'%s/Forimputation/local_methFeature_for_imputation/localRegion_%d' % (gse,region)
    filenames=os.listdir(meragefiledir)

    path = r'%s/Forimputation/local_methFeature_for_imputation/localRegion_%d/%s' % (gse,region,filenames[0])
    data = pd.read_csv(path,header=0,sep='\t')
    data[list(data)[2:]] = data[list(data)[2:]].astype('float16')
    data[list(data)[2:]] = data[list(data)[2:]].round(4)
    data = data[data.chrom.isin(chh[i1])]
    
    print('here2')    
    
    for i in range(1,len(filenames)):
        path = r'%s/%s' % (meragefiledir,filenames[i])
        df = pd.read_csv(path,header=0,sep='\t')
        df[list(df)[2:]] = df[list(df)[2:]].astype('float16')
        df[list(df)[2:]] = df[list(df)[2:]].round(4)
        df = df[df.chrom.isin(chh[i1])]
        data = pd.merge(data,df,how='outer',on=['chrom','location'])
    print('here3')        
    file_dir_1 = r'%s/region%d/corr' % (file_dir_Block,region)
    if not os.path.exists(file_dir_1):
        os.makedirs(file_dir_1)

    file_dir_2 = r'%s/region%d/methy' % (file_dir_Block,region)
    if not os.path.exists(file_dir_2):
        os.makedirs(file_dir_2)
    a =  list(data)   
    for i in range(len(cell_num)): 
        cell_r = []
        cell_meth = []
        for j in range(2,len(a)):
            if (a[j].split('_')[0] == cell_num[i]) or (a[j].split('_')[1] == cell_num[i]):
                if a[j].split('_')[-1] == 'r':
                    cell_r.append(a[j])
                if a[j].split('_')[-1] == 'methy':
                    cell_meth.append(a[j])            
        
        df_r =data[ ['chrom','location'] + cell_r ] 
        df_m = data[ ['chrom','location'] + cell_meth ] 
        
        df_r = df_r.dropna(axis=0,subset=list(df_r)[2:],how='all') 
        df_m = df_m.dropna(axis=0,subset=list(df_m)[2:],how='all') 
        
        df_r.to_csv(r'%s/%s_r.txt' % (file_dir_1,cell_num[i])  ,sep='\t',header=True,index=False)    
        df_m.to_csv(r'%s/%s_m.txt' % (file_dir_2,cell_num[i])  ,sep='\t',header=True,index=False)   
    return ('Block%d:Done!' % i1)
        
def unionfinalfile(file_dir,meragefiledir,filenames,filenames1,filenames2,i1,region,file_dir_1,file_dir_2):
	
    if i1+4 >= len(filenames1):
        i2 = len(filenames1)
    else:
        i2 = i1+4   
    print(i1,i2)

    for j in range(i1,i2):
        path = r'%s/%s/region%d/corr/%s' %(meragefiledir,filenames[0],region,filenames1[j])
        df_r = pd.read_csv(path,header=0,sep='\t')
        path = r'%s/%s/region%d/methy/%s' %(meragefiledir,filenames[0],region,filenames2[j])
        df_m = pd.read_csv(path,header=0,sep='\t')
             
        for i in range(1,len(filenames)):
            path = r'%s/%s/region%d/corr/%s' %(meragefiledir,filenames[i],region,filenames1[j])
            data = pd.read_csv(path,header=0,sep='\t')
            df_r = pd.concat([df_r,data])
            
            path = r'%s/%s/region%d/methy/%s' %(meragefiledir,filenames[i],region,filenames2[j])
            data = pd.read_csv(path,header=0,sep='\t')
            df_m = pd.concat([df_m,data])            
            
        df_r = df_r.drop_duplicates(['chrom','location'])

 
        df_r.to_csv(r'%s/%s' % (file_dir_1,filenames1[j]),sep='\t',header=True,index=False) 
        
        df_m = df_m.drop_duplicates(['chrom','location'])

        df_m.to_csv(r'%s/%s' % (file_dir_2,filenames2[j]),sep='\t',header=True,index=False)         
        
    return ('Block%d:Done!' % i1)
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
            data[list(data)[2:]] = data[list(data)[2:]].astype('float16')
            data[list(data)[2:]] = data[list(data)[2:]].round(4)
            if i==0:
                df = pd.merge(df,data,how='outer')
            else:
                df = pd.concat([df,data])
            #print (j,i,len(df))
        df = df.drop_duplicates(['chrom','location'])
        df.to_csv(r'%s/%s' % (file_dir,filenames1[j]),sep='\t',header=True,index=False)     
    return ('jincheng%d:Done!' % i1)
    
#get local_methFeature
def test(data,value,neighbor_region,file_dir,gse,bocknum,i1):    
    if i1+2 >= ( len(list(data))-1 ):
        i2 = len(list(data))-1 
    else:
        i2 = i1+2
    for i in range(i1,i2):
        print ('%s: start!' % list(data)[i])       
        #data_all = data[['chrom','location']] 
        data_all = pd.DataFrame(np.random.randn(0, 2), columns=['chrom','location'])
        #¶ÁÈ¡Ïà¹ØÐÔ±í¸ñ
        path = r'%s/for_impu_r/Block%d/%s_local_r.txt' % (gse,bocknum,list(data)[i]) 
        data_r_for_imputation = pd.read_csv(path,header=0,sep='\t')
        
        
        for k in range(i+1,len(list(data))):
            data1 = data[['%s' % list(data)[0],'%s' % list(data)[1],'%s' % list(data)[i],'%s' % list(data)[k]]]
            data4 = data_r_for_imputation[ ['%s' % list(data_r_for_imputation)[0],'%s' % list(data_r_for_imputation)[1],'%s_%s_r' % (list(data1)[2],list(data1)[3]) ] ]
            #
            data_chr = pd.DataFrame(np.random.randn(0, 4), columns=['chrom','location','%s_%s_r' % (list(data1)[2],list(data1)[3]),'%s_%s_methy' % (list(data1)[2],list(data1)[3])])
            for j in range(len(value)):                                                 
                #ÓÅ»¯-zl
                data2 = data1[data1['chrom'] ==  value[j]] 
                data2 = data2.dropna(subset=[list(data2)[3]])
                data2 = data2.sort_values(by='location') 
                
                serchtable = list(data2[list(data2)[2]])
                serchtable_m = list(data2[list(data2)[3]])
                serchtable_l = list(data2[list(data2)[1]])
                
                data2 = data2.dropna(subset=[list(data2)[2]]) 
                #print('here?')
                data2 = pd.merge(data2,data4,how='left',on=['chrom','location'])
                                
                data2_value =data2.values
                                                
                loca_for_imputation = []
                imputation_r = []
                imputation_m = []
                
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
                                loca_for_imputation.append(serchtable_l[l])
                                imputation_r.append( data2_value[ss[-1],-1]  )
                                imputation_m.append( serchtable_m[l] )    
                            
                data_for_impu = pd.DataFrame(np.random.randn(len(loca_for_imputation), 4), columns=['chrom','location','%s_%s_r' % (list(data1)[2],list(data1)[3]),'%s_%s_methy' % (list(data1)[2],list(data1)[3])])        
                data_for_impu['chrom'] =  value[j]                     
                data_for_impu['location'] =   loca_for_imputation                       
                data_for_impu['%s_%s_r' % (list(data1)[2],list(data1)[3])] =  imputation_r                      
                data_for_impu['%s_%s_methy' % (list(data1)[2],list(data1)[3])] =  imputation_m                        
                data_for_impu['%s_%s_methy' % (list(data1)[2],list(data1)[3])] = np.log2(data_for_impu['%s_%s_methy' % (list(data1)[2],list(data1)[3])]+1.01)
                data_for_impu['%s_%s_methy' % (list(data1)[2],list(data1)[3])] = data_for_impu['%s_%s_methy' % (list(data1)[2],list(data1)[3])].round(4)
                if len(data_for_impu) != 0:
                    data_chr = pd.merge(data_chr,data_for_impu,how='outer')               
            if len(data_chr) != 0:
                data_all = pd.merge(data_all,data_chr,how='outer',on=['chrom','location'])
                
            print ('%s-%s: Done!' % (list(data1)[2],list(data1)[3]))

        data_all.to_csv(r'%s/%s_local_methFeature.txt' %  (file_dir,list(data)[i]) ,sep='\t',header=True,index=False)  
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
    #gse = r'/home/tjx/SinglecellData/test_MPP_0702'
    #ff = r'GSE87197_mpp_d4.txt'    
    #neighbor_region = 10    
    #threshold =0.8
    
    file_dir = r'%s/Forimputation/local_methFeature' %  gse
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
    
    
    data = data[data.chrom.isin(ch)]  
    data = data.drop_duplicates(['chrom','location'])
     
    print('start!')

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
               
        with ProcessPoolExecutor() as pool:
            #´«µÝdata-·Ö¿é
            futures = [pool.submit(test,data1,value,neighbor_region,file_dir_Block,gse,k,i) for i in ll]
                           
            for j in as_completed(futures):
                print(j.result()) # ¶ÔÓ¦½ø³ÌÍê³ÉË³ÐòÊä³ö£¬Êä³ö2£¬3£¬5£¬6£¬10    
    
    elapsed = (time.clock() - start)
    print("Time used: %d s" % round(elapsed,4))
    del data
    del data1
    print("get local methFeature: Done!" )

####################################    
    #

    print("union local methFeature: Start!" )

    meragefiledir = r'%s' % file_dir
    filenames=os.listdir(meragefiledir)
    
    meragefiledir1 = r'%s/%s' % (file_dir,filenames[0])
    filenames1=os.listdir(meragefiledir1)

    file_dir = r'%s/Forimputation/local_methFeature_for_imputation/localRegion_%d' % (gse,neighbor_region)
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


    print("union large file")
####################################  
    region = neighbor_region 

    ll = range(len(chh))
    with ProcessPoolExecutor() as pool:
        futures = [pool.submit(unionlargefile,i,gse,region,chh,cell_num) for i in ll]
        for j in as_completed(futures):
            print(j.result())       
        
    print("union large file: start!") 
    print("-----------------------") 
                      
    ############################################################################       
    print("union final file: start!") 
    file_dir = r'%s/Forimputation/local_methFeature_cellbycell' % gse
    
    meragefiledir = r'%s' % file_dir
    filenames=os.listdir(meragefiledir)     
    
    meragefiledir1 = r'%s/%s/region%d/corr' % (file_dir,filenames[0],region)
    filenames1=os.listdir(meragefiledir1)
    meragefiledir2 = r'%s/%s/region%d/methy' % (file_dir,filenames[0],region)
    filenames2=os.listdir(meragefiledir2)    
    
    file_dir_1 = r'%s/Forimputation/local_methFeature_cellbycell/region%d/corr' % (gse,region)
    if not os.path.exists(file_dir_1):
        os.makedirs(file_dir_1)

    file_dir_2 = r'%s/Forimputation/local_methFeature_cellbycell/region%d/methy' % (gse,region)
    if not os.path.exists(file_dir_2):
        os.makedirs(file_dir_2)
        
    a = len(filenames1)
    ll = [0]
    k = 0
    for i in range(int(math.floor(a/4))):
        k += 4
        ll.append(k)
    if a%2 != 0:
        ll.append(k+1)
    print(ll)
    with ProcessPoolExecutor() as pool:
        futures = [pool.submit(unionfinalfile,file_dir,meragefiledir,filenames,filenames1,filenames2,i,region,file_dir_1,file_dir_2) for i in ll]
        for j in as_completed(futures):
            print(j.result())  
             
    
    print("union final file:Done!")    
    ############################################################################
    print("split: start!")
   
    for i in range(len(cell_num)):
        #corr
        path = r'%s/%s_r.txt' % (file_dir_1,cell_num[i])
        data_r = pd.read_csv(path,header=0,sep='\t')
        #methy
        path = r'%s/%s_m.txt' % (file_dir_2,cell_num[i])
        data_m = pd.read_csv(path,header=0,sep='\t')
        
        df_r = data_r[['chrom','location']]
    
        name = list(data_r)[2:]
                       
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
                
        file_dir = r'%s/Forimputation/region%d_localmatched_morethan08' % (gse,region)
        if not os.path.exists(file_dir):
            os.makedirs(file_dir)
        
        path = r'%s/%s.txt' % (file_dir,cell_num[i])
        
        df_r.to_csv(path, sep='\t', header=True, index=False)    
    print("split: Done!")
    
    
    
    
    
    
    
    
    
    
