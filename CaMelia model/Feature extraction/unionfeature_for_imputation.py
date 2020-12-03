# -*- coding: utf-8 -*-
from __future__ import division
from sys import argv
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')


###########################################
if __name__ == '__main__': 
    #All data located in the same directory
    DataPath = r'%s' % argv[1]
    #Input data     
    InputDataName = '%s' % argv[2]  
    #neighboring range
    neighbor_region = int(argv[3])
    
   
   
    gse = DataPath
    ff  = InputDataName	    
    region = neighbor_region
    
    
    name = ff.split('.')[0].split('_')[-1]
    
    
    path = r'%s/%s' % (gse,ff)
    data = pd.read_csv(path,header=0,sep='\t')
    if list(data)[0] != 'chrom':
        del data['%s' % list(data)[0]]

    data[list(data)[2:]] = data[list(data)[2:]].astype('float16')   
    cell_num = list(data)[2:( len(list(data))-1 )]
    
                
    file_dir = r'%s/Available_Imputation_dataset/region%d' % (gse,region)
    if not os.path.exists(file_dir):
        os.makedirs(file_dir)   
        
    for i in range(len(cell_num)):            
        #local_methFeature
        path = r'%s/Forimputation/region10_localmatched_morethan08/%s.txt' % (gse,cell_num[i])
        local_methFeature = pd.read_csv(path,header=0,sep='\t')
        local_methFeature = local_methFeature.rename(columns={'aver_meth':'local_methFeature'})
        local_methFeature = local_methFeature[['chrom','location','local_methFeature']]
        local_methFeature[list(local_methFeature)[2:]] = local_methFeature[list(local_methFeature)[2:]].astype('float16')

        #neighbor_methFeature
        path = r'%s/Forimputation/neighbor_methFeature_%d/localRegion_%s/%s_neighbor_methFeature.txt' % (gse,region,region,cell_num[i])
        neighbor_methFeature = pd.read_csv(path,header=0,sep='\t')
        neighbor_methFeature[list(neighbor_methFeature)[2:]] = neighbor_methFeature[list(neighbor_methFeature)[2:]].astype('float16')

        #merge-[neighbor,local,global]
        data_all = pd.merge(data[['chrom','location','%s' % cell_num[i] ]],local_methFeature)
        data_all = pd.merge(data_all,neighbor_methFeature,how='inner',on=['chrom','location'])
        
        data_all.to_csv(r'%s/%s.txt' % (file_dir,cell_num[i]),sep='\t',header=True,index=False)    
        del data_all