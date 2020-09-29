# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 09:16:08 2020

@author: tjx
"""
from __future__ import division
from sys import argv
import pandas as pd
import numpy as np
import os


from numba import jit

import matplotlib.pyplot as plt 



#»®·ÖÊý¾Ý¼¯¡¢¼ÆËãÖ¸±ê
from sklearn import metrics
#5ÕÛ-½»²æÑéÖ¤
from sklearn.model_selection import KFold

import catboost as cb

import warnings
warnings.filterwarnings('ignore')

##################################
#¼ÆËãÖ¸±ê
@jit
def auc(m, train, test):
    return (metrics.roc_auc_score(y_train,m.predict_proba(X_train)[:,1]),metrics.roc_auc_score(y_test,m.predict_proba(X_test)[:,1]))
@jit   
def acc(m, train, test):
    a = m.predict_proba(X_train)[:,1]
    a = np.int64(a>=0.5)
    b = m.predict_proba(X_test)[:,1]
    b = np.int64(b>=0.5)    
    return (metrics.accuracy_score(y_train,a),metrics.accuracy_score(y_test,b))

@jit
def spsemcc(m,X_train,y_true):                   
    TP = 0
    FP = 0
    FN = 0
    TN = 0 
    
    y_pre = m.predict_proba(X_train)[:,1]
    y_pre = np.int64(y_pre>0.5)
    
    for j in range(len(y_true)):
        if y_true[j] == 1:
            if y_pre[j] == 1:
                TP += 1
            else:
                FN += 1
        if y_true[j] == 0:
            if y_pre[j] == 0:
                TN += 1
            else:
                FP += 1  

    sp =  (1 - FP / (TN + FP))  
    se = TP/(TP + FN) 
    mcc = (TP*TN-FP*FN)/(((TN+FN)*(TN+FP)*(TP+FN)*(TP+FP))**0.5)
         
    return 'sp: ',sp,'se: ',se,'mcc: ',mcc
################################################
if __name__ == '__main__':
    
    #All data located in the same directory
    DataPath = r'%s' % argv[1]
    #Input data     
    InputDataName = '%s' % argv[2]  
    #neighboring range
    t_t = '%s' % argv[3]
   
    
    gse = DataPath
    #gse = r'/home/tjx/SinglecellData/test_MPP_0702'
    
    meragefiledir_out = r'%s/Available_Train_dataset/region10' % gse
    filenames_out=os.listdir(meragefiledir_out)
    
    meragefiledir_impu = r'%s/Available_Imputation_dataset/region10' % gse
    filenames_impu=os.listdir(meragefiledir_impu)
    
    
    file_dir_png = r'%s/model_feature_importance' % gse
    if not os.path.exists(file_dir_png):
        os.makedirs(file_dir_png) 
        
    file_dir_pre = r'%s/imputation_data' % gse
    if not os.path.exists(file_dir_pre):
        os.makedirs(file_dir_pre) 
        
    file_dir_model = r'%s/model_save' % gse
    if not os.path.exists(file_dir_model):
        os.makedirs(file_dir_model)    
    
    ACC =[]
    AUC = []
    SP =[]
    SE =[]
    MCC =[]
    sample =[]    
    
    for i in range(len(filenames_out)):          
        #Train data
        path = r'%s/%s' % (meragefiledir_out,filenames_out[i])
        dataset_train = pd.read_csv(path,header=0,sep='\t')
        dataset_train = dataset_train.dropna(axis=1,how='all')   
        dataset_train = dataset_train.dropna(axis=0,how='any')
        
        #Imputation data
        path = r'%s/%s' % (meragefiledir_impu,filenames_out[i])
        dataset_impu = pd.read_csv(path,header=0,sep='\t')
        dataset_impu = dataset_impu.dropna(axis=1,how='all')   
        dataset_impu = dataset_impu.dropna(axis=0,how='any') 
    
        ####
        dataset_ = dataset_train.values
        dataset_im = dataset_impu.values
              
        # split data into X and y
        X = dataset_[:,3:]
        Y = dataset_[:,2]
        Y = np.int64(Y>=0.5)
        
        
        seed = 514
        kf = KFold(n_splits=5,random_state=seed,shuffle=True)

        
        kf_ACC = []
        kf_AUC = []
        kf_SP = []
        kf_SE = []
        kf_MCC = []
        for train_index, test_index in kf.split(X, Y):
            X_train, X_test = X[train_index,:], X[test_index,:]
            y_train, y_test = Y[train_index], Y[test_index] 
            
            #ÑµÁ·
            model = cb.CatBoostClassifier(random_state=seed,learning_rate=0.1,max_depth=7,verbose=1000,eval_metric='AUC',task_type=t_t)
            model.fit(X_train,y_train)
                    
            #ÆÀ¹ÀÐÔÄÜ    
            kf_ACC.append( acc(model, X_train, X_test)[1] )
            kf_AUC.append( auc(model, X_train, X_test)[1] )
            
            a =spsemcc(model,X_test,y_test)
            
            kf_SP.append( a[1] )
            kf_SE.append( a[3] )
            kf_MCC.append( a[5] ) 
        
        print ('%s: 5-fold Crossvalidation Done!' % filenames_out[i].split('.')[0])
        ACC.append( round( np.mean(kf_ACC),4) ) 
        AUC.append( round(np.mean(kf_AUC),4) )
        SP.append(round( np.mean(kf_SP),4) )
        SE.append(round( np.mean(kf_SE),4) )
        MCC.append(round( np.mean(kf_MCC),4) ) 
        #all data fit
        model = cb.CatBoostClassifier(learning_rate=0.1,max_depth=7,verbose=1000,eval_metric='AUC',task_type=t_t)
        model.fit(X,Y)   
        #save model
        model.save_model(r'%s/%s.model' % (file_dir_model,filenames_out[i].split('.')[0]))
        print('%s: training Done!' % filenames_out[i].split('.')[0])
        print('----------------')                
        print('%s: imputation start!' % filenames_out[i].split('.')[0])    
        #ÖØÒªÐÔ
        #print(model.feature_importances_)        
        plt.figure(figsize=(10,6))
        ax = plt.bar(range(len(model.feature_importances_)), model.feature_importances_)
        plt.rcParams['savefig.dpi'] = 300 
        plt.savefig(r'%s/%s_feature_importance.png' % (file_dir_png,filenames_out[i].split('.')[0]), dpi=300)
        #pre 
        prediction = model.predict_proba(dataset_im[:,2:])[:,1]
        pre = []
        for m in range(len(prediction)):
            if prediction[m] >= 0.6:
                pre.append(1)
            elif prediction[m] <= 0.4: 
                pre.append(0)
            else: 
                pre.append(np.nan)
            
        data_pre = pd.DataFrame(np.random.randn(0, 3), columns=['chrom','location','pre_meth']) 
        data_pre['chrom'] = dataset_im[:,0]
        data_pre['location'] = dataset_im[:,1]  
        data_pre['pre_meth'] = pre
        #save
        data_pre.to_csv(r'%s/%s.txt' % (file_dir_pre,filenames_out[i].split('.')[0]),sep='\t',header=True,index=False)     
        print('%s: imputation Done!' % filenames_out[i].split('.')[0])
            
    cell = []
    for i in range(len(filenames_out)):
        cell.append( filenames_out[i].split('.')[0] )
           
    data = pd.DataFrame(np.random.randn(0,6), columns=['cell','SP','SE','ACC','AUC','MCC'])    
    data['cell'] = cell
    data['SP'] = SP
    data['SE'] = SE
    data['ACC'] = ACC
    data['AUC'] = AUC
    data['MCC'] = MCC   
    data.to_csv(r'%s/5fold_crossvalidation_catboost.csv' % gse  ,sep=',',header=True,index=False)     