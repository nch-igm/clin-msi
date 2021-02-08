#!/usr/bin/env python
# coding: utf-8

# In[12]:


## IMPORT
import pandas as pd
import numpy as np
import pickle
import xgboost as xgb


# In[ ]:





# In[56]:


## FUNCTIONS
## INTERNAL TRAINING FUNCTION
def modfit(df,curfeats):
    X=df[curfeats]
    Y=df['y']
    model = xgb.XGBClassifier(max_depth=1,
                              colsample_bytree=0.001,
                              n_estimators=500,
                              eval_metric='error',
                              use_label_encoder=False ## to prevent a warning
                             )
    model.fit(X,Y)
    return model

## MAIN FUNCTION
def train_models(df,moddir):
    curfeats=[x for x in df if x != 'y']
    moddict={}
    nrun=500
    modfiles=[moddir+'/xgb_'+str(j+1)+'.pkl' for j in range(nrun)]
    for j in range(nrun):
        if (j % 10 == 9):
            ## DINKY PROGRESS REPORT
            print("Training model "+str(j+1))
        traininds_loc=df.sample(frac=0.8).index
        mymod=modfit(df.loc[traininds_loc],curfeats)
        pkl_connect = open(modfiles[j], 'wb')
        pickle.dump(mymod, pkl_connect)
        pkl_connect.close()
    return(modfiles)


# In[ ]:





# In[57]:


## EXE

#infile='/igm/home/jbg001/git/msi/dat/sample_patricks_msi_df.csv'
#rawdat=pd.read_csv(infile)
#moddir_temp='/igm/home/jbg001/git/msi/tm/temp_mods_020521'
#modfiles=train_models(rawdat,moddir_temp)


# In[ ]:





# In[59]:


## OUTPUT
#modfiles[:5]


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




