#!/usr/bin/env python
# coding: utf-8

# In[ ]:


## IMPORT
import pandas as pd
import matplotlib
import numpy as np
import random
import imp
import sys
import pickle
import matplotlib.patches as mpatches

#from pyliftover import LiftOver
from matplotlib.cm import ScalarMappable
import math  
import string
import f
import shap
import matplotlib.pyplot as plt
pd.set_option("display.max_rows", 100)
pd.set_option('display.max_columns', 500)
pd.options.mode.chained_assignment = None  # default='warn'


if 1==1:
    infile=sys.argv[sys.argv.index('--infile')+1]
    outfile=sys.argv[sys.argv.index('--outfile')+1]
    plotfile=sys.argv[sys.argv.index('--plotfile')+1]
    moddir=sys.argv[sys.argv.index('--moddir')+1]
else:
    infile='/igm/home/jbg001/git/msi/dat/exp8_endcapped_zscore_normalized.csv'
    outfile='/igm/home/jbg001/git/msi/toyout.csv'
    plotfile='/igm/home/jbg001/git/msi/toyout.csv'
    moddir='/igm/home/jbg001/git/msi/mod_071019_frozen'


# In[ ]:





# In[ ]:





# In[ ]:


### MODEL - TEST
imp.reload(f)
shapdict={}
nrun=500
## DF ..combdat is 49, comb is 30
#df=pd.read_csv(infile)
df=pd.read_csv(infile).head(1)

incols=f.grab_feats(df)
for j in range(nrun):
    if (j % 50) ==0 :
        xx=1
        #print(j)
    modfile=moddir+'/xgb_'+str(j+1)+'.pkl'
    pkl_connect = open(modfile, 'rb')
    mymod=pickle.load(pkl_connect)
    pkl_connect.close()
    incols=mymod.get_booster().feature_names
    df.loc[:,'yprob_'+str(j+1)]=mymod.predict_proba(df[incols])[:,1]
    shapflag=1
    if shapflag==1:
        explainer=shap.TreeExplainer(mymod)
        shapdat_loc=pd.DataFrame(explainer.shap_values(df[incols]),columns=['shap_'+x+'_' + str(j+1) for x in incols])
        shapdict[j]=shapdat_loc
df['yprob']=df[['yprob_'+str(j+1) for j in range(nrun)]].sum(axis=1)/nrun
df['ypred']=1*(df['yprob']>=.5)
df['MSI_prediction']=df['ypred'].apply(lambda x: 'MSI' if x==1 else 'MSS')


# In[ ]:


# EXPORT 
df=df.rename(columns={df.columns[0]:'samp'})
df['y']=df['yprob'].apply(lambda x: 0 if x<.5 else 1)
df[['samp','yprob']].to_csv(outfile,index=False)


# In[ ]:





# In[ ]:


## SHAP
shapdat=pd.concat(shapdict,axis=1)
shapdat.columns=[x[1] for x in shapdat.columns]
hf=df.copy()
for x in incols:
    temp=shapdat[['shap_'+x+'_'+str(y+1) for y in range(nrun)]]
    hf['shap_'+x+'_mean']=temp.mean(axis=1)
    hf['shap_'+x+'_std']=temp.std(axis=1)
markers=[x[2:] for x in hf.columns if x[:2]=='0_']
bounddict={}
bounddict[0]=[0,10]
bounddict[1]=[10,20]
bounddict[2]=[20,30]
for ibin in bounddict.keys():
    sumset=range(bounddict[ibin][0],bounddict[ibin][1])
    for marker in markers:
        hf['binnedshap_'+str(ibin)+'_'+marker+'_mean']=hf[['shap_'+str(x)+'_'+marker+'_mean' for x in sumset]].mean(axis=1)
        hf['binnedshap_'+str(ibin)+'_'+marker+'_std']=hf[['shap_'+str(x)+'_'+marker+'_mean' for x in sumset]].std(axis=1)
        hf['binnedval_'+str(ibin)+'_'+marker+'_mean']=hf[[str(x)+'_'+marker for x in sumset]].mean(axis=1)
        hf['binnedval_'+str(ibin)+'_'+marker+'_std']=hf[[str(x)+'_'+marker for x in sumset]].std(axis=1)
bmdict={}
for i in hf.index:
    currec=hf.loc[i]
    mdict={}
    for marker in markers:
        for ibin in bounddict.keys():
            bindesc=str(bounddict[ibin][0])+'...'+str(bounddict[ibin][1])
            keydesc=marker+"_"+bindesc
            locdict={}
            for stat in ['mean','std']:
                for dattype in ['shap','val']:
                    locdict[dattype+'_'+stat]=currec['binned'+dattype+'_'+str(ibin)+'_'+marker+'_'+stat]
            mdict[keydesc]=locdict
    bmdict[i]=mdict


# In[ ]:





# In[ ]:


## PLOT
imp.reload(f)
testdat=hf
cv=0
for i in testdat.index:
    currec=testdat.loc[i]
    mdict=bmdict[i]
    sdat=pd.DataFrame(mdict).transpose().reset_index().rename(columns={'index':'markerbin'})
    sdat['shapabs']=sdat['shap_mean'].apply(abs)
    sdat=sdat.sort_values('shapabs',ascending=True)
    sdat['val_mean_norm']=sdat['val_mean'].apply(lambda x: math.sqrt(x-sdat.val_mean.min()+.00001))
    #sdat['val_mean_norm']=sdat['val_mean'].apply(lambda x: min(x,1))
    sdat['val_mean_discrete']=sdat['val_mean'].apply(lambda x: -1 if x <= 0 else 0 if x<=3 else 1)
    sdat['color']=sdat['val_mean_discrete'].apply(lambda x: 'green' if x==-1 else 'orange' if x==0 else 'red')
    sdat=sdat[sdat.shapabs > 0.005]
    colguide=sdat.val_mean_discrete
    #sdat=sdat[sdat.val_mean < 2]
    tempdict={};tempdict[0]='MSS';tempdict[1]='MSI'
    data_x=sdat.markerbin
    data_y=sdat.shap_mean
    data_color=sdat.val_mean_norm*1.0
    my_cmap=plt.cm.get_cmap('RdYlGn_r')
    colors=my_cmap(data_color)

    fig, ax = plt.subplots()
    rects=ax.barh(data_x, data_y, color=colors)
    sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(0,max(data_color)))
    sm.set_array([])
    cbar = plt.colorbar(sm)
    cbar.set_label('Normalized read-count score', rotation=270,labelpad=25)
    plt.title('; '.join([currec.samp,tempdict[currec.y], 'msi_conflevel='+str(round(currec.yprob,3))]))
    plt.savefig(plotfile,bbox_inches='tight')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




