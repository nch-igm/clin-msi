import os
import math  
import glob
import pickle
import logging

import shap
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable

pd.options.mode.chained_assignment = None  # default='warn'

## AUX FUNCS
def grab_marker_dict(df):
    split_col_set=[y for y in [x.split("_") for x in df.columns] if len(y)>1]
    coords=list(set([y[1] for y in split_col_set if ':' in y[1] and '-' in y[1]]))
    marker_dict={x:[y[0] for y in split_col_set if y[1]==x] for x in coords}
    return(marker_dict)
def grab_marker_int(x):
    if ('-' not in x and '+' not in x):
        return(0)
    else:
        return(int(x[1:]))
def apply_mod_to_dataframe(df,moddir):
    shapdict={}
    nrun=500
    ## DF ..combdat is 49, comb is 30
    #df=pd.read_csv(infile)
    j=0
    modfiles = glob.glob(os.path.join(moddir, "*pkl"))
    for modfile in modfiles:
        if (j % 50) ==0 :
            xx=1
            #logging.info(j)
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
        j+=1
    df['yprob']=df[['yprob_'+str(j+1) for j in range(nrun)]].sum(axis=1)/nrun
    df['ypred']=1*(df['yprob']>=.5)
    df['MSI_prediction']=df['ypred'].apply(lambda x: 'MSI' if x==1 else 'MSS')
    df['y']=df['yprob'].apply(lambda x: 0 if x<.5 else 1)
    return(df,shapdict,incols)

## GATHER SHAP DATA
def grab_shap_data(df,shapdict,incols):
    nrun=len(shapdict)
    shapdat=pd.concat(shapdict,axis=1)
    shapdat.columns=[x[1] for x in shapdat.columns]
    hf=df.copy()
    for x in incols:
        temp=shapdat[['shap_'+x+'_'+str(y+1) for y in range(nrun)]]
        hf['shap_'+x+'_mean']=temp.mean(axis=1)
        hf['shap_'+x+'_std']=temp.std(axis=1)
    marker_dict=grab_marker_dict(hf)
    bounddict={}
    bounddict[0]=[-10,0]
    bounddict[1]=[0,10]
    for ibin in bounddict.keys():
        sumset=range(bounddict[ibin][0],bounddict[ibin][1])
        for marker in marker_dict:
            marker_strings=[x for x in marker_dict[marker] if grab_marker_int(x) in sumset]
            hf['binnedshap_'+str(ibin)+'_'+marker+'_mean']=hf[['shap_'+x+'_'+marker+'_mean' for x in marker_strings]].mean(axis=1)
            hf['binnedshap_'+str(ibin)+'_'+marker+'_std']=hf[['shap_'+x+'_'+marker+'_mean' for x in marker_strings]].std(axis=1)
            hf['binnedval_'+str(ibin)+'_'+marker+'_mean']=hf[[x+'_'+marker for x in marker_strings]].mean(axis=1)
            hf['binnedval_'+str(ibin)+'_'+marker+'_std']=hf[[x+'_'+marker for x in marker_strings]].std(axis=1)
    bmdict={}
    for i in hf.index:
        logging.info(i)
        currec=hf.loc[i]
        mdict={}
        for marker in marker_dict:
            for ibin in bounddict.keys():
                bindesc=str(bounddict[ibin][0])+'...'+str(bounddict[ibin][1])
                keydesc=marker+"_"+bindesc
                locdict={}
                for stat in ['mean','std']:
                    for dattype in ['shap','val']:
                        locdict[dattype+'_'+stat]=currec['binned'+dattype+'_'+str(ibin)+'_'+marker+'_'+stat]
                mdict[keydesc]=locdict
        bmdict[i]=mdict.copy()
    return(hf,bmdict)

## BUILD SHAP PLOTS
def build_n_save_shap_plot(currec,mdict,shap_plot_file):
    sdat=pd.DataFrame(mdict).transpose().reset_index().rename(columns={'index':'markerbin'})
    sdat['shapabs']=sdat['shap_mean'].apply(abs)
    sdat=sdat.sort_values('shapabs',ascending=True)
    sdat['val_mean_norm']=sdat['val_mean'].apply(lambda x: math.sqrt(x-sdat.val_mean.min()+.00001))
    #sdat['val_mean_norm']=sdat['val_mean'].apply(lambda x: min(x,1))
    sdat['val_mean_discrete']=sdat['val_mean'].apply(lambda x: -1 if x <= 0 else 0 if x<=3 else 1)
    sdat['color']=sdat['val_mean_discrete'].apply(lambda x: 'green' if x==-1 else 'orange' if x==0 else 'red')
    if ((sdat.shapabs>0.005).mean() > 0):
        sdat=sdat[sdat.shapabs > 0.005]
    else:
        sdat=sdat.head(10)
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
    plt.title('; '.join([currec.samp,'msi_conflevel='+str(round(currec.yprob,3))]))
    plt.savefig(shap_plot_file,bbox_inches='tight')


## MAIN RUNNER
def apply_model(infile,moddir,outfile,shap_plot_dir=None):
    df=pd.read_csv(infile)
    sampcol='SAMPLE_NAME'
    if sampcol not in df:
        df[sampcol]=['SAMPLE_' + str(i+1) for i in range(len(df))]
    df=df.rename(columns={sampcol:'samp'})
    dfnew,shapdict,curfeats=apply_mod_to_dataframe(df,moddir)
    dfnew[['samp','yprob']].to_csv(outfile,index=False)
    ## GRAB  SHAP DATA
    testdat,bmdict=grab_shap_data(dfnew,shapdict,curfeats)
    ## EXPORT SHAP PLOTS
    for i in range(len(testdat)):
        shap_outfile=shap_plot_dir+"/shap_plot_"+str(i+1)+'.png'
        shaprec=testdat.iloc[i]
        shapdict_loc=bmdict[i]
        build_n_save_shap_plot(shaprec,shapdict_loc,shap_outfile)
