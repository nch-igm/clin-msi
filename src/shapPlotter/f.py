#!/usr/bin/env python
# coding: utf-8

# In[1]:


## IMPORT
import pandas as pd
import matplotlib
import numpy as np
import random
import imp
#from pyliftover import LiftOver
import math
from sklearn import svm
import string
import xgboost as xgb
from sklearn.linear_model import LogisticRegression

pd.set_option("display.max_rows", 100)
pd.set_option('display.max_columns', 500)

reffile="/igm/apps/genomes/Homo_sapiens/human_g1k_v37_decoy/human_g1k_v37_decoy.fasta"
basedir='/igm/home/jbg001/git/msi'
datdir=basedir+'/dat'
rawdir=datdir+'/raw_excel'


# In[ ]:





# In[2]:


## READ

def grab_rawfiles():
    varlist=['OPTEC_Exp' + str(x) +'_Variants_ALL' for x in [4,5]]
    otherlist=['OPTEC_Exp5_Testing_MSI']+['optec_msi_hg_countv3_updated_'+x for x in ['markertables','summarysheet']
                                         ]+['OPTEC_Exp'+str(i)+'_MSI_TRUTH' for i in [4,5]]+['OPTEC_Exp5_rawtable']
    exp3list=['Exp3_'+x for x in ['data_21loci','data_21loci_oldpreds','Truth']]
    totlist=varlist+otherlist+exp3list
    filelist=[rawdir+'/'+x+'.tsv' for x in totlist]

    filelist
    locdict={}
    for x in filelist:
        locdict[x]=pd.read_csv(x,sep="\t")
    return(locdict)
def build_rawdat(rawdict):
    pretraindat=rawdict['/igm/home/jbg001/git/msi/dat/raw_excel/optec_msi_hg_countv3_updated_summarysheet.tsv']
    truthdat=rawdict['/igm/home/jbg001/git/msi/dat/raw_excel/OPTEC_Exp4_MSI_TRUTH.tsv']
    traindat=pretraindat.drop(columns='nplus_type').set_index('marker').transpose().reset_index().rename(columns={'index':'samp'})
    truthdat['samp']=truthdat['ID'].apply(lambda x:x.replace(".","-"))
    truedat=pd.merge(right=truthdat,left=traindat,on=['samp'],how='inner')
    pretestdat=rawdict['/igm/home/jbg001/git/msi/dat/raw_excel/OPTEC_Exp5_Testing_MSI.tsv']
    truthdat_new=rawdict['/igm/home/jbg001/git/msi/dat/raw_excel/OPTEC_Exp5_MSI_TRUTH.tsv']
    return(truedat,pretestdat,truthdat_new)


# In[3]:


## FORMAT
def grab_rawdat_new():
    trainfile=datdir+'/msi_training_loci_expanded_ar.csv'
    testfile=datdir+'/test_matrix_axr.csv'
    test2file=datdir+'/NA_70_test_marker_matrix_z-normalized.csv'
    pretraindat=pd.read_csv(trainfile).assign(expnum=4)
    valdat=pd.read_csv(testfile).assign(expnum=5)
    ambsamps=['Exp5A_D2_1045_normal_tissue','Exp5A_C3_1056_tumor_tissue']
    valdat.loc[valdat.samp.isin(ambsamps),'y']=-1
    val2dat=pd.read_csv(test2file).assign(expnum=3)
    val2dat['y']=val2dat['TRUTH'].apply(lambda x: 1 if x=='MSI' else 0 if x=='MSS' else -1)
    val2dat=val2dat.drop(columns='TRUTH')
    curdat=pd.concat([pretraindat,valdat,val2dat]).reset_index(drop=True)
    rawfeats=[x for x in curdat.columns if x not in ['samp','y','valflag']]
    return(curdat,rawfeats)

def grab_rawdat():
    trainfile=datdir+'/msi_training_loci_expanded_ar.csv'
    testfile=datdir+'/test_matrix_axr.csv'
    pretraindat=pd.read_csv(trainfile)
    valdat=pd.read_csv(testfile)
    curdat=pd.concat([pretraindat.assign(expnum=4),valdat.assign(expnum=5)]).reset_index(drop=True)
    rawfeats=[x for x in curdat.columns if x not in ['samp','y','valflag']]
    return(curdat,rawfeats)
def grab_varhouse():
    vhdict={}
    vhdict['pretrain']=pd.read_csv('/igm/home/jbg001/git/msi/dat/raw_excel/OPTEC_Exp4_Variants_ALL.tsv',sep="\t")
    vhdict['val']=pd.read_csv('/igm/home/jbg001/git/msi/dat/raw_excel/OPTEC_Exp5_Variants_ALL.tsv',sep="\t")
    tempcols=[x for x in vhdict['pretrain'].columns if x in vhdict['val'].columns]
    outdat=pd.concat([vhdict['pretrain'][tempcols],vhdict['val'][tempcols]])
    outdat['samp']=outdat['Sample'].apply(lambda x: x.replace('.tsv',''))
    return(outdat)
def grab_tables():
    listdict={'OPTEC_Exp5_rawtable':'test','optec_msi_hg_countv3_updated_markertables':'train'}
    locdict={}
    for x in listdict.keys():
        curf=rawdir+'/'+x+'.tsv'
        locdict[listdict[x]]=pd.read_csv(curf,sep="\t",header=None)
    truthdict={}
    truthdict['train']=pd.read_csv('/igm/home/jbg001/git/msi/dat/raw_excel/OPTEC_Exp4_MSI_TRUTH.tsv',sep="\t")
    truthdict['test']=pd.read_csv('/igm/home/jbg001/git/msi/dat/raw_excel/OPTEC_Exp5_MSI_TRUTH.tsv',sep="\t")
    vhdict={}
    vhdict['train']=pd.read_csv('/igm/home/jbg001/git/msi/dat/raw_excel/OPTEC_Exp4_Variants_ALL.tsv',sep="\t")
    vhdict['test']=pd.read_csv('/igm/home/jbg001/git/msi/dat/raw_excel/OPTEC_Exp5_Variants_ALL.tsv',sep="\t")
    return(locdict,truthdict,vhdict)
def format_table(rawdat):
    tabdict={}
    for i in range(len(rawdat)):
        if rawdat[1].isnull()[i]:
            j=i+1
            while (j<(len(rawdat)-1)) & (~rawdat[1].isnull()[j]):
                j=j+1
            if j==(len(rawdat)-1):
                j=j+1
            marker=rawdat[0][i]
            tempdat=rawdat[(i+2):j].copy()
            tempdat.columns=rawdat.iloc[1,:]
            xdat=tempdat.transpose()
            xdat.columns=['run_'+str(k+1) for k in range(30)]+['ave_loose','ave_tight']
            ydat=xdat.reset_index().rename(columns={1:'presamp'})
            tabdict[marker]=ydat
    tabdat=pd.concat(tabdict).reset_index().drop(columns='level_1').rename(columns={'level_0':'marker'}).pivot(index='presamp',columns='marker').reset_index()
    newdat=tabdat.copy()
    newdat.columns=['_'.join(col) for col in newdat.columns]
    return(newdat.rename(columns={'presamp_':'presamp'}))
def add_yvals(featdat,truthdat):
    newdat=pd.merge(left=truthdat.assign(flag1=1),right=featdat.assign(flag2=1),on='presamp',how='outer')
    newdat['Tumor_Type']=newdat['Tumor_Type'].fillna("X")
    newdat['y']=0
    newdat.loc[newdat['Tumor_Type']=="MSI",'y']=1
    return(newdat)
def format_testdat(estdat_pre):
    estdat=estdat_pre.copy()
    estdat['samp']=estdat.presamp.apply(lambda x: x.split("_")[2])
    estdat.loc[estdat.presamp.str.contains('HeLa'),'samp']=estdat[estdat.presamp.str.contains('HeLa')]['presamp'].apply(
        lambda x: 'HeLa_control_'+x.split("_")[4])
    estdat.loc[estdat.presamp.str.contains('G4_NA12878'),'samp']='NA12878_G4'
    estdat.loc[estdat.presamp.str.contains('D4_NA12878'),'samp']='NA12878_D4'
    #testdat=estdat[((estdat.presamp.str.contains("tumor")) | (estdat.presamp.str.contains("NA12878")))]
    testdat=estdat
    testdat['normflag']=testdat['presamp'].apply(lambda x: 1*("orma" in x))
    return(testdat)
def format_truthdat_test(temp):
    pretruthdat=temp.copy()
    print(len(pretruthdat))
    truthdat=pretruthdat[(~pretruthdat['NCH MSI RESULTS'].isnull())]
    print(len(truthdat))
    truthdat['samp']=truthdat['OPTEC Patient Code ']
    #truthdat=truthdat[~truthdat['Paul Goodfellow Sample ID'].str.contains("ormal")]
    truthdat['goodN']=truthdat['Paul Goodfellow Sample ID'].apply(lambda x: x[len(x)-1]=='N')
    #truthdat=truthdat[~truthdat['goodN']].drop(columns='goodN')
    truthdat=truthdat.drop(columns='goodN')
    truthdat['normflag']=truthdat['Paul Goodfellow Sample ID'].apply(lambda x: 1*(("norm" in x) or (x[len(x)-1]=='N')))
    return(truthdat)
def add_yvals_testdat(pretestdat,pretruthdat):
    testdat=pretestdat.copy()
    truthdat=pretruthdat.copy()
    testdat=testdat[~(testdat['presamp'].str.contains("Water"))]
    newdat=pd.merge(left=truthdat.assign(flag1=1),right=testdat.assign(flag2=1),on=['samp','normflag'],how='inner')
    newdat['y']=0
    newdat['MSI_status']=newdat['MSI_status'].fillna("X")
    newdat.loc[newdat['MSI_status'].str.contains("MSI"),'y']=1
    newdat.loc[newdat['MSI_status'].str.contains("MSS"),'y']=0
    newdat.loc[(newdat['MSI_status'].str.contains("MSI")) & ((newdat['MSI_status'].str.contains("MSS"))) ,'y']=-1
    ## Assign yvals

    return(newdat)
def add_y(curdat_pre):
    curdat=curdat_pre.copy()
    curdat.loc[curdat['Tumor_Type'].isnull(),'Tumor_Type']='none'
    curdat['y']=curdat['Tumor_Type'].apply(lambda x: 1 if "MSI" in x else 0 if "MSS" in x else -1)
    curdat.loc[curdat['Tumor_Type'].contains("quivocal"),'y']=-1
    curdat.loc[(curdat['Tumor_Type'].str.contains("MSI")) & (curdat['Tumor_Type'].str.contains("MSS"))  ,'y']=-1
    curdat.loc[(curdat['presamp'].str.contains("orma")) ,'y']=0
    return(curdat)
def format_groupcols(tempdat):
    tempdat['ClinVar Significance']=tempdat['ClinVar Significance'].fillna("none")
    tempdat['Impact']=tempdat['Impact'].fillna("none")

    tempdat['clinvar_score']=tempdat['ClinVar Significance'].apply(
        lambda x: 3 if "athogenic" in x else 2 if 'onflicting' in x else 1 if 'enign' in x else 0)
    tempdat['impact_score']=tempdat['Impact'].apply(
        lambda x: 3 if x=="HIGH" else 2 if x=="MODERATE" else 1 if x=="LOW" else 0)
    return(tempdat)
def distill_vhdat(df):
    maxscores=['impact_score','clinvar_score','REVEL Score','AF PopMax gnomAD EX','Max Polyphen2 HVAR score',
              'Max_Polyphen2 HDIV score','MetaLR_score','Homozygous Count gnomAD EX ']
    aggdict={x:'max' for x in maxscores}
    heydat=df.groupby(['samp','Gene']).agg(aggdict)
    eydat=heydat.reset_index().pivot(index='samp',columns='Gene').reset_index()
    eydat.columns=['_'.join(col) for col in eydat.columns]
    eydat=eydat.rename(columns={'samp_':'samp'})
    return(eydat)
def group_vardat(tempdat):
    tempdat['presamp']=tempdat['Sample'].apply(lambda x: x.replace('.tsv',''))
    tempdat['ClinVar Significance']=tempdat['ClinVar Significance'].fillna("none")
    tempdat['Impact']=tempdat['Impact'].fillna("none")

    tempdat['clinvar_score']=tempdat['ClinVar Significance'].apply(
        lambda x: 3 if "athogenic" in x else 2 if 'onflicting' in x else 1 if 'enign' in x else 0)
    tempdat['impact_score']=tempdat['Impact'].apply(
        lambda x: 3 if x=="HIGH" else 2 if x=="MODERATE" else 1 if x=="LOW" else 0)
    maxscores=['impact_score','clinvar_score','REVEL Score','AF PopMax gnomAD EX','Max Polyphen2 HVAR score',
              'Max_Polyphen2 HDIV score','MetaLR_score','Homozygous Count gnomAD EX ']
    aggdict={x:'max' for x in maxscores}
    heydat=tempdat.groupby(['presamp','Gene']).agg(aggdict)
    eydat=heydat.reset_index().pivot(index='presamp',columns='Gene').reset_index()
    eydat.columns=['_'.join(col) for col in eydat.columns]
    eydat=eydat.rename(columns={'presamp_':'presamp'})
    #newdat=pd.merge(left=ydat.assign(flag3=1),right=eydat.assign(flag4=1),on='presamp',how='outer')
    return(eydat)
def process_vardat(tempdat,ydat):
    tempdat['presamp']=tempdat['Sample'].apply(lambda x: x.replace('.tsv',''))
    tempdat['ClinVar Significance']=tempdat['ClinVar Significance'].fillna("none")
    tempdat['Impact']=tempdat['Impact'].fillna("none")

    tempdat['clinvar_score']=tempdat['ClinVar Significance'].apply(
        lambda x: 3 if "athogenic" in x else 2 if 'onflicting' in x else 1 if 'enign' in x else 0)
    tempdat['impact_score']=tempdat['Impact'].apply(
        lambda x: 3 if x=="HIGH" else 2 if x=="MODERATE" else 1 if x=="LOW" else 0)
    maxscores=['impact_score','clinvar_score','REVEL Score','AF PopMax gnomAD EX','Max Polyphen2 HVAR score',
              'Max_Polyphen2 HDIV score','MetaLR_score','Homozygous Count gnomAD EX ']
    aggdict={x:'max' for x in maxscores}
    heydat=tempdat.groupby(['presamp','Gene']).agg(aggdict)
    eydat=heydat.reset_index().pivot(index='presamp',columns='Gene').reset_index()
    eydat.columns=['_'.join(col) for col in eydat.columns]
    eydat=eydat.rename(columns={'presamp_':'presamp'})
    newdat=pd.merge(left=ydat.assign(flag3=1),right=eydat.assign(flag4=1),on='presamp',how='outer')
    return(newdat)
def process_vardat_test(tempdat,ydat):
    tempdat['presamp']=tempdat['Sample'].apply(lambda x: x.replace('.tsv',''))
    tempdat['normflag']=tempdat['Sample'].apply(lambda x: 1*('orma' in x))
    tempdat['samp']=tempdat['presamp']+'_'+tempdat['normflag'].apply(str)
    ydat['samp']=ydat['presamp']+'_'+ydat['normflag'].apply(str)
    tempdat['ClinVar Significance']=tempdat['ClinVar Significance'].fillna("none")
    tempdat['Impact']=tempdat['Impact'].fillna("none")

    tempdat['clinvar_score']=tempdat['ClinVar Significance'].apply(
        lambda x: 3 if "athogenic" in x else 2 if 'onflicting' in x else 1 if 'enign' in x else 0)
    tempdat['impact_score']=tempdat['Impact'].apply(
        lambda x: 3 if x=="HIGH" else 2 if x=="MODERATE" else 1 if x=="LOW" else 0)
    maxscores=['impact_score','clinvar_score','REVEL Score','AF PopMax gnomAD EX','Max Polyphen2 HVAR score',
              'Max_Polyphen2 HDIV score','MetaLR_score','Homozygous Count gnomAD EX ']
    aggdict={x:'max' for x in maxscores}
    heydat=tempdat.groupby(['samp','Gene']).agg(aggdict)
    eydat=heydat.reset_index().pivot(index='samp',columns='Gene').reset_index()
    eydat.columns=['_'.join(col) for col in eydat.columns]
    eydat=eydat.rename(columns={'samp_':'samp'})
    newdat=pd.merge(left=ydat,right=eydat,on='samp',how='inner')
    return(newdat)


# In[4]:


## MODEL
from sklearn.ensemble import RandomForestClassifier
def grab_featimp(mymod):
    locdict={}
    locdict['feature']=mymod.get_booster().feature_names
    locdict['importance']=mymod.feature_importances_
    outdat=pd.DataFrame(locdict).sort_values('importance',ascending=False)
    return(outdat)
def modfit(df,curfeats,modtype):
    X=df[curfeats]
    Y=df['y']
    if modtype=='xgb':
        model = xgb.XGBClassifier(max_depth=1,
                                  colsample_bytree=0.001,
                                  n_estimators=500,
                                  seed=np.random.randint(0,100000)
                                 )
        # old model
        #model = xgb.XGBClassifier(max_depth=1,
        #                          colsample_bytree=0.001,
        #                          n_estimators=500,
        #                          seed=np.random.randint(0,100000)
        #                         )

    if modtype=='glm':
        model = LogisticRegression(solver="liblinear")
    if modtype=='rf':
        model=RandomForestClassifier(n_estimators=100,
                                     #criterion='entropy',
                                     warm_start =True,
                                    # max_leaf_nodes =4,
                                    #max_depth =4,
                                     random_state=np.random.randint(0,100000)
                                    )
    if modtype=='svm':
        model=svm.SVC(probability=True)
    model.fit(X,Y)
    return model

def get_preds(testdat,curfeats,mod):
    X=testdat[curfeats]
    xpred_pre=mod.predict_proba(X)
    return([x[1] for x in xpred_pre])
def impute_means(dfpre,feats):
    df=dfpre.copy()
    for x in feats:
        df[x]=df[x].fillna(df[x].mean())
    return(df)
def onehot_encode_cold(dfpre,catcols):
    df=dfpre.copy()
    ohcols=[]
    for x in catcols:
        catvals=list(df[x].unique())
        for i in range(len(catvals)-1):
            curval=catvals[i]
            df[x+'_'+str(curval)]=0
            df.loc[df[x]==curval,x+'_'+str(curval)]=1
        ohcols=ohcols+[x+'_'+str(y) for y in catvals[:(len(catvals)-1)]]
    return(df,ohcols)
def eval_model_new(df,mymod):
    df.loc[:,'yprob']=mod.predict_proba(df)[:,1]
    return(df)

def eval_model(cestdat_pre,incols,mymod):
    cestdat=cestdat_pre.copy()
    cestdat.loc[:,'yprob']=get_preds(cestdat,incols,mymod)
    cestdat.loc[:,'ypred']=1*(cestdat['yprob']>=.3)
    outdat=cestdat
    if 'NCH MSI RESULTS' not in outdat.columns:
        outdat['NCH MSI RESULTS']='dummyvar'
    outdat['NCH MSI RESULTS']=outdat['NCH MSI RESULTS'].fillna('dummyvar')
    outdat['y_nch']=outdat['NCH MSI RESULTS'].apply(lambda x: 1 if "MSI" in x else 0 if "MSS" in x else -1)
    #outdat[['yprob','ypred','y_nch', 'y','Tumor_Type', 'NCH MSI RESULTS','MSI Probability NCH Sequencing Results']]
    return(outdat)
def grab_feats(df):
    markers=[x.replace('0_','') for x in df.columns if x[:2]=='0_']
    outcols=[x for x in df.columns for y in markers if y in x]
    return(list(set(outcols)))


# In[1]:


import matplotlib.pyplot as plt
def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


# In[5]:


def stack_matrices_to_df(df):
    temp=df
    nullinds=temp[temp[0].isnull()].index
    if (len(nullinds) % 4) != 0:
        print('# null inds must be multiple of 4')
        return(None)
    nmark=int(len(nullinds)/4)
    dfdict={}
    for i in range(nmark):
        startind=nullinds[4*i+1]+1
        openstopind=nullinds[4*i+2]-1
        locdat=df.loc[startind:openstopind]
        marker=df.loc[nullinds[4*i]][1]
        locdat=locdat.drop(columns=0)
        newcols=df.drop(columns=0).loc[startind-1].values
        #print(len(newcols))
        locdat.columns=newcols
        transdat=locdat.transpose()
        transdat.columns=[str(x)+'_'+marker for x in range(len(locdat))]
        utdat=transdat.transpose()
        dfdict[i]=utdat
    catdat=pd.concat(dfdict).reset_index().drop(columns='level_0').rename(columns={'level_1':'marker'})
    newdat=catdat.drop(columns='marker').transpose()
    newdat.columns=catdat.marker
    hewdat=newdat.reset_index().rename(columns={'index':'samp'})
    return(hewdat)


# In[ ]:





# In[ ]:





# In[ ]:





