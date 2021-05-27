import pickle
import logging

import xgboost as xgb

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
            logging.info("Training model "+str(j+1))
        traininds_loc=df.sample(frac=0.8).index
        mymod=modfit(df.loc[traininds_loc],curfeats)
        pkl_connect = open(modfiles[j], 'wb')
        pickle.dump(mymod, pkl_connect)
        pkl_connect.close()
    return(modfiles)
