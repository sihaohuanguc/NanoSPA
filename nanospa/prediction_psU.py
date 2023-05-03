#!/usr/bin/env python
# -*- coding: utf-8

import os,numpy,pickle
import pandas as pd
from pandas.core.frame import DataFrame
from sklearn.ensemble import ExtraTreesClassifier

__author__ = "Sihao Huang"
__copyright__ = ""
__credits__ = []
__license__ = "GPL 2.0"
__version__ = "1.0"
__maintainer__ = "Sihao Huang"
__email__ = "sihaohuang1024@gmail.com"
__status__ = "Development"

def pred(in_name,out_name):
    working_path=os.getcwd()
    if not in_name:
        in_folder=working_path+"/alignment"
        in_file=in_folder+"/features.csv"
        if not out_name:
            out_prediction_file=working_path+"/alignment/prediction_psU.csv"
        else:
            out_prediction_file=out_name
    else:
        if "/" in in_name:
            in_folder="/".join(in_name.split("/")[:-1])
        else:
            in_folder=working_path
        in_file=in_name
        if not out_name:
            out_prediction_file=in_folder+"/prediction_psU.csv"
        else:
            out_prediction_file=out_name
    
    # print(sys.path[0])
    base_path=os.path.abspath(__file__)
    folder=os.path.dirname(base_path)
    in_model=folder+"/data/model/psU_model_3_1.pkl"

    numpy.set_printoptions(suppress=True)

    with open(in_model,"rb") as f:
        ET1=pickle.load(f)

    df=pd.read_csv(in_file,header=None,names=["transcript_ID","position","base_type","coverage","ins","ins_len","del","del_len","fuzzy","mis","misA","misC","misG","misT","base_qual_mean","base_qual_STD","base_qual_count_0"])
    df=df[df["base_type"]=="T"]
    cols=["transcript_ID","position","base_type","coverage"]
    site_info=df[cols]
    df.drop(["transcript_ID","position","base_type","coverage","misT"],axis=1,inplace=True)

    ET1_predict_proba=ET1.predict_proba(df)
    pred_df=pd.DataFrame(ET1_predict_proba,columns=["unmodi","modi"])
    site_info=site_info.reset_index(drop=True).join(pred_df["modi"])
    site_info.to_csv(out_prediction_file,index=False,header=None)



