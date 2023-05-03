#!/usr/bin/env python
# -*- coding: utf-8

import tensorflow as tf
import os
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

__author__ = "Sihao Huang"
__copyright__ = ""
__credits__ = []
__license__ = "GPL 2.0"
__version__ = "1.0"
__maintainer__ = "Sihao Huang"
__email__ = "sihaohuang1024@gmail.com"
__status__ = "Development"

def pred_m6A(in_name,out_name):
    working_path=os.getcwd()
    if not in_name:
        in_folder=working_path+"/alignment"
        in_file=in_folder+"/features_GGACU.csv"
        if not out_name:
            out_prediction_file=working_path+"/alignment/prediction_m6A.csv"
        else:
            out_prediction_file=out_name
    else:
        if "/" in in_name:
            in_folder="/".join(in_name.split("/")[:-1])
        else:
            in_folder=working_path
        in_file=in_name
        if not out_name:
            out_prediction_file=in_folder+"/prediction_m6A.csv"
        else:
            out_prediction_file=out_name

    df=pd.read_csv(in_file) # this table has header line
    cols=["gene_ID","position","base_type","coverage"]
    site_info=df[cols]
    df.drop(["gene_ID","position","base_type","coverage","mis_1","mis_2","mis_3","mis_4","mis_5"],axis=1,inplace=True)  # remove all overall mismatch features

    scaler=MinMaxScaler()
    scaler.fit(df)
    df_standarized=pd.DataFrame(scaler.transform(df),columns=df.columns)

    model=tf.keras.models.Sequential([
        tf.keras.layers.Dense(128,activation="relu",input_shape=(55,)),
        tf.keras.layers.Dropout(0.2),
        tf.keras.layers.Dense(64,activation="relu"),
        tf.keras.layers.Dropout(0.1),
        tf.keras.layers.Dense(2,activation="softmax"),   # this is the output layer
    ])

    base_path=os.path.abspath(__file__)
    folder=os.path.dirname(base_path)
    model_path=folder+"/data/model/m6A_model_1_0_GGACU/cp.ckpt"
    if os.path.exists(model_path+".index"):
        print("Load the model.")
        model.load_weights(model_path).expect_partial()
    else:
        print("Model for m6A prediction not exist.")
        exit(100)

    prediction=model.predict(df_standarized)
    pred_df=pd.DataFrame(prediction,columns=["unmodi","modi"])
    site_info=site_info.reset_index(drop=True).join(pred_df["modi"]) # you need to reset the index of site_info df
    site_info.to_csv(out_prediction_file,index=False,header=None)









