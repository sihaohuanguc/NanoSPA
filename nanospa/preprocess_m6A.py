#!/usr/bin/env python
# -*- coding: utf-8

import os,csv
import pandas as pd

__author__ = "Sihao Huang"
__copyright__ = ""
__credits__ = []
__license__ = "GPL 2.0"
__version__ = "1.0"
__maintainer__ = "Sihao Huang"
__email__ = "sihaohuang1024@gmail.com"
__status__ = "Development"

def pre_m6A(in_name,out_name):
    working_path=os.getcwd()
    if not in_name:
        in_folder=working_path+"/alignment"
        in_file=in_folder+"/features.csv"
        if not out_name:
            out_file=working_path+"/alignment/features_GGACU.csv"
        else:
            out_file=out_name
    else:
        if "/" in in_name:
            in_folder="/".join(in_name.split("/")[:-1])
        else:
            in_folder=working_path
        in_file=in_name
        if not out_name:
            out_file=in_folder+"/features_GGACU.csv"
        else:
            out_file=out_name
    temp_file=in_folder+"/features_GGACU_temp.csv" 

    motif="GGACT" # use T instead of U !!!
    modification_index=2 # this is the position-1 of the modification site A in the motif GGACT

    queue_base=[]  # for bases
    queue_position=[]  # for the positions
    queue_info=[]  # for all the info
    have_write_something=False
    with open(temp_file,"w+") as out_f:
        writer=csv.writer(out_f)
        with open(in_file,"r") as inf:
            reader=csv.reader(inf)
            for line in reader:        
                queue_base.append(line[2])
                queue_position.append(int(line[1]))
                queue_info.append(line)
                if len(queue_base)==len(motif): # now we have enough bases
                    if "".join(queue_base)==motif:  # it is a motif needed
                        if [x-queue_position[0] for x in queue_position]==[i for i in range(len(motif))]: # the positions are neighbors, and the difference will be [0,1,2,3,4] for 5 mers
                            # print(queue_info)
                            out_line=queue_info[modification_index][:4] # the base info of the modification position
                            for info in queue_info:
                                out_line.extend(info[4:])
                            writer.writerow(out_line)
                            have_write_something=True
                            
                    queue_base.pop(0)  # pop out the first site
                    queue_position.pop(0)
                    queue_info.pop(0)

    # there are too many columns, let's generate the name of the columns
    cols=["gene_ID","position","base_type","coverage"]
    base_col=["ins","ins_len","del","del_len","del_site","mis","mis_A","mis_C","mis_G","mis_T","base_qual_mean","base_qual_STD","base_qual_count_0"]
    for i in range(1,6):
        for item in base_col:
            cols.append(f"{item}_{i}")
    # print(cols)
    if have_write_something:
        df=pd.read_csv(temp_file,header=None)
        df.columns=cols

        feature_order=["gene_ID","position","base_type","coverage",'base_qual_mean_3', 'base_qual_mean_2', 'base_qual_mean_4', 'del_site_3', 'del_2', 'base_qual_mean_1', 'del_len_2', 'mis_3', 'mis_4', 'base_qual_mean_5', 'base_qual_STD_4', 'mis_T_4', 'mis_2', 'mis_5', 'mis_G_3', 'base_qual_STD_3', 'base_qual_STD_5', 'mis_A_2', 'base_qual_count_0_3', 'ins_2', 'mis_1', 'base_qual_count_0_2', 'base_qual_count_0_4', 'mis_T_3', 'del_site_2', 'mis_A_5', 'base_qual_STD_2', 'mis_C_2', 'mis_C_3', 'ins_len_4', 'base_qual_count_0_5', 'ins_4', 'mis_A_1', 'ins_len_2', 'del_1', 'base_qual_count_0_1', 'mis_A_4', 'del_site_4', 'del_len_1', 'base_qual_STD_1', 'ins_3', 'del_site_5', 'mis_C_5', 'mis_C_1', 'del_len_4', 'del_4', 'del_site_1', 'ins_len_3', 'ins_1', 'del_len_3', 'del_3', 'ins_len_1', 'mis_T_2', 'mis_G_4', 'mis_T_1', 'del_5', 'del_len_5', 'ins_5', 'ins_len_5', 'mis_G_5']
        df=df.reindex(columns=feature_order)
        df.to_csv(out_file,header=True,index=False)
    else:
        print("No A site available.")

    os.system("rm "+temp_file)
