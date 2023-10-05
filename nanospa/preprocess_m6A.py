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
            out_folder=working_path+"/alignment/features"
        else:
            out_folder=out_name.rstrip("/")
    else:
        if "/" in in_name:
            in_folder="/".join(in_name.split("/")[:-1])  # to remove file name
        else:
            in_folder=working_path
        in_file=in_name
        if not out_name:
            out_folder=in_folder+"/features"
        else:
            out_folder=out_name.rstrip("/")
    temp_folder=in_folder+"/features_temp"

    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    motif_list=["GGACT","AGACT","GGACA","TGACT","GGACC","AGACA","AGACC","GGACG"]
    
    for motif in motif_list:
        temp_file=temp_folder+f"/features_temp_{motif}.csv"
        out_file=out_folder+f"/features_{motif}.csv"
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

            df.to_csv(out_file,header=True,index=False)
        else:
            print(f"No A site available in {motif}.")

    os.system("rm -r "+temp_folder)
