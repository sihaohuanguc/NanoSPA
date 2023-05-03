#!/usr/bin/env python
# -*- coding: utf-8

import os,re

__author__ = "Sihao Huang"
__copyright__ = ""
__credits__ = []
__license__ = "GPL 2.0"
__version__ = "1.0"
__maintainer__ = "Sihao Huang"
__email__ = "sihaohuang1024@gmail.com"
__status__ = "Development"

def re_in(in_name):

    working_path=os.getcwd()
    if not in_name:
        in_folder=working_path+"/alignment"
    else:
        in_folder=in_name

    for item in os.listdir(in_folder):
        sub_folder=in_folder+"/"+item
        if os.path.isdir(sub_folder) and "collect_pile.txt" in os.listdir(sub_folder):
            in_file=sub_folder+"/collect_pile.txt"
            out_file=sub_folder+"/collect_pile_no_intron.txt"
            with open(out_file,"w+") as out_f:
                for line in open(in_file,"r"):
                    site=line.strip("\n").split("\t")
                    out_line=site[:3]
                    num_intron_base=site[4].count(">")+site[4].count("<")
                    out_line.append(str(int(site[3])-num_intron_base))
                    
                    pattern=re.compile("\\^.")
                    alignment=re.sub(pattern,"",site[4]) 
                    qual=list(site[5])
                    site_discard_switch=0 
                    p=0
                    for i in range(int(site[3])):  
                        if alignment[p]=="." or alignment[p]==",":
                            if len(alignment)-p>1 and alignment[p+1]=="+":
                                pattern2=re.compile("\\+[0-9]+")
                                sub_str=re.search(pattern2,alignment[p:]).group()
                                num_base=sub_str.strip("+")
                                p+=(int(num_base)+len(num_base)+2)
                            elif len(alignment)-p>1 and alignment[p+1]=="-":
                                pattern2=re.compile("\\-[0-9]+")
                                sub_str=re.search(pattern2,alignment[p:]).group()
                                num_base=sub_str.strip("-")
                                p+=(int(num_base)+len(num_base)+2)
                            elif len(alignment)-p>1 and alignment[p+1]=="$":
                                p+=2                     
                            else:
                                p+=1
                        elif alignment[p] in "*ACGTacgt":
                            p+=1
                        elif alignment[p]==">" or alignment[p]=="<":
                            qual[i]=""
                            if len(alignment)-p>1 and alignment[p+1]=="+":
                                pattern2=re.compile("\\+[0-9]+")
                                sub_str=re.search(pattern2,alignment[p:]).group()
                                num_base=sub_str.strip("+")
                                p+=(int(num_base)+len(num_base)+2)
                            elif len(alignment)-p>1 and alignment[p+1]=="-":
                                pattern2=re.compile("\\-[0-9]+")
                                sub_str=re.search(pattern2,alignment[p:]).group()
                                num_base=sub_str.strip("-")
                                p+=(int(num_base)+len(num_base)+2)
                            else:
                                p+=1
                        else:
                            site_discard_switch=1 
                            break
                    if site_discard_switch==0:        
                        pattern_intron_and_indel=re.compile("[<>][-+][0-9]+[ATCGNatcgn]+") 
                        alignment_state=site[4]
                        results=re.search(pattern_intron_and_indel,alignment_state)
                        while not results==None:
                            pattern_num=re.compile("[-+][0-9]+")
                            sub_str=re.search(pattern_num,results.group()).group()
                            num_base=sub_str.strip("-").strip("+")
                            delete_part=alignment_state[results.start():(results.start()+int(num_base)+len(num_base)+2)]
                            alignment_state=alignment_state.replace(delete_part,"",1) 
                            results=re.search(pattern_intron_and_indel,alignment_state)
                        alignment_state=alignment_state.replace(">","").replace("<","")
                        out_line.append(alignment_state)
                        out_line.append("".join(qual))
                        out_f.write("\t".join(out_line)+"\n")    
                        
