#!/usr/bin/env python
# -*- coding: utf-8

import os

__author__ = "Sihao Huang"
__copyright__ = ""
__credits__ = []
__license__ = "GPL 2.0"
__version__ = "1.0"
__maintainer__ = "Sihao Huang"
__email__ = "sihaohuang1024@gmail.com"
__status__ = "Development"

def align(in_folder,ref_fa,out_name):

    working_path=os.getcwd()
    temp_folder=working_path+"/temp"
    if not out_name:
        out_folder=working_path+"/alignment"
    else:
        out_folder=out_name

    if not os.path.exists(temp_folder):
        os.mkdir(temp_folder)
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
        
    in_temp_fastq=temp_folder+"/all_pass_fastq.fastq"
    with open(in_temp_fastq,"w+") as out_f:
        for item in os.listdir(in_folder):
            in_file=in_folder+"/"+item
            if os.path.isfile(in_file) and not item[0]=="." and item.split(".")[-1]=="fastq":
                with open(in_file,"r") as f:
                    switch=0
                    for line in f:
                        if switch==0 and line[0]=="@":
                            switch=1
                        elif switch==1 and line[0] in ["A","C","G","U"]:
                            switch=0
                            line=line.replace("U","T")
                        out_f.write(line)

    def rev_comp(ref):
        rev_comp=[]
        rev_comp_dict={"A":"T","a":"T","C":"G","c":"G","G":"C","g":"C","T":"A","t":"A","N":"N","n":"N"}
        for item in ref:
            if item in rev_comp_dict.keys():
                rev_comp.append(rev_comp_dict[item])
            else:
                print(item+": Wrong input!")
        rev_comp.reverse()
        return "".join(rev_comp)

    def rev_score(str1):
        list1=list(str1)
        list1.reverse()
        return "".join(list1)

    sam_file=temp_folder+"/all_pass_fastq.sam"
    os.system("minimap2 -ax splice -uf -k14 "+ref_fa+" "+in_temp_fastq+" > "+sam_file)

    out_sub_folder=out_folder+"/plus_strand"
    out_sub_folder1=out_folder+"/minus_strand"
    if not os.path.exists(out_sub_folder):
        os.mkdir(out_sub_folder)
    if not os.path.exists(out_sub_folder1):
        os.mkdir(out_sub_folder1)
    out_file=out_sub_folder+"/collect.fastq"
    out_file1=out_sub_folder1+"/collect.fastq"
    with open(out_file,"w+") as out_f:
        with open(out_file1,"w+") as out_f1:
            with open(sam_file,"r") as in_sam:
                for line in in_sam:
                    if not line[0]=="@" and len(line)>1:
                        elements=line.split("\t")
                        aligned_gene=elements[2]
                        aligned_strand=elements[1]
                        if not aligned_gene=="*":
                            if aligned_strand=="0":
                                out_line=["@"+elements[0],elements[9],"+",elements[10]]
                                for k in out_line:
                                    out_f.write(k+"\n")
                            elif aligned_strand=="16":
                                out_line=["@"+elements[0],rev_comp(elements[9]),"+",rev_score(elements[10])]
                                for k in out_line:
                                    out_f1.write(k+"\n")
                                    
    os.system("rm -r "+temp_folder)
    if os.path.getsize(out_file)==0:
        os.system("rm -r "+out_sub_folder)
    if os.path.getsize(out_file1)==0:
        os.system("rm -r "+out_sub_folder1) 

    with open(ref_fa,"r") as f:
        all_ref=[]
        one_strand=[]
        all_text=[]
        for line in f:
            if line[0]==">":
                all_text="".join(all_text)
                one_strand.append(all_text)
                all_ref.append(one_strand)
                one_strand=[]
                one_strand.append(line.strip("\n").strip(">").split(" ")[0]) 
                all_text=[]
            else:
                all_text.append(line.strip("\n"))
        all_text="".join(all_text)
        one_strand.append(all_text)
        all_ref.append(one_strand)
        all_ref=all_ref[1:]
        
        if "plus_strand" in os.listdir(out_folder):
            sub_folder=out_folder+"/plus_strand"
            out_file=sub_folder+"/reference.fa"
            in_col=sub_folder+"/collect"
            with open(out_file,"w+") as out_f:
                for i in range(len(all_ref)):
                    out_f.write(">"+all_ref[i][0]+"_F\n")
                    out_f.write(all_ref[i][1]+"\n")
            os.system("minimap2 -ax splice -uf -k14 "+out_file+" "+in_col+".fastq > "+in_col+".sam")
            os.system("samtools view -bS "+in_col+".sam > "+in_col+".bam")
            os.system("samtools sort "+in_col+".bam -o "+in_col+".sorted.bam")
            os.system("samtools mpileup --max-depth 1000000 -Q 0 -f "+out_file+" "+in_col+".sorted.bam > "+in_col+"_pile.txt")
        if "minus_strand" in os.listdir(out_folder):
            sub_folder=out_folder+"/minus_strand"
            out_file=sub_folder+"/reference.fa"
            in_col=sub_folder+"/collect"
            with open(out_file,"w+") as out_f:
                for i in range(len(all_ref)):
                    out_f.write(">"+all_ref[i][0]+"_R\n")
                    out_f.write(rev_comp(all_ref[i][1])+"\n")
            os.system("minimap2 -ax splice -uf -k14 "+out_file+" "+in_col+".fastq > "+in_col+".sam")
            os.system("samtools view -bS "+in_col+".sam > "+in_col+".bam")
            os.system("samtools sort "+in_col+".bam -o "+in_col+".sorted.bam")
            os.system("samtools mpileup --max-depth 1000000 -Q 0 -f "+out_file+" "+in_col+".sorted.bam > "+in_col+"_pile.txt")

