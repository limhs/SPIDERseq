import os
from multiprocessing import Process, Queue
import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO

def set_pos(poslist):
	if ";" in poslist:
		poslist=poslist.split(";")
	elif "-" in poslist:
		poslist=range(int(poslist.split("-")[0]),int(poslist.split("-")[1])+1)
	elif ".." in poslist:
		poslist=range(int(poslist.split("..")[0]),int(poslist.split("..")[1])+1)
	else:
		poslist=[poslist]
	poslist=[int(pos) for pos in poslist]
	return poslist
	
def get_loci(target_info):
	targets=target_info.index
	chrs=target_info["chromosome"]
	poslist=target_info["pos"]
	loci=zip(targets, chrs, poslist)
	return loci

def run_multi(thread_function, n_core, df, WD, input,sample, target_info, min_base_ratio, other_args=() ):
	Ps=[]
	unit=len(df.index)/n_core
	

	for n in range(1,n_core+1):
	    stt=(n-1)*unit
	    end=n*unit-1
	    if n==n_core:
	        end=(n+1)*unit-1

	    procs=Process(target=thread_function, args=(WD, input, sample,n, target_info, df.loc[df.index[stt:end]], min_base_ratio)+other_args )
	    Ps.append(procs)
	    procs.start()

	for p in Ps:
	    p.join()

	return "complete"

def output_merge(target,chr, poslist, dfdic, out_df, min_base_ratio):
	dff_tmp=pd.DataFrame()
	for pos in poslist:
		locus_id=target+"_"+str(pos)
		df=dfdic[locus_id]

		df["total"]=0
		for b in "ATCGN":
			df["total"]+=df[b]
		df=df[df["total"]>0]
		
		df["max_base"]=df[list("ATCGN")].idxmax(axis=1)
		df["max_base_count"]=df[list("ATCGN")].max(axis=1)
		df["max_base_ratio"]=df["max_base_count"].astype("float")/df["total"]
		df["position"]=chr+":"+str(pos)

		df_filtered=df.loc[df["max_base_ratio"]>min_base_ratio, :]
		dff_tmp=dff_tmp.append(df_filtered)
		max_base_list=list(df_filtered["max_base"])
		tmp=[locus_id]
		for b in "ATCGN":
			tmp.append(max_base_list.count(b))
		out_df.append(tmp)
	return (out_df, dff_tmp)

def get_vcf_info(vcfpath, refdic):
	vcf=open(vcfpath)
	vcf_info=[]
	for line in vcf:
		if line.startswith("#"): continue
		tmp=line.split('\t')
		chr=tmp[0].strip()
		pos=tmp[1].strip()
		ref=tmp[3].strip()
		alt=tmp[4].strip()
		mut="{}:{}/{}>{}".format(chr,pos,ref,alt)
		left_flk=refdic[chr][int(pos)-10:int(pos)]
		right_flk=refdic[chr][int(pos)+len(ref):int(pos)+len(ref)+10]
		ref_query=left_flk+ref+right_flk
		alt_query=left_flk+alt+right_flk
		form=[chr,int(pos),ref,alt, mut, ref_query, alt_query]
		vcf_info.append(form)
	return pd.DataFrame(vcf_info, columns=["chr", "pos", "ref", "alt", "mut", "ref_query", "alt_query"])

def get_refseq(refpath):
	reffasta=SeqIO.parse(open(refpath), "fasta")
	ref={}
	for chrom in reffasta:
		ref[chrom.id]="S"+str(chrom.seq)
	return ref
