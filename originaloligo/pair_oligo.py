import pysam
import pandas as pd
import numpy as np
import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
import os

def pair_save(WD, input, sample, target_info, cutoff=1):
	if not os.path.exists("{}/pair".format(WD)):
			os.mkdir("{}/pair".format(WD))

	targets=target_info.index

	fastq_data={}
	for seg in ["R1", "R2"]:
		seqr=SeqIO.parse(open("{}/trimmed/{}_{}.trimmed.fastq".format(WD, input, seg)), "fastq")
		for read in seqr:
			read_id=read.id.split("_")[0]
			seq=str(read.seq)
			uid=read.id.split("_")[-1]
			target=read.id.split("_")[-2]
			chk_seq=target_info.loc[target, "check_sequence"]
			N1reg=target_info.loc[target, "N1reg"]
			N2reg=target_info.loc[target, "N2reg"]
			if seg=="R2":
				chk_seq=str(Seq(chk_seq).reverse_complement())
				N1reg=str(Seq(N1reg).reverse_complement())
				N2reg=str(Seq(N2reg).reverse_complement())

			N1search=re.findall(N1reg, seq)
			N2search=re.findall(N2reg, seq)
			if len(N1search)!=1 or len(N2search)!=1 : continue
			if chk_seq not in seq: continue
			try:
				chk=fastq_data[read_id]
			except KeyError:
				fastq_data[read_id]=[]
			fastq_data[read_id].append([uid, target])

	df=[]
	for id in fastq_data.keys():
		if len(fastq_data[id])!=2:
			continue
		tmp=fastq_data[id]
		df.append([id, tmp[0][0], tmp[1][0], tmp[0][0]+"_"+tmp[1][0], tmp[0][1]])

	df=pd.DataFrame(df, columns=["id", "R1_UID", "R2_UID", "pair", "target"])

	## filter by pair redundancy	
	uniq=np.unique(df["pair"], return_counts=True)
	count_dic=dict(zip(uniq[0], uniq[1]))
	df["count"]=df["pair"].apply(lambda x: count_dic[x])
	filtered=df.loc[df["count"]>=cutoff,df.columns.isin(["id", "R1_UID", "R2_UID", "pair", "target"])]
	filtered.to_csv("{}/pair/{}_barcode_pair.txt".format(WD, sample), sep="\t", index=False)
