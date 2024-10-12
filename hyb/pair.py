#!/usr/bin/env python
import pysam
import pandas as pd
import numpy as np
import argparse
from utils import set_pos
import os

def pair_save(WD, input, sample, target_info):
	if not os.path.exists("{}/pair".format(WD)):
		os.mkdir("{}/pair".format(WD))

	targets=target_info.index
	chrs=target_info["chromosome"]
	poslist=target_info["pos"]
	loci=zip(targets, chrs, poslist)

	bam_data={}
	for seg in ["R1", "R2"]:
		bamfile=pysam.AlignmentFile("{}/align/{}_{}.sorted.bam".format(WD, input, seg), "rb")
		for target, chr, poslist in loci:
			poslist=set_pos(poslist)
			minpos=min(poslist)-10
			maxpos=max(poslist)+10
			for read in bamfile.fetch(chr, minpos, maxpos ):
				id=read.qname
				read_id=id.split("_")[0]
				UID=id.split("_")[-1]
				chr=read.reference_id
				cigar=read.cigarstring
				mapq=read.mapping_quality
				if mapq<55 or "S" in cigar: continue
				if min(read.get_reference_positions())-1>min(poslist): continue
				stt=min(read.get_reference_positions())+1
				end=max(read.get_reference_positions())+1
				try:
					chk=bam_data[read_id+"_"+target]
				except KeyError:
					bam_data[read_id+"_"+target]=[]
				bam_data[read_id+"_"+target].append([UID, target, stt, end])
	
	df=[]
	for k in bam_data.keys():
		if len(bam_data[k])!=2: continue
		tmp=bam_data[k]
		id=k.split("_")[0]
		fragmentpos=tmp[0][2:]
		fragmentpos.extend(tmp[1][2:])
		df.append([id, tmp[0][0], tmp[1][0],  tmp[0][1], min(fragmentpos), max(fragmentpos), str(min(fragmentpos))+"-"+str(max(fragmentpos))])

	df=pd.DataFrame(df, columns=["id", "R1_UID", "R2_UID",  "target", "stt", "end", "stt-end"])
	df["frag-R1U"]=df["target"]+":"+df["stt-end"]+":"+df["R1_UID"]
	df["frag-R2U"]=df["target"]+":"+df["stt-end"]+":"+df["R2_UID"]
	df["pair"]=df["frag-R1U"]+"_"+df["frag-R2U"]
	df.to_csv("{}/pair/{}_barcode_pair.txt".format(WD, sample), sep="\t", index=False)

