import pysam
import pandas as pd
import numpy as np
import argparse
import os
def get_pos_list(poslist):
	if ";" in poslist:
		poslist=poslist.split(";")
	elif "-" in poslist:
		poslist=range(int(poslist.split("-")[0]), int(poslist.split("-")[1])+1)
	else:
		poslist=[poslist]
	return [int(pos) for pos in poslist]

def pair_save(WD, input, sample, target_info, cutoff=1):
	if not os.path.exists("{}/pair".format(WD)):
		os.mkdir("{}/pair".format(WD))

	targets=target_info.index
	chrs=target_info["chromosome"]
	poslist=target_info["pos"]
	loci=zip(targets, chrs, poslist)

	## read UID from bam
	bam_data={}
	for seg in ["R1", "R2"]:
		bamfile=pysam.AlignmentFile("{}/align/{}_{}.sorted.bam".format(WD, input, seg), "rb")
		for target, chr, poslist in loci:
			poslist=get_pos_list(poslist)
			for read in bamfile.fetch(chr, min(poslist), max(poslist) ):
				id=read.qname
				read_id=id.split("_")[0]
				UID=id.split("_")[-1]
				target_in_bam=id.split("_")[-2]
				if target_in_bam!=target: continue
				if min(read.get_reference_positions())-1>min(poslist): continue
				chr=read.reference_id
				cigar=read.cigarstring
				mapq=read.mapping_quality
				if mapq<55 or "S" in cigar: continue
				try:
					chk=bam_data[read_id+"_"+target]
				except KeyError:
					bam_data[read_id+"_"+target]=[]
				bam_data[read_id+"_"+target].append([UID, target])
	
	df=[]
	for k in bam_data.keys():
		if len(bam_data[k])!=2:
			continue
		id=k.split("_")[0]
		tmp=bam_data[k]
		df.append([id, tmp[0][0], tmp[1][0], tmp[0][0]+"_"+tmp[1][0], tmp[0][1]])
	df=pd.DataFrame(df, columns=["id", "R1_UID", "R2_UID", "pair", "target"])
	## filter by pair redundancy
	uniq=np.unique(df["pair"], return_counts=True)
	count_dic=dict(zip(uniq[0], uniq[1]))
	df["count"]=df["pair"].apply(lambda x: count_dic[x])
	filtered=df.loc[df["count"]>=cutoff,["id", "R1_UID", "R2_UID", "pair", "target"]]
	filtered.to_csv("{}/pair/{}_barcode_pair.txt".format(WD, sample), sep="\t", index=False)

