#!/usr/bin/env python
import argparse
import gzip
import os
import sys

import align
import base_count
import network_construction
import numpy as np
import pair
import pandas as pd
import pysam
import trim
from Bio import SeqIO


def add_args():

	usage = "python RUN-all.py [options]\n"
	usage += "to see details, python RUN-all.py -h\n"
	parser = argparse.ArgumentParser(description = usage)

	parser.add_argument("-I" , "--input_prefix", type=str, help="input fastq is must be named as <prefix>_<R1 or R2>.fastq.gz")
	parser.add_argument("-O", "--output_prefix", type=str, help="prefix of output files")
	parser.add_argument("-T" ,"--target_info", type=str, help="information of targeted region (tsv format)")
	parser.add_argument("-t" , "--threads", type=int, default=1, help="threads (CPU, default : 1)")
	parser.add_argument("-W" , "--work_dir", type=str, default="./", help="Working directory (base dir)")
	parser.add_argument("-y" ,"--type", type=str, default="c", help="type for basecall [c/o/u/b], default : c")
	
    ## arguments for alignment
	parser.add_argument("-b" ,"--bwa_index", type=str, help="<dir>/prefix of bwa index")
	### arguments for depth-filtering
	parser.add_argument("-C" ,"--cycle", type=int ,default=8, help="Number of cycle [int]")
	parser.add_argument("-u", "--min_UID_pair", type=int, default=2, help="at least depth of UID")
	parser.add_argument("-m" ,"--min_base_ratio", type=float ,default=0.7, help="minimum frequency of major base per CID")
	
	### arguments for indel analysis
	parser.add_argument("-M" ,"--indel_mode", type=bool ,default=False, help="if indel mode : -M")
	parser.add_argument("-v" ,"--indel_vcf", type=str ,default="", help="vcf path for target indel")

	### skip
	parser.add_argument("-1" ,"--skip_preprocess", type=bool ,default=False, help="True if trimming has been proceeded already")
	parser.add_argument("-2" ,"--skip_network_build", type=bool ,default=False, help="True if trimming and network building has been proceeded already")
	args = parser.parse_args()
	return args


def main():
    ### initialization of parameters
	args=add_args()

	input=args.input_prefix
	sample=args.output_prefix
	target_path=args.target_info
	threads=args.threads
	WD=args.work_dir
	type=args.type

	align_index=args.bwa_index

	cycle=args.cycle
	min_nUID=args.min_UID_pair
	min_base_ratio=args.min_base_ratio

	indel_mode=args.indel_mode
	indelvcf=args.indel_vcf
	skip1=args.skip_preprocess
	skip2=args.skip_network_build

    ### preprocessing
	print "-------target_info save"
	target_info=pd.read_csv(target_path, sep="\t", index_col="target", keep_default_na=False)
	targets=list(target_info.index)

	if skip1==False:
		print "-------Trimming"
		# trim.trim_sequence(WD, input)
		print "-------Alignment"
		# align.bwa(WD, input, threads, align_index)

	elif skip1==True:
		print "Skipping trimming step"

    ### build peer-to-peer network
	if type=="c" and skip2!=True:
		print "-------pair_save"
		# pair.pair_save(WD, input, sample, target_info)
		print "-------network_construction"
		# network_construction.network_const(WD, sample, target_info, threads, cycle)

	else:
		print "Skip building peer-to-peer network"


    ### base call
	if indel_mode:
		if not os.path.exists("{}/indel_call".format(WD)):
			os.mkdir("{}/indel_call".format(WD))
		if type=="c":
			print "-------cluster based indel count"
			base_count.cluster_based_indel(WD, input, sample, threads, target_info, min_base_ratio, min_nUID, indelvcf, align_index)

	else:
		if not os.path.exists("{}/base_call".format(WD)):
			os.mkdir("{}/base_call".format(WD))
		# cluster-based
		if type=="c":
			print "-------cluster based base count"
			base_count.cluster_based_call(WD, input, sample, threads, target_info, min_base_ratio, min_nUID)

		# UID-pair based base call
		elif type=="u":
			print "-------UID-pair based deduplication"
			base_count.uid_based_call(WD,input, sample, target_info, threads, min_base_ratio)

		# basic base count
		elif type=="b":
			print "-------Base call with conventional base counting"
			base_count.basic_base_count(WD,input, sample, target_info)

		else:
			print "type option Error"

if __name__=="__main__":
	main()
