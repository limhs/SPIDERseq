#!/usr/bin/env python
import argparse
import gzip, os, sys

import pair_oligo, trim_oligo, network_construction, specificity
import numpy as np
import pandas as pd
import pysam
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
	
	### arguments for depth-filtering
	parser.add_argument("-C" ,"--cycle", type=int ,default=8, help="Number of cycle [int]")
	parser.add_argument("-u", "--min_UID_pair", type=int, default=2, help="at least depth of UID")
	parser.add_argument("-m" ,"--min_base_ratio", type=float ,default=0.7, help="minimum frequency of major base per CID")

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

	cycle=args.cycle
	min_nUID=args.min_UID_pair
	min_base_ratio=args.min_base_ratio

	skip1=args.skip_preprocess
	skip2=args.skip_network_build

    ### preprocessing
	print "-------target_info save"
	target_info=pd.read_csv(target_path, sep="\t", index_col="target", keep_default_na=False)
	targets=list(target_info.index)

	# if skip1==True or skip2==True:
	# 	print "Skipping trimming step"

	# else:
	# 	print "-------Trimming"
	# 	trim_oligo.trim_sequence(WD, input, target_info)

    # ### build peer-to-peer network
	# if skip2==True:
	# 	print "Skip building peer-to-peer network"

	# else:
	# 	print "-------pair_save"
	# 	pair_oligo.pair_save(WD,input, sample, target_info)
	# 	print "-------network_construction"
	# 	network_construction.network_const(WD, sample, target_info, threads, cycle)

	### caculate specificity
	print "-------specificity"
	specificity.specificity(WD,input, sample, threads, target_info)

if __name__=="__main__":
	main()
