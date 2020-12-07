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

	usage = "python RUN-BASECALL.py [options]\n"
	usage += "to see details, python RUN-BASECALL.py -h\n"
	parser = argparse.ArgumentParser(description = usage)

	parser.add_argument("-o", "--output_prefix", type=str, help="prefix of output files")
	parser.add_argument("-I" , "--input_prefix", type=str, help="input fastq is must be named as <prefix>_<R1 or R2>.fastq.gz")
	parser.add_argument("-T" ,"--target_info", type=str, help="information of targeted region")
	parser.add_argument("-W" , "--work_dir", type=str, default="./", help="Working directory (base dir of output files)")
	parser.add_argument("-y" ,"--type", type=str, default="c", help="type for basecall [c/o/u/b]")
	parser.add_argument("-t" , "--threads", type=int, default=1, help="Number of threads")
	
	### arguments for depth-filtering
	parser.add_argument("-u", "--min_umi_pair", type=int, default=1, help="minimum depth of umi-pair")
	parser.add_argument("-m" ,"--min_base_ratio", type=float ,default=0.7, help="cutoff frequency for major base")
	parser.add_argument("-n", "--min_nread", type=int, default=3, help="minimum nRead per cluster")
	parser.add_argument("-C" ,"--cycle", type=int ,default=8, help="Number of cycle [int]")
	
	### arguments for indel analysis
	parser.add_argument("-b" ,"--bwa_index", type=str , help="<dir>/prefix of bowtie index")
	parser.add_argument("-M" ,"--indel_mode", type=bool ,default=False, help="if indel mode : -M")
	parser.add_argument("-v" ,"--indel_vcf", type=str ,default="", help="vcf path for target indel")
	args = parser.parse_args()

	return args


def main():
	args=add_args()
	n_core=args.threads
	sample=args.output_prefix
	WD=args.work_dir
	target_path=args.target_info
	min_pair_depth=args.min_umi_pair
	align_index=args.bwa_index
	min_nread=args.min_nread
	type=args.type
	min_base_ratio=args.min_base_ratio
	input=args.input_prefix
	cycle=args.cycle
	indel_mode=args.indel_mode
	indelvcf=args.indel_vcf

	print sample
	print "type is : "+type
	print "-------target_info save"
	target_info=pd.read_csv(target_path, sep="\t", index_col="target", keep_default_na=False)
	targets=list(target_info.index)


	if indel_mode:
		if not os.path.exists("{}/indel_call".format(WD)):
			os.mkdir("{}/indel_call".format(WD))
		if type=="c":
			print "-------cluster based indel count"
			base_count.cluster_based_indel(WD, input,sample,n_core, target_info, min_base_ratio, min_nread,  indelvcf, align_index)

	else:
		if not os.path.exists("{}/base_call".format(WD)):
			os.mkdir("{}/base_call".format(WD))
		# 3-1 cluster-based
		if type=="c":
			print "-------cluster based base count"
			base_count.cluster_based_call(WD, input,sample,n_core, target_info, min_base_ratio, min_nread)
		# 3-2 origin-based
		elif type=="o":
			print "-------origin based base count"
			base_count.origin_based_call(WD, input,sample,n_core, target_info, min_base_ratio, min_nread)

		# 3-3. UMI-pair based base call
		elif type=="u":
			print "-------umi-pair based deduplication"
			base_count.umi_based_call(WD,input, sample, target_info, n_core, min_base_ratio)

		# 3-4. basic base count
		elif type=="b":
			print "-------Base call with conventional base counting"
			base_count.basic_base_count(WD,input, sample, target_info)

		else:
			print "type option Error"

if __name__=="__main__":
	main()
