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

	usage = "python run.py [options]\n"
	usage += "to see details, python run.py -h\n"
	parser = argparse.ArgumentParser(description = usage)

	parser.add_argument("-I" , "--input_prefix", type=str, help="input fastq is must be named as <prefix>_<R1 or R2>.fastq.gz")
	parser.add_argument("-T" ,"--target_info", type=str, help="information of targeted region")
	parser.add_argument("-W" , "--work_dir", type=str, default="./", help="Working directory (base dir of output files)")
	parser.add_argument("-b" ,"--bwa_index", type=str , help="<dir>/prefix of bowtie index")
	parser.add_argument("-t" , "--threads", type=int, default=1, help="Number of threads")

	args = parser.parse_args()

	return args


def main():
	args=add_args()
	n_core=args.threads
	WD=args.work_dir
	target_path=args.target_info
	input=args.input_prefix
	align_index=args.bwa_index

	print "-------target_info save"
	target_info=pd.read_csv(target_path, sep="\t", index_col="target", keep_default_na=False)
	targets=list(target_info.index)
	# 1. trim flanking, extract UMI, filter incorrect UMI
	print "-------Trimming"
	trim.trim_sequence(WD, input, target_info)

	# 2. align
	print "-------Alignment"
	align.bwa(WD, input, n_core, align_index)


if __name__=="__main__":
	main()
