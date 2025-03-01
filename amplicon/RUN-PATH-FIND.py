from Bio import SeqIO
import gzip, pysam, sys, os
import pandas as pd
import numpy as np
import argparse
import path_find


def add_args():

	usage = "python run.py [options]\n"
	usage += "to see details, python run.py -h\n"
	parser = argparse.ArgumentParser(description = usage)
	parser.add_argument("-o", "--output_prefix", type=str, help="prefix of output files")
	parser.add_argument("-I" , "--input_prefix", type=str, help="input fastq is must be named as <prefix>_<R1 or R2>.fastq.gz")
	parser.add_argument("-T" ,"--target_info", type=str, help="information of targeted region")
	parser.add_argument("-W" , "--work_dir", type=str, default="./", help="Working directory (base dir of output files)")
	parser.add_argument("-y" ,"--type", type=str, default="c", help="type for basecall [c/o/u/b]")
	parser.add_argument("-t" , "--threads", type=int, default=1, help="Number of threads")
	parser.add_argument("-u", "--min_umi_pair", type=int, default=1, help="minimum depth of umi-pair")
	parser.add_argument("-m" ,"--min_base_ratio", type=float ,default=0.7, help="cutoff frequency for major base")
	parser.add_argument("-n", "--min_nread", type=int, default=3, help="minimum nRead per cluster")
	parser.add_argument("-g" ,"--min_cluster_size", type=int ,default=3, help="lower boundary of cluster size")
	parser.add_argument("-C" ,"--cycle", type=int ,default=8, help="Number of cycle [int]")
	args = parser.parse_args()

	return args


def main():
	args=add_args()
	n_core=args.threads
	sample=args.output_prefix
	WD=args.work_dir
	target_path=args.target_info
	cutoff=args.min_umi_pair
	type=args.type
	min_base_ratio=args.min_base_ratio
	min_nread=args.min_nread
	input=args.input_prefix
	cycle=args.cycle
	min_cluster_size=args.min_cluster_size
	print(sample)
	print("type is : "+type)
	print("-------load target_info")
	target_info=pd.read_csv(target_path, sep="\t", index_col="target", keep_default_na=False)
	targets=list(target_info.index)



	# 3-1 cluster-based
	if type=="c":
		print("-------run path find")
		path_find.path_find(WD,input, sample,n_core, min_nread,target_info)
	else:
		print("type option Error")

if __name__=="__main__":
	main()
