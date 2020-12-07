#for capture

import os
import argparse

def bwa(WD, sample, n_core, bwa_ref):

	if not os.path.exists("{}/align".format(WD)):
		os.mkdir("{}/align".format(WD))


	os.system("bwa mem -t {n_core} {bw_idx} {WD}/trimmed/{sample}_R1.trimmed.fastq > {WD}/align/{sample}_R1.sam".format(WD=WD, sample=sample, bw_idx=bwa_ref, n_core=n_core))
	os.system( "samtools import {bw_idx} {WD}/align/{sample}_R1.sam {WD}/align/{sample}_R1.bam".format(WD=WD, sample=sample, bw_idx=bwa_ref) )
	os.system( "samtools sort {WD}/align/{sample}_R1.bam -o {WD}/align/{sample}_R1.sorted.bam".format(WD=WD, sample=sample, bw_idx=bwa_ref) )
	os.system( "samtools index {WD}/align/{sample}_R1.sorted.bam".format(WD=WD, sample=sample) )

	os.system("bwa mem -t {n_core} {bw_idx} {WD}/trimmed/{sample}_R2.trimmed.fastq > {WD}/align/{sample}_R2.sam".format(WD=WD, sample=sample, bw_idx=bwa_ref, n_core=n_core))
	os.system( "samtools import {bw_idx} {WD}/align/{sample}_R2.sam {WD}/align/{sample}_R2.bam".format(WD=WD, sample=sample, bw_idx=bwa_ref) )
	os.system( "samtools sort {WD}/align/{sample}_R2.bam -o {WD}/align/{sample}_R2.sorted.bam".format(WD=WD, sample=sample, bw_idx=bwa_ref) )
	os.system( "samtools index {WD}/align/{sample}_R2.sorted.bam".format(WD=WD, sample=sample) )

if __name__=="__main__":
	usage = "python align.py [options]\n"
	parser = argparse.ArgumentParser(description = usage)
	parser.add_argument("-s", "--sample_prefix", type=str, help="prefix of input & output files, input fastq is must be named as <prefix>_<R1 or R2>.fastq.gz")
	parser.add_argument("-W" , "--work_dir", type=str, default="./", help="Working directory (base of output files)")
	parser.add_argument("-b" ,"--bwa_ref", type=str , help="<dir>/prefix of bwa reference")
	parser.add_argument("-t" , "--threads", type=int, default=1, help="Number of threads")
	args = parser.parse_args()
	n_core=args.threads
	sample=args.sample_prefix
	WD=args.work_dir
	bwa_ref=args.bwa_ref
	bwa(WD, sample, n_core, bwa_ref)
