#!/usr/bin/env python

import os
import argparse

def bwa(WD, sample, threads, bwa_ref):

	if not os.path.exists("{}/align".format(WD)):
		os.mkdir("{}/align".format(WD))

	os.system( "bwa mem -t {threads} {bw_idx} {WD}/trimmed/{sample}_R1.trimmed.fastq > {WD}/align/{sample}_R1.sam".format(WD=WD, sample=sample, bw_idx=bwa_ref, threads=threads))
	os.system( "samtools import {bw_idx} {WD}/align/{sample}_R1.sam {WD}/align/{sample}_R1.bam".format(WD=WD, sample=sample, bw_idx=bwa_ref) )
	os.system( "samtools sort {WD}/align/{sample}_R1.bam -o {WD}/align/{sample}_R1.sorted.bam".format(WD=WD, sample=sample, bw_idx=bwa_ref) )
	os.system( "samtools index {WD}/align/{sample}_R1.sorted.bam".format(WD=WD, sample=sample) )

	os.system( "bwa mem -t {threads} {bw_idx} {WD}/trimmed/{sample}_R2.trimmed.fastq > {WD}/align/{sample}_R2.sam".format(WD=WD, sample=sample, bw_idx=bwa_ref, threads=threads))
	os.system( "samtools import {bw_idx} {WD}/align/{sample}_R2.sam {WD}/align/{sample}_R2.bam".format(WD=WD, sample=sample, bw_idx=bwa_ref) )
	os.system( "samtools sort {WD}/align/{sample}_R2.bam -o {WD}/align/{sample}_R2.sorted.bam".format(WD=WD, sample=sample, bw_idx=bwa_ref) )
	os.system( "samtools index {WD}/align/{sample}_R2.sorted.bam".format(WD=WD, sample=sample) )

