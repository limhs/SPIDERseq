#!/usr/bin/env python
import argparse
import gzip
import os
import re

from Bio import SeqIO
from Bio.Seq import Seq

def UID_extract(WD, input):

	R1_path="{}/raw/{}_R1.fastq.gz".format(WD, input)
	R2_path="{}/raw/{}_R2.fastq.gz".format(WD, input)
	fastq_1 = SeqIO.parse(gzip.open(R1_path, "rt"), "fastq")
	fastq_2 = SeqIO.parse(gzip.open(R2_path, "rt"), "fastq")

	out_R1 = open("{}/trimmed/{}_R1.trimmed.fastq".format(WD, input), "w")
	out_R2 = open("{}/trimmed/{}_R2.trimmed.fastq".format(WD, input), "w")

	for read1 in fastq_1:
		read2=next(fastq_2)
		spl=read1.description.split(':')[-1].split("+")
		uid1=spl[1][3:]
		uid2=spl[0][:-3]
		if "N" in uid1+uid2: continue

		read1.id=read1.id+"_"+uid1
		read2.id=read2.id+"_"+uid2

		r1_qs=read1.letter_annotations["phred_quality"]
		r2_qs=read2.letter_annotations["phred_quality"]
		r1_qavg=sum(r1_qs)/len(r1_qs)
		r2_qavg=sum(r2_qs)/len(r2_qs)
		if r1_qavg<30 or r2_qavg<30: continue

		SeqIO.write(read1, out_R1, "fastq")
		SeqIO.write(read2, out_R2, "fastq")


def trim_sequence(WD, input):

	if not os.path.exists("{}/trimmed".format(WD)):
		os.mkdir("{}/trimmed".format(WD))

	print("Extract UID sequence")
	UID_extract(WD, input)



