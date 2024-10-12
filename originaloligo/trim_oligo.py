import argparse
import gzip
import os
import re

from Bio import SeqIO


def flk_trim(WD, input, target_info, targets):

	R1_path="{}/raw/{}_R1.fastq.gz".format(WD, input)
	R2_path="{}/raw/{}_R2.fastq.gz".format(WD, input)

	fastq_1 = SeqIO.parse(gzip.open(R1_path, "rt"), "fastq")
	fastq_2 = SeqIO.parse(gzip.open(R2_path, "rt"), "fastq")

	flks=zip(target_info.index, target_info["flanking_R1"],target_info["flanking_R2"])

	out={}
	out["R1"]={}
	out["R2"]={}
	for target in targets:
		out["R1"][target] = open("{}/trimmed/{}_R1_flkrm.{}.fastq".format(WD, input, target), "w")
		out["R2"][target] = open("{}/trimmed/{}_R2_flkrm.{}.fastq".format(WD, input, target), "w")

	for read1 in fastq_1:
		read2=next(fastq_2)
		read1seq=str(read1.seq)
		read2seq=str(read2.seq)
		
		for target, R1_flk, R2_flk in flks:
			chk_r1=read1seq.startswith(R1_flk)
			if R1_flk=="NA":
				chk_r1=True
				R1_flk=''
			chk_r2=read2seq.startswith(R2_flk)
			if R2_flk=="NA":
				chk_r2=True
				R2_flk=''
			if chk_r1 and chk_r2:
				read1.id="{}_{}".format(read1.id, target)
				read2.id="{}_{}".format(read2.id, target)
				SeqIO.write(read1[len(R1_flk):], out["R1"][target], "fastq")
				SeqIO.write(read2[len(R2_flk):], out["R2"][target], "fastq")

def UID_extract(WD, input, target_info, targets):

		for target in targets:
			fastq_1 = SeqIO.parse(open("{}/trimmed/{}_R1_flkrm.{}.fastq".format(WD, input, target)), "fastq")
			fastq_2 = SeqIO.parse(open("{}/trimmed/{}_R2_flkrm.{}.fastq".format(WD, input, target)), "fastq")
			r1_pattern=target_info["UID_fwd"][target]
			r2_pattern=target_info["UID_rev"][target]

			amp_size=target_info["amplicon_size"][target]
			out_R1 = open("{}/trimmed/{}_R1_processed.{}.fastq".format(WD, input, target), "w")
			out_R2 = open("{}/trimmed/{}_R2_processed.{}.fastq".format(WD, input, target), "w")

			for read1 in fastq_1:
				read2=next(fastq_2)
				filter_check=False
				uids=[]
				for read in [read1, read2]:
					pattern=r1_pattern
					if read.description.split(" ")[-1].split(":")[0]=="2":
						pattern=r2_pattern

					pattern=pattern.replace("N", ".")
					uid=str(read[:len(pattern)].seq)
					uid_qs=read[:len(pattern)].letter_annotations["phred_quality"]
					uids.append(uid)
					if not re.match(pattern, uid) or min(uid_qs)<25:
						filter_check=True
						break

				if filter_check: continue

				read1.id=read1.id+"_"+uids[0]
				read2.id=read2.id+"_"+uids[1]
				R1stt=len(r1_pattern)
				R1end=len(r1_pattern)+amp_size
				R2stt=len(r2_pattern)
				R2end=len(r2_pattern)+amp_size

				read1=read1[R1stt:R1end]
				read2=read2[R2stt:R2end]
				r1_qs=read1.letter_annotations["phred_quality"]
				r2_qs=read2.letter_annotations["phred_quality"]
				r1_qavg=sum(r1_qs)/len(r1_qs)
				r2_qavg=sum(r2_qs)/len(r2_qs)
				if r1_qavg<30 or r2_qavg<30: continue

				SeqIO.write(read1, out_R1, "fastq")
				SeqIO.write(read2, out_R2, "fastq")


def trim_sequence(WD, input, target_info):
	targets=target_info.index

	if not os.path.exists("{}/trimmed".format(WD)):
		os.mkdir("{}/trimmed".format(WD))

	## Trim flanking sequence
	print("Trim flanking sequence")	
	flk_trim(WD, input, target_info, targets)

	## Extract UID sequence
	print("Extract UID sequence")
	UID_extract(WD, input, target_info, targets)

	## merge & cleaning
	os.system("cat {WD}/trimmed/{input}_R1_processed.*.fastq > {WD}/trimmed/{input}_R1.trimmed.fastq".format(WD=WD, input=input))
	os.system("cat {WD}/trimmed/{input}_R2_processed.*.fastq > {WD}/trimmed/{input}_R2.trimmed.fastq".format(WD=WD, input=input))
	for target in targets:
		os.system("rm {}/trimmed/{}_R1_flkrm.{}.fastq".format(WD, input, target))
		os.system("rm {}/trimmed/{}_R2_flkrm.{}.fastq".format(WD, input, target))
		os.system("rm {}/trimmed/{}_R1_processed.{}.fastq".format(WD, input, target))
		os.system("rm {}/trimmed/{}_R2_processed.{}.fastq".format(WD, input, target))
