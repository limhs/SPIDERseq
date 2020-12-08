import numpy as np
import pandas as pd
import pysam
import os
from utils import get_loci, get_vcf_info, output_merge, run_multi, set_pos
import time

def uid_based_call(WD,input, sample, target_info,n_core, min_base_ratio):
	def thread_function(WD, input, sample,n, target_info, df_chunk, min_base_ratio, idTOpair):
		out_df=[]
		for target, chr, poslist in loci:
			poslist=set_pos(poslist)

			filtered=df_chunk[df_chunk["target"]==target]
			if filtered.shape[0]==0: continue
			pairchk=dict(zip(filtered["pair_id"], filtered["target"]))
			dfdic={}
			for pos in poslist:
				locus_id=target+"_"+str(pos)
				dfdic[locus_id]=pd.DataFrame({"pair_id":list(filtered["pair_id"])})
				for b in "ATCGN":
					dfdic[locus_id][b]=0
				dfdic[locus_id]=dfdic[locus_id].set_index("pair_id")

			for seg in ["R1", "R2"]:
				bamfile=pysam.AlignmentFile("{}/align/{}_{}.sorted.bam".format(WD, input, seg), "rb")
				for read in bamfile.fetch(chr, pos, pos+1):
					id=read.qname.split("_")[0]

					try:
						pair_id=idTOpair[id]
						tmp=pairchk[pair_id]
					except KeyError:
						continue
					for pos in poslist:
						locus_id=target+"_"+str(pos)
						try:
							idx=read.get_reference_positions(full_length=True).index(pos-1)
						except ValueError:
							continue
						base=read.query_sequence[idx]
						dfdic[locus_id].loc[pair_id,base]+=1
				
			merged=output_merge(target,chr, poslist, dfdic, out_df, min_base_ratio)
			out_df=merged[0]
		out_df=pd.DataFrame(out_df, columns="locus_id\tA\tT\tC\tG\tN".split('\t')).set_index("locus_id")
		out_df.to_csv("{}/base_call/{}.uid.{}.base_count.txt".format(WD, sample, str(n)), sep='\t')

	loci=get_loci(target_info)
	pair_info=pd.read_csv("{}/pair/{}_barcode_pair.txt".format(WD, input), sep='\t')
	pair_info["pair_id"]=pair_info["pair"]+"_"+pair_info["target"]
	idTOpair=dict(zip(pair_info["id"], pair_info["pair_id"]))
	pair_dedup=pair_info.loc[:, ["pair_id", "target"]].drop_duplicates().reset_index()

	run_multi(thread_function, n_core, pair_dedup, WD, input,sample, target_info, min_base_ratio, other_args=(idTOpair,))
	out_df=pd.DataFrame()
	for n in range(1, n_core+1):
		while True:
			try:
				out_chunk=pd.read_csv("{}/base_call/{}.uid.{}.base_count.txt".format(WD, sample, str(n)), sep='\t', index_col="locus_id")
			except IOError:
				time.sleep(30)
				continue
			break
		out_df=out_df.add(out_chunk, fill_value=0)
		os.system("rm {}/base_call/{}.uid.{}.base_count.txt".format(WD, sample, str(n)) )
	out_df.astype("int").to_csv("{}/base_call/{}.uid.base_count.txt".format(WD, sample), sep='\t')


def basic_base_count(WD, input, sample, target_info):
	loci=get_loci(target_info)
	out_df=[]
	for target, chr, poslist in loci:
		poslist=set_pos(poslist)
		base_count={}
		for pos in poslist:
			locus_id=target+"_"+str(pos)
			base_count[locus_id]={}
			for b in "ATCGN":
				base_count[locus_id][b]=0

		for seg in ["R1", "R2"]:
			bamfile=pysam.AlignmentFile("{}/align/{}_{}.sorted.bam".format(WD, input, seg), "rb")
			for read in bamfile.fetch(chr, min(poslist), max(poslist)+1):
				cigar=read.cigarstring
				mapq=read.mapping_quality
				if mapq<55 or "S" in cigar: continue
				for pos in poslist:
					locus_id=target+"_"+str(pos)
					try:
						idx=read.get_reference_positions(full_length=True).index(pos-1)
					except ValueError:
						continue
					base=read.query_sequence[idx]
					base_count[locus_id][base]+=1

		for pos in poslist:
			locus_id=target+"_"+str(pos)
			tmp=[locus_id]
			for b in "ATCGN":
				tmp.append(base_count[locus_id][b])
			if sum(tmp[1:])==0:continue
			out_df.append(tmp)
	out_df=pd.DataFrame(out_df, columns="locus_id\tA\tT\tC\tG\tN".split('\t')).set_index("locus_id")
	out_df.astype("int").to_csv("{}/base_call/{}.basic.base_count.txt".format(WD, sample), sep='\t')