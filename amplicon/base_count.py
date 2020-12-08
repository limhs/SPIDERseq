import os
from multiprocessing import Process, Queue

import numpy as np
import time
import pandas as pd
import pysam
from conventional_call import basic_base_count, uid_based_call
from utils import get_loci, get_vcf_info, output_merge, run_multi, set_pos, get_refseq
import re

def cluster_based_call(WD,input, sample,n_core, target_info, min_base_ratio,  min_nread):

	def thread_function(WD, input, sample, n, target_info, df_chunk, min_base_ratio, pair_info):
		out_df=[]
		df_chunk_out=pd.DataFrame()
		tmp_out=open("{}/base_call/{}.{}.{}.tmp.txt".format(WD, sample, str(min_nread),str(n)), "w")
		tmp_out.write("cluster_id\tpair\tpos\tbase\n")
		for target, chr, poslist in loci:
			poslist=set_pos(poslist)

			cluster_filtered=df_chunk.loc[df_chunk["target"]==target, :]
			cluster_filtered.loc[:,"cluster_id"]=target+"_"+str(n)+":"+(cluster_filtered.index+1).astype("str")
			if cluster_filtered.shape[0]==0: continue
			pair_to_cid={}
			for cid, pairlist in zip(cluster_filtered["cluster_id"], cluster_filtered["pair_list"]):
				for pair in pairlist.split(";"):
					pair_to_cid[pair]=cid

			pair_filtered=pair_info.loc[pair_info["pair"].isin(pair_to_cid.keys()), :]
			id_to_pair={}
			id_df=pair_filtered.groupby("pair")["id"].apply(list)

			for pair, idlist in zip(id_df.index, id_df):
				for id in idlist:
					id_to_pair[id]=pair
			
			dfdic={}
			for pos in poslist:
				locus_id=target+"_"+str(pos)
				dfdic[locus_id]=pd.DataFrame(cluster_filtered["cluster_id"])
				for b in "ATCGN":
					dfdic[locus_id][b]=0
				dfdic[locus_id]=dfdic[locus_id].set_index("cluster_id")

			for seg in ["R1", "R2"]:
				bamfile=pysam.AlignmentFile("{}/align/{}_{}.sorted.bam".format(WD, input, seg), "rb")
				for read in bamfile.fetch(chr, min(poslist), max(poslist)+1):
					id=read.qname.split("_")[0]
					
					try:
						pair=id_to_pair[id]
						cid=pair_to_cid[pair]
					except KeyError:
						continue
					
					for pos in poslist:
						locus_id=target+"_"+str(pos)
						try:
							idx=read.get_reference_positions(full_length=True).index(pos-1)
						except ValueError:
							continue
					
						base=read.query_sequence[idx]
						dfdic[locus_id].loc[cid,base]+=1
						tmp_out.write(cid+'\t'+pair+'\t'+chr+":"+str(pos)+'\t'+base+'\n')
			merged=output_merge(target,chr, poslist, dfdic, out_df, min_base_ratio)
			out_df=merged[0]
			df_chunk_out=df_chunk_out.append(merged[1])
		out_df=pd.DataFrame(out_df, columns="locus_id\tA\tT\tC\tG\tN".split('\t')).set_index("locus_id")
		out_df.to_csv("{}/base_call/{}.cluster.{}.{}.base_count.txt".format(WD, sample, str(min_nread),str(n)), sep='\t')
		df_chunk_out.to_csv("{}/base_call/{}.cluster.{}.{}.rawchunk.txt".format(WD, sample, str(min_nread),str(n)), sep='\t')

	loci=get_loci(target_info)
	cluster_info=pd.read_csv("{}/cluster/{}.cluster.txt".format(WD, input), sep='\t')
	cluster_info=cluster_info.loc[cluster_info["nPairs"]>=2, :]
	cluster_info=cluster_info.loc[cluster_info["nReads"]>=min_nread, :]
	pair_info=pd.read_csv("{}/pair/{}_barcode_pair.txt".format(WD, input), sep="\t")
	
	run_multi(thread_function, n_core, cluster_info, WD, input,sample, target_info, min_base_ratio,other_args=(pair_info,))
	out_df=pd.DataFrame()
	tmp_df=pd.DataFrame()
	raw_df=pd.DataFrame()
	for n in range(1, n_core+1):
		while True:
			try:
				out_chunk=pd.read_csv("{}/base_call/{}.cluster.{}.{}.base_count.txt".format(WD, sample, str(min_nread), str(n)), sep='\t', index_col="locus_id")
				tmp_chunk=pd.read_csv("{}/base_call/{}.{}.{}.tmp.txt".format(WD, sample, str(min_nread),str(n)), sep='\t')
				raw_chunk=pd.read_csv("{}/base_call/{}.cluster.{}.{}.rawchunk.txt".format(WD, sample, str(min_nread),str(n)), sep='\t')
			except IOError:
				time.sleep(30)
				continue
			break
		out_df=out_df.add(out_chunk, fill_value=0)
		tmp_df=tmp_df.append(tmp_chunk)
		raw_df=raw_df.append(raw_chunk)
		os.system("rm {}/base_call/{}.cluster.{}.{}.base_count.txt".format(WD, sample, str(min_nread), str(n)) )
		os.system("rm {}/base_call/{}.{}.{}.tmp.txt".format(WD, sample, str(min_nread),str(n)) )
		os.system("rm {}/base_call/{}.cluster.{}.{}.rawchunk.txt".format(WD, sample, str(min_nread),str(n)) )

	out_df.astype("int").to_csv("{}/base_call/{}.cluster.{}.base_count.txt".format(WD, sample, str(min_nread)), sep='\t')
	tmp_df.to_csv("{}/base_call/{}.{}.tmp.txt".format(WD, sample, str(min_nread)  ), sep='\t')
	raw_df.to_csv("{}/base_call/{}.cluster.{}.rawchunk.txt".format(WD, sample, str(min_nread)), sep='\t')

def origin_based_call(WD,input, sample,n_core, target_info, min_base_ratio, min_nread):

	def thread_function(WD, input, sample, target_info, Q, df_chunk, min_base_ratio, pair_info):
		out_df=[]
		for target, chr, poslist in loci:
			poslist=set_pos(poslist)

			cluster_filtered=df_chunk.loc[df_chunk["target"]==target, :]
			cluster_filtered.loc[:,"cluster_id"]=target+(cluster_filtered.index+1).astype("str")
			if cluster_filtered.shape[0]==0: continue
			
			pair_to_cid={}
			for cid, origin, pairlist in zip(cluster_filtered["cluster_id"], cluster_filtered["origin"], cluster_filtered["pair_list"]):
				for pair in pairlist.split(";"):
					if origin in pair:
						pair_to_cid[pair]=cid
			pair_filtered=pair_info.loc[pair_info["pair"].isin(pair_to_cid.keys()), :]

			id_to_pair={}
			id_df=pair_filtered.groupby("pair")["id"].apply(list)
			for pair, idlist in zip(id_df.index, id_df):
				for id in idlist:
					id_to_pair[id]=pair
			
			dfdic={}
			for pos in poslist:
				locus_id=target+"_"+str(pos)
				dfdic[locus_id]=pd.DataFrame(cluster_filtered["cluster_id"])
				for b in "ATCGN":
					dfdic[locus_id][b]=0
				dfdic[locus_id]=dfdic[locus_id].set_index("cluster_id")

			for seg in ["R1", "R2"]:
				bamfile=pysam.AlignmentFile("{}/align/{}_{}.sorted.bam".format(WD, input, seg), "rb")
				for read in bamfile.fetch(chr, min(poslist), max(poslist)+1):
					id=read.qname.split("_")[0]
					
					try:
						pair=id_to_pair[id]
						cid=pair_to_cid[pair]
					except KeyError:
						continue
					
					for pos in poslist:
						locus_id=target+"_"+str(pos)
						try:
							idx=read.get_reference_positions(full_length=True).index(pos-1)
						except ValueError:
							continue
					
						base=read.query_sequence[idx]
						dfdic[locus_id].loc[cid,base]+=1

			out_df=output_merge(target,chr, poslist, dfdic, out_df, min_base_ratio)
		out_df=pd.DataFrame(out_df, columns="locus_id\tA\tT\tC\tG\tN".split('\t')).set_index("locus_id")
		Q.put(out_df)

	loci=get_loci(target_info)
	cluster_info=pd.read_csv("{}/cluster/{}.cluster.txt".format(WD, input), sep='\t')
	cluster_info=cluster_info.loc[cluster_info["nReads"]>=min_nread, :]
	pair_info=pd.read_csv("{}/pair/{}_barcode_pair.txt".format(WD, input), sep="\t")
	out_df=run_multi(thread_function, n_core, cluster_info, WD, input,sample, target_info, min_base_ratio,other_args=(pair_info,))
	out_df.to_csv("{}/base_call/{}.origin.{}.base_count.txt".format(WD, sample, str(min_nread)), sep='\t')


def cluster_based_indel(WD,input, sample,n_core, target_info, min_base_ratio, min_nread,vcfpath, refpath):

	def thread_function(WD, input, sample, n, target_info, df_chunk, min_base_ratio, vcf_info, pair_info):
		print n
		out_df=[]
		for target, chr, poslist in loci:
			poslist=set_pos(poslist)

			cluster_filtered=df_chunk[df_chunk.loc[:,"target"]==target]
			cluster_filtered.loc[:,"cluster_id"]=target+"_"+str(n)+":"+(cluster_filtered.index+1).astype("str")
			if cluster_filtered.shape[0]==0: continue

			pair_to_cid={}
			for cid, pairlist in zip(cluster_filtered["cluster_id"], cluster_filtered["pair_list"]):
				for pair in pairlist.split(";"):
					pair_to_cid[pair]=cid

			pair_filtered=pair_info.loc[pair_info["pair"].isin(pair_to_cid.keys()), :]

			id_to_pair={}
			id_df=pair_filtered.groupby("pair")["id"].apply(list)

			for pair, idlist in zip(id_df.index, id_df):
				for id in idlist:
					id_to_pair[id]=pair

			vcf_filtered=vcf_info[ (vcf_info["chr"]==chr) &  (vcf_info["pos"].isin( range( min(poslist) - 200, max(poslist) + 200 )  ))]
			vcf_filtered=vcf_filtered.set_index("mut")
			if vcf_filtered.shape[0]==0: continue
			dfdic={}
			for mut in vcf_filtered.index:
				dfdic[mut]=pd.DataFrame(cluster_filtered["cluster_id"])
				for muttype in ["WT", "INDEL"]:
					dfdic[mut][muttype]=0
				dfdic[mut]=dfdic[mut].set_index("cluster_id")

			
			for seg in ["R1", "R2"]:
				bamfile=pysam.AlignmentFile("{}/align/{}_{}.sorted.bam".format(WD, input, seg), "rb")
				for read in bamfile.fetch(chr, min(poslist), max(poslist)+1):
					id=read.qname.split("_")[0]
					
					try:
						pair=id_to_pair[id]
						cid=pair_to_cid[pair]
					except KeyError:
						continue
					
					for mut in dfdic.keys():
						vcfchr=vcf_filtered.loc[mut, "chr"]
						if vcfchr!=chr: continue
						refquery=vcf_filtered.loc[mut, "ref_query"]
						altquery=vcf_filtered.loc[mut, "alt_query"]
						qseq=read.query_sequence
						if altquery in qseq:
							muttype="INDEL"
						elif refquery in qseq:
							muttype="WT"
						else:
							continue
						dfdic[mut].loc[cid, muttype]+=1
						
						
			
			for mut in dfdic.keys():
				df=dfdic[mut]
				df["total"]=0
				for muttype in ["WT", "INDEL"]:
					df["total"]+=df[muttype]
				df=df[df["total"]>0]

				df.loc[:,"max_type"]=df.loc[:,df.columns.isin(["WT", "INDEL"])].idxmax(axis=1)
				df.loc[:,"max_type_count"]=df.loc[:,df.columns.isin(["WT", "INDEL"])].max(axis=1)
				df.loc[:,"max_type_ratio"]=df.loc[:,"max_type_count"].astype("float")/df.loc[:,"total"]

				df_filtered=df.loc[df["max_type_ratio"] > min_base_ratio, :]

				max_type_list=list(df_filtered["max_type"])
				tmp=[mut]
				for b in ["WT", "INDEL"]:
					tmp.append(max_type_list.count(b))
				out_df.append(tmp)
		out_df=pd.DataFrame(out_df, columns="mutation\tWT\tINDEL".split('\t')).set_index("mutation")
		out_df.to_csv("{}/indel_call/{}.cluster.{}.{}.base_count.txt".format(WD, sample, str(min_nread),str(n)), sep='\t')
		
	print refpath
	loci=get_loci(target_info)
	cluster_info=pd.read_csv("{}/cluster/{}.cluster.txt".format(WD, sample), sep='\t')
	cluster_info=cluster_info.loc[cluster_info["nReads"]>=min_nread, :]
	pair_info=pd.read_csv("{}/pair/{}_barcode_pair.txt".format(WD, sample), sep="\t")

	refseq=get_refseq(refpath)
	vcf_info=get_vcf_info(vcfpath, refseq)

	run_multi(thread_function, n_core, cluster_info, WD, input,sample, target_info, min_base_ratio,other_args=(vcf_info, pair_info))
	out_df=pd.DataFrame()
	
	for n in range(1, n_core+1):
		while True:
			try:
				out_chunk=pd.read_csv("{}/indel_call/{}.cluster.{}.{}.indel.txt".format(WD, sample, str(min_nread), str(n)), sep='\t', index_col="mutation")
				
			except IOError:
				time.sleep(30)
				continue
			break
		out_df=out_df.add(out_chunk, fill_value=0)
		os.system("rm {}/indel_call/{}.cluster.{}.{}.indel.txt".format(WD, sample, str(min_nread),str(n)) )
		
	out_df.astype("int").to_csv("{}/indel_call/{}.cluster.{}.indel.txt".format(WD, sample, min_nread), sep='\t')
	