from Bio import SeqIO
import pandas as pd
import numpy as np
from multiprocessing import Process, Queue
import os
import re
import numpy

def HD(a,b):
	d=0
	for n in range(len(a)):
		if a[n]!=b[n]:
			d+=1
	return d
	


def specificity(WD,input, sample,n_core, target_info):
	if not os.path.exists("{}/specificity".format(WD)):
		os.mkdir("{}/specificity".format(WD))
	def function(WD,loci, input, sample, cnumber, target_info, df_chunk, pair_info):
		out_df=pd.DataFrame()
		tmp_out=open("{}/specificity/{}.{}.pair-barcode.txt".format(WD, sample, str(cnumber)), "w")
		tmp_out.write('cluster_id\tpair\tbarcode\n')
		for target, N1reg, N2reg in loci:
			cluster_filtered=df_chunk.loc[df_chunk["target"]==target,:]
			cluster_filtered["cluster_id"]=target+"_"+cnumber+"_"+(cluster_filtered.index+1).astype("str")
			cluster_filtered["size"]=cluster_filtered["UID_R1"].apply(lambda x: len(x.split(";"))) + cluster_filtered["UID_R2"].apply(lambda x: len(x.split(";")))
			cluster_filtered=cluster_filtered.set_index("cluster_id")


			## for pair to cid convert
			pair_to_cid={}
			for cid, pairlist in zip(cluster_filtered.index, cluster_filtered["pair_list"]):
				for pair in pairlist.split(";"):
					pair_to_cid[pair]=cid

			id_to_pair={}
			pair_filtered=pair_info.loc[pair_info["pair"].isin(pair_to_cid.keys()), :]
			id_df=pair_filtered.groupby("pair")["id"].apply(list)
			## for id to pair convert
			for pair, idlist in zip(id_df.index, id_df):
				for id in idlist:
					id_to_pair[id]=pair

			barcode_count={}
			for cid in cluster_filtered.index:
				barcode_count[cid]={}

			for seg in ["R1", "R2"]:
				fastq=SeqIO.parse(open("{}/trimmed/{}_{}.trimmed.fastq".format(WD, input, seg)), "fastq")
				for read in fastq:
					id=read.id.split("_")[0]
					seq=read.seq
					if seg=="R2":
						seq=read.seq.reverse_complement()
					seq=str(seq)
					
					try:
						pair=id_to_pair[id]
						cid=pair_to_cid[pair]
					except KeyError:
						continue
					N1search=re.findall(N1reg, seq)
					N2search=re.findall(N2reg, seq)
					if len(N1search)!=1 or len(N2search)!=1 : continue
					bcd=N1search[0][4:-4]+N2search[0][4:-4]
					try:
						barcode_count[cid][bcd]+=1
					except KeyError:
						barcode_count[cid][bcd]=1
					tmp_out.write(cid+'\t'+pair+'\t'+bcd+'\n')

			barcode_count_corrected={}
			mismatch={}
			major_barcode_dic={}
			for cid in barcode_count.keys():
				if len(barcode_count[cid])==0: continue
				count_df=pd.DataFrame(zip(barcode_count[cid].keys(), barcode_count[cid].values()), columns=["barcode", "count"]).set_index("barcode")
				major_barcode=count_df["count"].idxmax()
				major_barcode_dic[cid]=major_barcode

				barcode_count_corrected[cid]={}
				mismatch[cid]={1:0, 2:0, "diff":0}
				barcode_count_corrected[cid][major_barcode]=barcode_count[cid][major_barcode]

				for bcd in count_df.index:
					if bcd==major_barcode: continue
					if HD(major_barcode,bcd) < 3:
						barcode_count_corrected[cid][major_barcode]+=barcode_count[cid][bcd]
						mismatch[cid][HD(major_barcode,bcd)]+=1
					else:
						barcode_count_corrected[cid][bcd]=barcode_count[cid][bcd]
						mismatch[cid]["diff"]+=1

			df=pd.DataFrame(cluster_filtered.index)
			df=df.loc[df["cluster_id"].isin(barcode_count_corrected.keys())]
			df["barcodes"]=df["cluster_id"].apply(lambda x : ";".join(barcode_count[x].keys()))
			df["barcodes_corrected"]=df["cluster_id"].apply(lambda x : ";".join(barcode_count_corrected[x].keys()))
			df["total"]=df["cluster_id"].apply(lambda x: sum(barcode_count[x].values()))
			df=df.loc[df["total"]>0,:]

			df["major_barcode"]=df["cluster_id"].apply(lambda x : major_barcode_dic[x])
			df["major_count"]=df["cluster_id"].apply(lambda x : barcode_count[x][major_barcode_dic[x]])
			df["major_count_corrected"]=df["cluster_id"].apply(lambda x : barcode_count_corrected[x][major_barcode_dic[x]])
			df["specificity"]=df["major_count"].astype("float")/df["total"]
			df["specificity_corrected"]=df["major_count_corrected"].astype("float")/df["total"]
			df["target"]=target
			for term in ["size", "nReads", "nPairs", "UID_R1", "UID_R2"]:
				df[term]=df["cluster_id"].apply(lambda x : cluster_filtered.loc[x, term])

			df["M1_corrected"]=df["cluster_id"].apply(lambda x : mismatch[x][1])
			df["M2_corrected"]=df["cluster_id"].apply(lambda x : mismatch[x][2])
			df["Not_corrected"]=df["cluster_id"].apply(lambda x : mismatch[x]["diff"])
			df=df.set_index("cluster_id")
			out_df=out_df.append(df)
		out_df.to_csv("{}/specificity/{}.{}.specificity.txt".format(WD, sample, str(cnumber)), sep='\t')
		# Q.put(out_df)	

	targets=target_info.index
	N1reg=target_info["N1reg"]
	N2reg=target_info["N2reg"]
	loci=zip(targets, N1reg, N2reg)
	cluster_info=pd.read_csv("{}/cluster/{}.cluster.txt".format(WD, sample), sep='\t')
	pair_info=pd.read_csv("{}/pair/{}_barcode_pair.txt".format(WD, sample), sep="\t")
	Ps=[]
	unit=len(cluster_info.index)/n_core
	for n in range(1,n_core+1):
		stt=(n-1)*unit
		end=n*unit-1
		if n==n_core:
			end=(n+1)*unit-1
		proc=Process(target=function, args=(WD,loci, input, sample, str(n),  target_info, cluster_info.loc[stt:end,], pair_info) )
		Ps.append(proc)
		proc.start()
	for p in Ps:
	    p.join()
	
	out_df=pd.DataFrame()
	tmp_df=pd.DataFrame()
	for n in range(1,n_core+1):
		df=pd.read_csv("{}/specificity/{}.{}.specificity.txt".format(WD, sample, str(n)), sep='\t')
		os.system("rm {}/specificity/{}.{}.specificity.txt".format(WD, sample, str(n)))
		out_df=out_df.append(df)

		tmp=pd.read_csv("{}/specificity/{}.{}.pair-barcode.txt".format(WD, sample, str(n)), sep='\t')
		os.system("rm {}/specificity/{}.{}.pair-barcode.txt".format(WD, sample, str(n)))
		tmp_df=tmp_df.append(tmp)
		
	out_df.to_csv("{}/specificity/{}.specificity.txt".format(WD, sample), sep='\t', index=False)
	tmp_df.to_csv("{}/specificity/{}.pair-barcode.txt".format(WD, sample), sep='\t', index=False)
