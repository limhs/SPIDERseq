from multiprocessing import Process
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import copy
import sys
import os
import time

def major_barcode(bcdlist):
	uniq=np.unique(bcdlist, return_counts=True)
	df=pd.DataFrame(zip(uniq[0], uniq[1]), columns=["barcode", "count"]).set_index("barcode")

	major_bcd=df["count"].idxmax()
	major_count=df.loc[major_bcd,"count"]
	total_count=df["count"].sum()

	return major_bcd,major_count, total_count

def maximum_streak(errorlist):
	errorlist=errorlist+[0]
	streaks=[]
	chk=False
	for idx, err in enumerate(errorlist):
		err=float(err)
		if err>0:
			if chk==False:
				chk=True
				stt=idx
				end=idx
			else:
				end=idx
		if err==0:
			if chk==True:
				streaks.append(end-stt+1)
				chk=False
			else:
				continue
	if len(streaks)!=0:	
		return max(streaks)
	else:
		return 0
	
	

def path_find(WD,input, sample,n_core,min_nread, target_info):
 	if not os.path.exists("{}/path/".format(WD)):
		os.mkdir("{}/path/".format(WD))
 	if not os.path.exists("{}/path/{}".format(WD, sample)):
		os.mkdir("{}/path/{}".format(WD, sample))

	filtered=target_info.loc[target_info["Mutation_Type"]=="Substitution", :]
	mut_info={}
	pos_chk={}
	for target, mut, strand, pos in zip(filtered.index, filtered["HGVS_Nomenclature"], filtered["strand"], filtered["mut_pos"]):
		pos=pos.split("..")[0]
		pos_chk[pos]=0
		ref=mut.split(">")[0][-1]
		alt=mut.split(">")[1]
		if strand=="-":
			ref=str(Seq(ref).reverse_complement())
			alt=str(Seq(alt).reverse_complement())
		mut_info[target]=[ref, alt]

	cluster_info=pd.read_csv("{}/cluster/{}.cluster.txt".format(WD, input), sep='\t')
	cids=pd.read_csv("{}/target_for_lineage/{}.txt".format(WD, input), header=None)[0].tolist()
	
	pairbcd_info=pd.read_csv("{}/base_call/{}.{}.tmp.txt".format(WD,sample, min_nread), sep='\t')
	pairbcd_info=pairbcd_info.loc[pairbcd_info["pos"].apply(lambda x: x.split(":")[1]).isin(pos_chk.keys()), :]
	pairbcd_info=pairbcd_info.loc[pairbcd_info["cluster_id"].isin(cids), :].reset_index()
	pairbcd_dedup=pairbcd_info.groupby(["cluster_id","pair"])["base"].apply(list).reset_index()	
	
	color_code=pd.read_csv("/Data2/JSY/PCR_dedup/color_reds.csv", index_col=0)
	color_code=dict(zip(color_code["a"], color_code["col"]))
	##### process start
	sideswap={}
	sideswap["UMI_R1"]="UMI_R2"
	sideswap["UMI_R2"]="UMI_R1"
	out_df=pd.DataFrame()
	count_label=open("{}/path/{}/count.tsv".format(WD, sample), "w")
	for idx in cluster_info.index:
		start_side=cluster_info.loc[idx, "Segment"]
		start_bcd=cluster_info.loc[idx, "origin"]
		
		pair_df=pd.DataFrame(cluster_info.loc[idx, "pair_list"].split(";"), columns=["pair_list"])
		pair_df["R1"]=pair_df["pair_list"].apply(lambda x: x.split("_")[0])
		pair_df["R2"]=pair_df["pair_list"].apply(lambda x: x.split("_")[1])
		ff=pair_df.groupby("R1")["R2"].apply(list).to_dict()
		rr=pair_df.groupby("R2")["R1"].apply(list).to_dict()
		pairbcd_tmp=pairbcd_dedup.loc[(pairbcd_dedup["pair"].isin(list(pair_df["pair_list"]))), :]
		pairbcd_tmp=pairbcd_tmp.set_index("pair")
		if pairbcd_tmp.shape[0]==0: continue
		cluster_id=list(pairbcd_tmp["cluster_id"])[0]
		target=cluster_id.split("_")[0]	
		d={}
		d["UMI_R1"]=ff
		d["UMI_R2"]=rr
		origin_base=mut_info[target][0]

		path_list=[]
		barcode_list=[]
		count_list=[]
		total_list=[]
		##
		pairs_list=[]

		while True:
			input=start_bcd
			side=start_side
			prev=''
			if len(d[side][input])==0:
				break
			path=[[start_bcd, start_side]]
			while True:
				next_list=d[side][input]
				next_side=sideswap[side]
				next_bcd=next_list[0]
				if next_bcd==prev:
					next_bcd=next_list[1]
				path.append([next_bcd,next_side])
				prev=input
				input=next_bcd
				side=next_side
				if len(d[next_side][input])==1:
					break

			for i in range(len(path)-1):
				ii=len(path)-i-1
				bcd=path[ii][0]
				side=path[ii][1]
				prev=path[ii-1][0]
				prev_side=d[path[ii-1][1]]
				if len(d[side][bcd])==1:
					prev_side[prev].remove(bcd)
				if prev_side[prev]==1:
					break
			route=[x[0] for x in path]
			path_list.append("-".join(route))
			barcodes=[]
			counts=[]
			totals=[]
			pairs=["{}({})".format(start_bcd, start_side)]
			for idx in range(len(route)-1):
				p1=route[idx]
				p2=route[idx+1]
				pair=p1+"_"+p2
				if (idx%2==1 and start_side=="UMI_R1") or (idx%2==0 and start_side=="UMI_R2"):
					pair=p2+"_"+p1
				pairs.append(pair)
				try:
					bcds=pairbcd_tmp.loc[pair, "base"]
					major_bcd=major_barcode(bcds)[0]
					cnt=bcds.count(origin_base)
					# print major_bcd
					# print origin_base
					# print cnt
					# print ""
					total=len(bcds)
					cnt=round(float(cnt)/total,2)
				except KeyError:
					major_bcd="unidentified"
					cnt=round(float(0),2)
					total=0
				barcodes.append(major_bcd)
				counts.append(str(1.00-cnt))
				totals.append(str(total))

			for p, c, b in zip(pairs[1:], counts, barcodes):
				try:
					count_label.write("{}\t{}\t{}\n".format(p, round(float(c),2), color_code[round(float(c),2)]))
				except KeyError: 
					count_label.write("{}\t{}\t{}\n".format(p, round(c,2), "NA"))
			pairs=[sample, target ]+pairs
			
			barcode_list.append("-".join(barcodes))
			count_list.append(";".join(counts))
			total_list.append(";".join(totals))
			pairs_list.append('='.join(pairs))

		path_df=pd.DataFrame(path_list, columns=["path"])
		path_df["start_segment"]=start_side
		path_df["cluster_id"]=cluster_id
		path_df["barcodes"]=barcode_list
		path_df["error_freq"]=count_list
		path_df["totals"]=total_list
		path_df["nNodes"]=[len(x.split("-")) for x in barcode_list]
		path_df["pairs"]=pairs_list
		out_df=out_df.append(path_df)
	

	out_df["number_of_contents"]=out_df["barcodes"].apply(lambda x: len(set(x.split("-"))-{"unidentified"}))
	out_df["unidentified"]=out_df["barcodes"].str.contains("unidentified")
	out_df=out_df.loc[out_df["unidentified"]==False, :]


	lineage=out_df["pairs"].apply(lambda x: '\t'.join(x.split("=")+["NA"]*(10-len(x.split('=')))))
	pd.DataFrame(lineage).to_csv("{}/path/{}/lineage.filtered.tsv".format(WD, sample), index=False, header=False)
	out_df["ErrGT0.5_ratio"]=out_df["error_freq"].apply(lambda x: np.count_nonzero(np.array(x.split(";")).astype(float)>0.5)/float(len(x.split(";")) ))
	out_df["maxerr_streak"]=out_df["error_freq"].apply(lambda x: maximum_streak(x.split(";")) )
	out_df["maximum_errorratio"]=out_df["error_freq"].apply(lambda x: max(np.array(x.split(";")).astype(float)))
	out_df[["cluster_id",  "start_segment", "path", "barcodes","error_freq", "totals", "number_of_contents","nNodes","ErrGT0.5_ratio", "maxerr_streak"]].to_csv("{}/path/{}/path.filtered.tsv".format(WD, sample), sep='\t', index=False)
