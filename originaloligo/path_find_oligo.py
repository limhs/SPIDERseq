from multiprocessing import Process
import pandas as pd
import numpy as np
import copy
import sys
import os
import random
import time

def major_barcode(major, bcdlist):
	uniq=np.unique(bcdlist, return_counts=True)
	df=pd.DataFrame(zip(uniq[0], uniq[1]), columns=["barcode", "count"]).set_index("barcode")

	major_bcd=df["count"].idxmax()
	major_count=bcdlist.count(major)
	total_count=len(bcdlist)

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


	

def path_find(WD,input, sample,n_core, target_info):
 	if not os.path.exists("{}/path/".format(WD)):
		os.mkdir("{}/path/".format(WD))
 	if not os.path.exists("{}/path/{}/".format(WD, sample)):
		os.mkdir("{}/path/{}/".format(WD, sample))

	specificity_info=pd.read_csv("{}/specificity/{}.specificity.txt".format(WD, sample), sep='\t')
	MB_convert=dict(zip(specificity_info["cluster_id"], specificity_info["major_barcode"]))
	specificity_info=specificity_info.loc[specificity_info["specificity_corrected"]==1, :]


	# for selecting random 10 barcodes
	# specificity_info=specificity_info.loc[specificity_info["nPairs"]>1, :]
	# filtered=specificity_info.loc[specificity_info["specificity"]<0.9, :]
	# barcode_list=list(set(filtered["major_barcode"]))
	# random.shuffle(barcode_list)
	# barcode_list=barcode_list[:20]
	# blist_out=open("{}/path/{}/barcode_list.txt".format(WD, sample), "w")
	# blist_out.write('\n'.join(barcode_list))
	# specificity_info=specificity_info.loc[specificity_info["major_barcode"].isin(barcode_list),:]

	# for selecting error
	# specificity_info=specificity_info.loc[specificity_info["nPairs"]!=1, :]
	# specificity_info=specificity_info.loc[specificity_info["specificity"]<0.9, :]
	# specificity_info=specificity_info.sort_values("specificity")

	# for selecting total
	specificity_info=specificity_info.loc[specificity_info["nPairs"]!=1, :]
	specificity_info=specificity_info.loc[specificity_info["specificity"]!=1, :]
	
	specificity_info=specificity_info.reset_index()

	cluster_info=pd.read_csv("{}/cluster/{}.cluster.txt".format(WD, sample), sep='\t')
	
	pairbcd_info=pd.read_csv("{}/specificity/{}.pair-barcode.txt".format(WD, sample), sep='\t')
	pairbcd_dedup=pairbcd_info.groupby(["cluster_id","pair"])["barcode"].apply(list).reset_index()
	pairbcd_dedup["major_barcode"]=pairbcd_dedup["cluster_id"].apply(lambda x : MB_convert[x])
	pairbcd_dedup["fn_input"]=zip(pairbcd_dedup["major_barcode"], pairbcd_dedup["barcode"])
	
	pairbcd_dedup["fn_output"]=pairbcd_dedup["fn_input"].apply(lambda x : major_barcode(x[0], x[1]))
	pairbcd_dedup["major_bcd"]=pairbcd_dedup["fn_output"].apply(lambda x : x[0])
	pairbcd_dedup["major_count"]=pairbcd_dedup["fn_output"].apply(lambda x : x[1])
	pairbcd_dedup["total_count"]=pairbcd_dedup["fn_output"].apply(lambda x : x[2])
	
	color_code_df=pd.read_csv("/Data/LHS/PCR_dedup/color_reds.csv", index_col=0)
	color_code={}
	for f, c in (zip(color_code_df["a"], color_code_df["col"])):
		color_code[round(f,4)]=c
	##### process start
	sideswap={}
	sideswap["UMI_R1"]="UMI_R2"
	sideswap["UMI_R2"]="UMI_R1"
	out_df=pd.DataFrame()
	count_label=open("{}/path/{}/count.tsv".format(WD, sample), "w")
	for idx in specificity_info.index:
		cluster_id=specificity_info.loc[idx, "cluster_id"]
		fwds=specificity_info.loc[idx, "UMI_R1"]
		revs=specificity_info.loc[idx, "UMI_R2"]
		cluster_df=cluster_info.loc[(cluster_info["UMI_R1"]==fwds) & (cluster_info["UMI_R2"]==revs), :]
		cluster_idx=cluster_df.index[0]
		start_side=cluster_df.loc[cluster_idx, "Segment"]
		start_bcd=cluster_df.loc[cluster_idx, "origin"]
		
		pair_df=pd.DataFrame(cluster_df.loc[cluster_idx, "pair_list"].split(";"), columns=["pair_list"])
		pair_df["R1"]=pair_df["pair_list"].apply(lambda x: x.split("_")[0])
		pair_df["R2"]=pair_df["pair_list"].apply(lambda x: x.split("_")[1])
		ff=pair_df.groupby("R1")["R2"].apply(list).to_dict()
		rr=pair_df.groupby("R2")["R1"].apply(list).to_dict()
		pairbcd_tmp=pairbcd_dedup.loc[(pairbcd_dedup["pair"].isin(list(pair_df["pair_list"]))) & (pairbcd_dedup["cluster_id"]==cluster_id), :]
		pairbcd_tmp=pairbcd_tmp.set_index("pair")
		
		d={}
		d["UMI_R1"]=ff
		d["UMI_R2"]=rr
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
			for idx2 in range(len(route)-1):
				p1=route[idx2]
				p2=route[idx2+1]
				pair=p1+"_"+p2
				if (idx2%2==1 and start_side=="UMI_R1") or (idx2%2==0 and start_side=="UMI_R2"):
					pair=p2+"_"+p1
				pairs.append(pair)
				try:
					major_bcd=pairbcd_tmp.loc[pair, "major_bcd"]
					cnt=pairbcd_tmp.loc[pair, "major_count"]
					total=pairbcd_tmp.loc[pair, "total_count"]
					cnt=round(float(cnt)/total,4)
				except KeyError:
					major_bcd="unidentified"
					cnt=0
					total=0
				barcodes.append(major_bcd)
				counts.append(str(1-cnt))
				totals.append(str(total))

			for p, c, b in zip(pairs[1:], counts, barcodes):
				if b=="unidentified": continue
				try:
					count_label.write("{}\t{}\t{}\n".format(p, round(float(c), 4), color_code[round(float(c), 4)]))
				except KeyError: 
					count_label.write("{}\t{}\t{}\n".format(p, round(float(c), 4), "NA"))
			pairs=[sample, specificity_info.loc[idx, "major_barcode"]]+pairs
			
			barcode_list.append("-".join(barcodes))
			count_list.append(";".join(counts))
			total_list.append(";".join(totals))
			pairs_list.append('-'.join(pairs))

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


	# lineage=out_df["pairs"].apply(lambda x: x.replace("-", '\t'))
	lineage=out_df["pairs"].apply(lambda x: ( '\t'.join(x.split("-")+["NA"]*(10-len(x.split('-'))))))
	lineage.to_csv("{}/path/{}/lineage.filtered.tsv".format(WD, sample), index=False)
	out_df["ErrGT0.5_ratio"]=out_df["error_freq"].apply(lambda x: np.count_nonzero(np.array(x.split(";")).astype(float)>0.5)/float(len(x.split(";")) ))
	out_df["maxerr_streak"]=out_df["error_freq"].apply(lambda x: maximum_streak(x.split(";")) )
	out_df["maximum_errorratio"]=out_df["error_freq"].apply(lambda x: max(np.array(x.split(";")).astype(float)))
	out_df[["cluster_id",  "start_segment", "path", "barcodes","error_freq", "totals", "number_of_contents","nNodes","ErrGT0.5_ratio", "maxerr_streak"]].to_csv("{}/path/{}/path.filtered.tsv".format(WD, sample), sep='\t', index=False)