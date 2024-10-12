#!/usr/bin/env python
import os
import pandas as pd
import glob
import numpy as np
import time
from multiprocessing import Process, Queue

def GC(a):
	return (a.count("G")+a.count("C"))/float(len(a))

def network_const(WD, sample, target_info,  threads, cycle=8):
	try:
		os.mkdir("{}/cluster/".format(WD))
	except OSError:
		pass

	pair_info=pd.read_csv("{}/pair/{}_barcode_pair.txt".format(WD, sample), sep="\t", index_col=0)
	uniq=np.unique(pair_info["pair"], return_counts=True)
	count_dic=dict(zip(uniq[0], uniq[1]))
	pair_info["GC_R1"]=pair_info["frag-R1U"].apply(GC)
	pair_info["GC_R2"]=pair_info["frag-R2U"].apply(GC)
	filtered=pair_info.loc[(pair_info["GC_R1"] <0.8) & (pair_info["GC_R2"] <0.8) , :]
	targets=target_info.index
	
	if threads>len(targets):
		threads=len(targets)
	Ps=[]
	for n in range(threads):
		idxlist=[True if x%threads==n else False for x in range(len(targets))]
		targets_chunk=targets[idxlist]
		procs=Process(target=thread_function, args=(WD, sample,n,  targets_chunk, cycle, filtered, count_dic) )
		Ps.append(procs)
		procs.start()

	for p in Ps:
	    p.join()
	
	out_df=pd.DataFrame()
	
	for n in range(threads):
		while True:
			try:
				out_chunk=pd.read_csv("{}/cluster/{}.{}.cluster.txt".format(WD, sample, str(n)), sep='\t')
			except IOError:
				time.sleep(30)
				continue
			break
		out_df=out_df.append(out_chunk)
		os.system("rm {}/cluster/{}.{}.cluster.txt".format(WD, sample,str(n))) 

	out_df.to_csv("{}/cluster/{}.cluster.txt".format(WD, sample), sep='\t', index=False)


def thread_function(WD, sample, n, targets_chunk, cycle, filtered, count_dic):
	out_cluster=[]
	for target in targets_chunk:
		max_pair=(2**(cycle-2))
		## collapse pairs

		U1_df=filtered.loc[filtered["target"]==target, :].groupby("frag-R1U")["frag-R2U"].apply(set).reset_index(name="Pairs")
		U1_df["U1_size"]=U1_df["Pairs"].apply(len)
		U1_ban=set(U1_df.loc[U1_df["U1_size"]>cycle, "frag-R1U"])
		
		U2_df=filtered.loc[filtered["target"]==target, :].groupby("frag-R2U")["frag-R1U"].apply(set).reset_index(name="Pairs")
		U2_df["U2_size"]=U2_df["Pairs"].apply(len)
		U2_ban=set(U2_df.loc[U2_df["U2_size"]>cycle, "frag-R2U"])
		
		U1_df["Pairs"]=U1_df["Pairs"]-U2_ban
		U1_df=U1_df.loc[~(U1_df["frag-R1U"].isin(U1_ban)), :]
		U1_df["U1_size"]=U1_df["Pairs"].apply(len)
		U1_df["linked_UID_R2"]=U1_df["Pairs"].apply(lambda x : ";".join(list(x)))
		U1_df["target"]=target

		U2_df["Pairs"]=U2_df["Pairs"]-U1_ban
		U2_df=U2_df.loc[~(U2_df["frag-R2U"].isin(U2_ban)), :]
		U2_df["U2_size"]=U2_df["Pairs"].apply(len)
		U2_df["linked_UID_R1"]=U2_df["Pairs"].apply(lambda x : ";".join(list(x)))
		U2_df["target"]=target
		
		U1_df=U1_df.set_index("frag-R1U")
		U2_df=U2_df.set_index("frag-R2U")
		if U1_df.shape[0]==0 or U2_df.shape[0]==0: continue

		## clustering
		seedlist=list(U1_df.index)
		for seed in seedlist:
			u1_cluster=[seed]
			u2_cluster=[]
			pair_list=[]
			for u1_el in u1_cluster:
				if len(u1_cluster)>max_pair:
					break
				for u2_el in list(U1_df.loc[u1_el, "Pairs"]):
					if not u2_el in u2_cluster:
						u2_cluster.append(u2_el)
						pair=u1_el+"_"+u2_el
						pair_list.append(pair)
						for u1_el2 in list(U2_df.loc[u2_el, "Pairs"]):
							if u1_el2 not in u1_cluster:
								u1_cluster.append(u1_el2)
								pair=u1_el2+"_"+u2_el
								pair_list.append(pair)
								try:
									seedlist.remove(u1_el2)
								except ValueError:
									continue
						
			### filter cluster
			if len(u1_cluster)>max_pair or len(u2_cluster)>max_pair: continue
			if len(u1_cluster)==0 or len(u2_cluster)==0: continue
			R2s=[set(x) for x in U1_df.loc[U1_df.index.isin(u1_cluster), "Pairs"] ]
			R2s_overlap=[ len(x.intersection(y)) for x in R2s for y in R2s[R2s.index(x)+1:]]
			try:
				overlap2_max=max(R2s_overlap)
			except ValueError:
				overlap2_max=0
			R1s=[set(x) for x in U2_df.loc[U2_df.index.isin(u2_cluster), "Pairs"] ]
			R1s_overlap=[ len(x.intersection(y)) for x in R1s for y in R1s[R1s.index(x)+1:]]
			try:
				overlap1_max=max(R1s_overlap)
			except ValueError:
				overlap1_max=0
			if overlap1_max>1 or overlap2_max>1: continue

			### numbef of reads
			nreads=0
			for p in pair_list:
				nreads+=count_dic[p]

			### check origin
			u1c_max    = U1_df.loc[U1_df.index.isin(u1_cluster), "U1_size"].max()
			u1c_maxuid = U1_df.loc[U1_df.index.isin(u1_cluster), "U1_size"].idxmax()
			u2c_max    = U2_df.loc[U2_df.index.isin(u2_cluster), "U2_size"].max()
			u2c_maxuid = U2_df.loc[U2_df.index.isin(u2_cluster), "U2_size"].idxmax()
			maxi=u1c_max
			origin_uid=u1c_maxuid
			origin_side="frag-R1U"
			if u2c_max>u1c_max:
				maxi=u2c_max
				origin_uid=u2c_maxuid
				origin_side="frag-R2U"
			if maxi>cycle: continue

			### save
			outform=[";".join(u1_cluster),str(len(u1_cluster)), ";".join(u2_cluster),str(len(u2_cluster)), target, origin_side, origin_uid, str(maxi), ";".join(pair_list), str(len(pair_list)), str(nreads)]
			out_cluster.append(outform)
			
	out_cluster=pd.DataFrame(out_cluster, columns="UID_R1\tR1_size\tUID_R2\tR2_size\ttarget\tSegment\torigin\tNumber_of_branch\tpair_list\tnPairs\tnReads".split('\t'))
	out_cluster.to_csv("{}/cluster/{}.{}.cluster.txt".format(WD, sample,str(n)),sep='\t', index=False )
			