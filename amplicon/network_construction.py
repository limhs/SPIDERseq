import os
import pandas as pd
import glob
import numpy as np
import time
from multiprocessing import Process, Queue

def GC(a):
	return (a.count("G")+a.count("C"))/float(len(a))

def network_const(WD, sample, target_info,  n_core, cycle=8):
	try:
		os.mkdir("{}/cluster/".format(WD))
	except OSError:
		pass

	pair_info=pd.read_csv("{}/pair/{}_barcode_pair.txt".format(WD, sample), sep="\t", index_col=0)
	uniq=np.unique(pair_info["pair"], return_counts=True)
	count_dic=dict(zip(uniq[0], uniq[1]))
	pair_info["GC_R1"]=pair_info["R1_UMI"].apply(GC)
	pair_info["GC_R2"]=pair_info["R2_UMI"].apply(GC)
	filtered=pair_info.loc[(pair_info["GC_R1"] <0.8) & (pair_info["GC_R2"] <0.8) , :]
	# filtered=pair_info.loc[(pair_info["GC_R1"] <0.8) & (pair_info["GC_R1"] >0.2) & (pair_info["GC_R2"] <0.8) & (pair_info["GC_R2"] >0.2), :]

	targets=target_info.index
	if n_core>len(targets):
		n_core=len(targets)
	Ps=[]
	for n in range(n_core):
		idxlist=[True if x%n_core==n else False for x in range(len(targets))]
		targets_chunk=targets[idxlist]
		procs=Process(target=thread_function, args=(WD, sample,n,  targets_chunk, cycle, filtered, count_dic) )
		Ps.append(procs)
		procs.start()

	for p in Ps:
	    p.join()
	
	out_df=pd.DataFrame()
	
	for n in range(n_core):
		while True:
			try:
				out_chunk=pd.read_csv("{}/cluster/{}.{}.cluster.txt".format(WD,  sample,str(n)), sep='\t')
			except IOError:
				time.sleep(30)
				continue
			break
		out_df=out_df.append(out_chunk)
		os.system("rm {}/cluster/{}.{}.cluster.txt".format(WD, sample,str(n))) 


	out_df.to_csv("{}/cluster/{}.cluster.txt".format(WD, sample), sep='\t', index=False)

	

def thread_function(WD, sample,n,  targets_chunk, cycle, filtered, count_dic):
	out_cluster=[]
	for target in targets_chunk:
		

		max_pair=(2**(cycle-2))
		## collapse pairs
		U1_df=filtered.loc[filtered["target"]==target, :].groupby("R1_UMI")["R2_UMI"].apply(set).reset_index(name="Pairs")
		U1_df["size"]=U1_df["Pairs"].apply(len)
		U1_ban=set(U1_df.loc[U1_df["size"]>cycle, "R1_UMI"])
		
		U2_df=filtered.loc[filtered["target"]==target, :].groupby("R2_UMI")["R1_UMI"].apply(set).reset_index(name="Pairs")
		U2_df["size"]=U2_df["Pairs"].apply(len)
		U2_ban=set(U2_df.loc[U2_df["size"]>cycle, "R2_UMI"])
		
		U1_df["Pairs"]=U1_df["Pairs"]-U2_ban
		U1_df=U1_df.loc[~(U1_df["R1_UMI"].isin(U1_ban)), :]
		U1_df["size"]=U1_df["Pairs"].apply(len)
		U1_df["linked_UMI_R2s"]=U1_df["Pairs"].apply(lambda x : ";".join(list(x)))
		U1_df["target"]=target

		U2_df["Pairs"]=U2_df["Pairs"]-U1_ban
		U2_df=U2_df.loc[~(U2_df["R2_UMI"].isin(U2_ban)), :]
		U2_df["size"]=U2_df["Pairs"].apply(len)
		U2_df["linked_UMI_R1s"]=U2_df["Pairs"].apply(lambda x : ";".join(list(x)))
		U2_df["target"]=target
		
		U1_df=U1_df.set_index("R1_UMI")
		U2_df=U2_df.set_index("R2_UMI")
		if U1_df.shape[0]==0 or U2_df.shape[0]==0: continue
		# U1_df.loc[:, U1_df.columns.isin(["R1_UMI", "linked_UMI_R2s", "size", "target"])].to_csv("{}/cluster/{}.{}.u1.tsv".format(WD,target, sample),  sep='\t')
		# U2_df.loc[:, U2_df.columns.isin(["R2_UMI", "linked_UMI_R1s", "size", "target"])].to_csv("{}/cluster/{}.{}.u2.tsv".format(WD,target, sample),  sep='\t')

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
						
			
			### cluster_size filter
			if len(u1_cluster)>max_pair or len(u2_cluster)>max_pair: continue
			if len(u1_cluster)==0 or len(u2_cluster)==0: continue
			### multibridge check
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

			### origin(which has largest number of branch) check
			u1c_max    = U1_df.loc[U1_df.index.isin(u1_cluster), "size"].max()
			u1c_maxumi = U1_df.loc[U1_df.index.isin(u1_cluster), "size"].idxmax()
			u2c_max    = U2_df.loc[U2_df.index.isin(u2_cluster), "size"].max()
			u2c_maxumi = U2_df.loc[U2_df.index.isin(u2_cluster), "size"].idxmax()
			maxi=u1c_max
			origin_umi=u1c_maxumi
			origin_side="UMI_R1"
			if u2c_max>u1c_max:
				maxi=u2c_max
				origin_umi=u2c_maxumi
				origin_side="UMI_R2"
			if maxi>cycle: continue

			### save
			outform=[";".join(u1_cluster),str(len(u1_cluster)), ";".join(u2_cluster),str(len(u2_cluster)), target, origin_side, origin_umi, str(maxi), ";".join(pair_list), str(len(pair_list)), str(nreads)]
			out_cluster.append(outform)
			
	
	out_cluster=pd.DataFrame(out_cluster, columns="UMI_R1\tR1_size\tUMI_R2\tR2_size\ttarget\tSegment\torigin\tNumber_of_branch\tpair_list\tnPairs\tnReads".split('\t'))
	out_cluster.to_csv("{}/cluster/{}.{}.cluster.txt".format(WD, sample,str(n)),sep='\t', index=False )
			