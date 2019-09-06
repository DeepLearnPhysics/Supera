from larcv import larcv
from ROOT import TChain
import numpy as np
import sys

reco_sparse3d = {}
for n in ["sparse3d_reco","sparse3d_semantics_reco"]: reco_sparse3d[n]=TChain("%s_tree" % n)
true_sparse3d = {}
for n in ["sparse3d_pcluster","sparse3d_pcluster_semantics"]: true_sparse3d[n]=TChain("%s_tree" % n)


names = [["sparse3d_reco","sparse3d_semantics_reco"],
         ["sparse3d_pcluster","sparse3d_pcluster_semantics"],
         ["sparse3d_pcluster_reco","cluster3d_pcluster_reco"],
         ["particle_pcluster","cluster3d_pcluster_reco"]]

def match_counts(names,brs):

    particle = True in ['particle' in n for n in names]
    counts = []
    for idx,br in enumerate(brs):
        n=names[idx]
        if n.startswith('sparse3d') or n.startswith('particle'): counts.append(br.as_vector().size())
        elif particle: counts.append(br.as_vector().size())
        else: counts.append(larcv.as_sparse_tensor3d(br).as_vector().size())

    if not len(np.unique(counts)) == 1:
        print('\nSize mismatch!')
        print(names)
        print(counts)
        return False
    return True
    
chains=[]    
for name_list in names:
    ch_list=[TChain('%s_tree'% n) for n in name_list]
    chains.append(ch_list)

import sys
for f in sys.argv:
    if not f.endswith('.root'): continue
    for ch_list in chains: 
        for ch in ch_list: ch.AddFile(f)

for entry in range(chains[0][0].GetEntries()):
    sys.stdout.write('%d/%d\r' % (entry,chains[0][0].GetEntries()))
    sys.stdout.flush()
    for ch_list in chains:
        for ch in ch_list: ch.GetEntry(entry)

    for idx0 in range(len(names)):
        nlist = names[idx0]
        brlist = [getattr(chains[idx0][idx1],'%s_branch' % names[idx0][idx1]) for idx1 in range(len(names[idx0]))]
        if not match_counts(nlist,brlist):
            raise Exception
