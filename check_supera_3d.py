import numpy as np
from larcv import larcv

def xyz_from_tensor3d(vs,meta):
    """
    Given a larcv::VoxelSet (vs) and larcv::Voxel3DMeta, constructs 3D index array of shape (N,3)
    """
    x = np.zeros(shape=(vs.size()),dtype=np.int32)
    y = np.zeros(shape=(vs.size()),dtype=np.int32)
    z = np.zeros(shape=(vs.size()),dtype=np.int32)
    v = np.zeros(shape=(vs.size()),dtype=np.float32)
    larcv.as_flat_arrays(vs,meta,x,y,z,v)
    return np.stack([x,y,z],axis=1)
    
def same_tensor3d(t0,t1,meta0=None,meta1=None,explicit=False):
    """
    Given two set of "tensors", compare to see if they are "identical"
    INPUT:
      - t0 ... either EventSparse3DTensor (then meta0 should be None) or VoxelSet (then meta0 needed)
      - t1 ... either EventSparse3DTensor (then meta1 should be None) or VoxelSet (then meta1 needed)
      - meta0 ... Voxel3DMeta, if t0 is VoxelSet which does not carry meta info
      - meta1 ... Voxel3DMeta, if t1 is VoxelSet which does not carry meta info
      - explicit ... if True, ALL elements are compared. If false, only start/end+randomly chosen elements compared
    """
    if not t0.size() == t1.size(): return False
    if t0.size() == 0: return True
    if meta0 is None: meta0 = t0.meta()
    if meta1 is None: meta1 = t1.meta()
    if explicit:
        a0 = xyz_from_tensor3d(t0,meta0)
        a1 = xyz_from_tensor3d(t1,meta1)
        return (a0 == a1).astype(np.int32).sum() == np.prod(a0.shape)
    else:
        # only check if voxel ID agrees at the start, end, and random min(10,t0.size()) elements
        elements = [0,1] + list(np.random.randint(0,t0.size()-1,size=(min(10,t0.size()))))
        elements = np.unique(elements).astype(np.int32)
        v0 = t0.as_vector()
        v1 = t1.as_vector()
        for e in elements:
            if not v0[int(e)].id() == v1[int(e)].id(): return False
    return True

def subset_tensor3d(t0,t1,meta0=None,meta1=None):
    """
    Given two set of "tensors", compare to see if one (t1) is a subset of the other (t0).
    INPUT:
      - t0 ... either EventSparse3DTensor (then meta0 should be None) or VoxelSet (then meta0 needed)
      - t1 ... either EventSparse3DTensor (then meta1 should be None) or VoxelSet (then meta1 needed)
      - meta0 ... Voxel3DMeta, if t0 is VoxelSet which does not carry meta info
      - meta1 ... Voxel3DMeta, if t1 is VoxelSet which does not carry meta info
    """
    if not t0.size() >= t1.size(): return False
    if t0.size() == 0: return True
    if t1.size() == 0: return True
    if meta0 is None: meta0 = t0.meta()
    if meta1 is None: meta1 = t1.meta()
    a0 = xyz_from_tensor3d(t0,meta0)
    a1 = xyz_from_tensor3d(t1,meta1)
    subset_not = np.array([not pt in a0 for pt in a1]).astype(np.int32).sum()
    return (subset_not == 0)

def same_cluster3d(c0,c1,explicit=False):
    """
    Given two larcv::EventCluserVoxel3D, compare if they are identical in voxel count and locations
    for each cluster element.
    """
    if not c0.size() == c1.size(): return False
    if c0.size() == 0: return True
    for i in range(c0.size()):
        if not same_tensor3d(c0.as_vector()[i],c1.as_vector()[i],c0.meta(),c1.meta(),explicit):
            return False
    return True

def subset_cluster3d(c0,c1):
    """
    Given two larcv::EventCluserVoxel3D, compare if one (c1) is a subset of the other (c0).
    Here, by "subset", it means individual clusters. So the length of c0 and c1 must match.
    """
    if not c0.size() == c1.size(): return False
    if c0.size() == 0: return True
    for i in range(c0.size()):
        if not subset_tensor3d(c0.as_vector()[i],c1.as_vector()[i],c0.meta(),c1.meta()):
            return False
    return True

class data_blob:
    """
    Utility class to handle arbitrary many larcv data products
    """
    def __init__(self,names):
        from ROOT import TChain
        self._names = names
        self._chains = {}
        for name in self._names: self._chains[name] = TChain('%s_tree' % name)
        
    def add_file(self,f):
        for name,ch in self._chains.items(): ch.AddFile(f)
    def get_entry(self,e):
        for name,ch in self._chains.items(): ch.GetEntry(e)
    def get_data(self,name):
        return getattr(self._chains[name],"%s_branch" % name)
    def get_entries(self):
        entry = self._chains[self._names[0]].GetEntries()
        for name,ch in self._chains.items():
            if not ch.GetEntries() == entry:
                print('Chain',name,'has different entry',ch.GetEntries(),'!=',entry)
                raise ValueError
        return entry

def check_supera():
    import sys
    names = ['sparse3d_reco',
             'sparse3d_pcluster',
             'sparse3d_pcluster_semantics',
             'sparse3d_pcluster_semantics_ghost',
             'cluster3d_pcluster',
             'cluster3d_pcluster_highE',
             'cluster3d_pcluster_lowE',
             'particle_pcluster',
             'particle_corrected']
    blob = data_blob(names)
    fs = [v for v in sys.argv if v.endswith('.root')]
    for f in fs: blob.add_file(f)

    entries = blob.get_entries()
    for entry in range(entries):
        sys.stdout.write('%d/%d\r' % (entry,entries))
        sys.stdout.flush()
        blob.get_entry(entry)
        error=False
        # sparse3d_reco and sparse3d_semantics_reco
        s1, s2 = 'sparse3d_reco', 'sparse3d_pcluster_semantics_ghost'
        if not same_tensor3d(blob.get_data(s1),blob.get_data(s2)):
            print('\nNot same tensor3d',s1,'(%d)' % blob.get_data(s1).size(),s2,'(%d)' % blob.get_data(s2).size())
            error=True
        # sparse3d_pcluster and sparse3d_pcluster_semantics
        s1, s2 = 'sparse3d_pcluster', 'sparse3d_pcluster_semantics'
        if not same_tensor3d(blob.get_data(s1),blob.get_data(s2)):
            print('\nNot same tensor3d',s1,'(%d)' % blob.get_data(s1).size(),s2,'(%d)' % blob.get_data(s2).size())
            error=True
        # cluster3d_pcluster and cluster3d_pcluster_highE
        c1, c2 = 'cluster3d_pcluster', 'cluster3d_pcluster_highE'
        if not subset_cluster3d(blob.get_data(c1),blob.get_data(c2)):
            print('\nNot a subset cluster3d',c1,c2);
            error=True
        # cluster3d_pcluster and cluster3d_pcluster_lowE
        c1, c2 = 'cluster3d_pcluster', 'cluster3d_pcluster_lowE'
        if not subset_cluster3d(blob.get_data(c1),blob.get_data(c2)):
            print('\nNot a subset cluster3d',c1,c2);
            error=True
        # cluster3d_pcluster_reco and cluster3d_pcluster_highE_reco
        c1, c2 = 'cluster3d_pcluster', 'cluster3d_pcluster_highE'
        if not subset_cluster3d(blob.get_data(c1),blob.get_data(c2)):
            print('\nNot a subset cluster3d',c1,c2);
            error=True
        # make sure the size is same
        p='particle_pcluster'
        part_count = blob.get_data(p).as_vector().size()
        for n in names:
            if not n.startswith('cluster3d_'): continue
            if not blob.get_data(n).as_vector().size() == part_count:
                print('\nParticle count mismatch',p,n)
                error=True
        if error:
            d = blob.get_data(p)
            print('Bad entry',entry,'run',d.run(),'subrun',d.subrun(),'event',d.event())


if __name__ == '__main__':
    check_supera()
    
