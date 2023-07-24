import numpy as np
import h5py as h5
from molmod.units import angstrom
from yaff import System 

h5_file = 'traj_1.h5'
indices_Zn = np.array([494, 1326, 1324, 778, 768, 500]) # indices of the 6 Zn atoms forming the 6MR ring, in neighboring order such that they define the normal in the correct sense
indices_S = np.array([2208,2215,2222,2229,2236,2243])
CV_array = np.linspace(-11.5,11.5,num=47)*angstrom # target array of CVs
CV_tol = 0.25*angstrom # maximum tolerance to move CV
target_dir = 'initial_structures'


with h5.File(h5_file, mode='r') as f:
    numbers = np.array(f['system/numbers'])[indices_Zn] # slicing h5 should be in increasing order; here the numpy array is sliced instead
    assert (numbers == 30).all() # quick check that all indices do correspond to zinc atoms
    pos_Zn = np.array(f['trajectory/pos'])[:,indices_Zn,:] # slicing h5 should be in increasing order; here the numpy array is sliced instead
    pos_S = np.array(f['trajectory/pos'])[:,indices_S,:] # slicing h5 should be in increasing order; here the numpy array is sliced instead

diffs = np.array([pos_Zn[:,(i+2)%6,:]-pos_Zn[:,i,:] for i in range(6)]) # array of r13, r24, r35, r46, r51, r62, where rij = r_j - r_i. Note: takes on form 6 x n_time x 3

# Calculate 6MR normals in the sense of movement of the SF6 molecules
normals = np.array([np.cross(diffs[i,:,:],diffs[(i+2)%6,:,:]) for i in range(6)]) # construct array of normals r13 x r35, r24 x r46, etc.
normals_av = np.mean(normals, axis=0)
normals_av = normals_av/ np.linalg.norm(normals_av,axis=1)[:,np.newaxis]

# Calculate the 6MR window centers
sixMR_centers = np.mean(pos_Zn,axis=1)

# Calculate the projected distances between each SF6 atom and the 6MR center and remember which molecule is the closest
CV_all = np.zeros((pos_S.shape[0], len(indices_S)))
for i in range(len(indices_S)):
    CV_all[:,i] = np.array([np.dot(pos_S[j,i,:] - sixMR_centers[j,:], normals_av[j,:]) for j in range(normals_av.shape[0])])
idx_closest = np.argmin(np.abs(CV_all), axis=1)
CV_closest = CV_all[np.arange(CV_all.shape[0]),idx_closest]

np.savetxt('closest_SF6.csv', np.vstack((indices_S[idx_closest].T,CV_closest.T/angstrom)).T, delimiter=',', header="index sulfur,CV (angstrom)")

print(CV_closest[np.argmax(np.abs(CV_closest))]/angstrom)
print(CV_closest[np.argmin(np.abs(CV_closest))]/angstrom)


with h5.File(h5_file, mode='r') as f:
    numbers=f['system/numbers'][:]
    ffatypes_old=f['system/ffatypes'][:]
    ffatypes_new = []
    for ffatype in ffatypes_old:
        if isinstance(ffatype, bytes):
            ffatypes_new.append(ffatype.decode('utf-8')) #decode byte-formatted ffatypes to strings
        else:
            ffatypes_new.append(ffatype)
    ffatype_ids=f['system/ffatype_ids'][:]
    bonds=f['system/bonds'][:]
    charges=f['system/charges'][:]
    radii=f['system/radii'][:]
    masses=f['system/masses'][:]

    for CV in CV_array:
        t_most_similar = np.argmin(np.abs(CV_closest-CV))
        CV_most_similar = CV_closest[t_most_similar]
        index_S_most_similar = indices_S[idx_closest[t_most_similar]]
        if np.abs(CV_most_similar-CV) > CV_tol:
            print('No suitable structure found for CV = %.5f A' %(CV/angstrom))
        else:
            pos = np.array(f['trajectory/pos'][t_most_similar,:,:])
            pos[np.arange(index_S_most_similar, index_S_most_similar+7),:] += (CV-CV_most_similar)*normals_av[t_most_similar,:]
            sys_new = System( numbers=numbers, pos=pos, ffatypes=ffatypes_new, ffatype_ids=ffatype_ids, bonds=bonds, rvecs=f['trajectory/cell'][t_most_similar,:,:], charges=charges, radii=radii, masses=masses)
            if CV < 0:
                sys_new.to_file('%s/init_m%d.chk' %(target_dir,10*np.abs(CV/angstrom)))
            else:
                sys_new.to_file('%s/init_%d.chk' %(target_dir,10*CV/angstrom))
