import numpy as np
from yaff import System
from molmod.units import angstrom

select = 1 # 0 = structure with most negative CV, 1 = structure with least negative CV
input_chks = ['furthest_6SF6.chk','closest_6SF6.chk']
input_chk = input_chks[select]
target_dir = 'initial_structures'
if select == 0:
    CV_array = np.linspace(-11.5,-6.0,num=12)*angstrom # target array of CVs
else:
    CV_array = np.linspace(-1.5,11.5,num=27)*angstrom # target array of CVs
CV_tol = 0.25*angstrom # maximum tolerance to move CV
indices_SF6s = np.array([[2236,2237,2238,2239,2240,2241,2242],[2215,2216,2217,2218,2219,2220,2221]]) 
indices_SF6 = indices_SF6s[select,:] # indices of the SF6 molecule to be moved, starting the sulphur atom
indices_Zn = np.array([494, 1326, 1324, 778, 768, 500]) # indices of the 6 Zn atoms forming the 6MR ring, in neighboring order such that they define the normal in the correct sense


sys0 = System.from_file(input_chk)

assert (sys0.numbers[indices_Zn] == 30).all() # quick check that all indices do correspond to zinc atoms
assert sys0.numbers[indices_SF6[0]] == 16  # quick check that all indices do correspond to sulphur or fluorine atoms
assert (sys0.numbers[indices_SF6[1:]] == 9).all()  # quick check that all indices do correspond to sulphur or fluorine atoms

pos0_SF6 = sys0.pos[indices_SF6]
pos0_Zn = sys0.pos[indices_Zn]

# Calculate 6MR normals in the sense of movement of the SF6 molecules
diffs = np.array([pos0_Zn[(i+2)%6,:]-pos0_Zn[i,:] for i in range(len(indices_Zn))]) # array of r13, r24, r35, r46, r51, r62, where rij = r_j - r_i. Note: takes on form 6 x n_time x 3
normals = np.array([np.cross(diffs[i,:],diffs[(i+2)%6,:]) for i in range(len(indices_Zn))]) # construct array of normals r13 x r35, r24 x r46, etc.
normal_av = np.mean(normals, axis=0)
normal_av = normal_av/ np.linalg.norm(normal_av)

# Calculate the 6MR window centers
sixMR_center = np.mean(pos0_Zn,axis=0)

# Calculate the original projected distance between the SF6 atom and the 6MR center
CV0 = np.dot(pos0_SF6[0,:] - sixMR_center, normal_av)

print(CV0/angstrom)
print(CV_array/angstrom)

# Generate new input files along which the one SF6 molecule is moved
for CV in CV_array:
    sys_new = System.from_file(input_chk)
    sys_new.pos[indices_SF6] += (CV-CV0)*normal_av
    if CV < 0:
        sys_new.to_file('%s/init_m%d.chk' %(target_dir,10*np.abs(CV/angstrom)))
    else:
        sys_new.to_file('%s/init_%d.chk' %(target_dir,10*CV/angstrom))
