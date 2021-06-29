import numpy as np
import h5py
from yaff import System, Cell
import matplotlib.pyplot as pt
from molmod.units import picosecond, angstrom

sixrings = [[223, 220, 226, 217, 224, 218],
            [217, 224, 219, 223, 220, 227],
            [216, 227, 220, 222, 219, 224],
            [216, 224, 218, 222, 220, 226],
            [217, 227, 221, 223, 219, 225],
            [226, 216, 225, 218, 222, 221],
            [223, 218, 225, 217, 226, 221],
            [219, 225, 216, 227, 221, 222]]

system = System.from_file('init.chk')

ffatype_ids = system.ffatype_ids
ffatype_c3 = np.where(system.ffatypes=='C3')[0]
ffatype_zn = np.where(system.ffatypes=='Zn')[0]

c3_ids = np.arange(len(ffatype_ids))[ffatype_ids==ffatype_c3]
zn_ids = np.arange(len(ffatype_ids))[ffatype_ids==ffatype_zn]

c3_zns = np.zeros((len(c3_ids), 3))

swing_angles = np.zeros((len(c3_ids)*2, 4))

for idx, c3_id in enumerate(c3_ids):
    c3_zns[idx,0] = c3_id
    # Determine third neighbors of methyl carbon
    neighs3 = system.neighs3[c3_id]
    # Extract the zinc atoms for these neighbors    
    i = 1
    for neigh3 in neighs3:
        if neigh3 in zn_ids:
            c3_zns[idx,i] = neigh3
            i += 1

i = 0
for c3_zn in c3_zns:
    c3 = c3_zn[0]
    zn1 = c3_zn[1]
    zn2 = c3_zn[2]
    for sixring in sixrings:
        if zn1 in sixring and zn2 in sixring:
            i1 = np.where(np.array(sixring)==zn1)[0]
            i2 = np.where(np.array(sixring)==zn2)[0]
            i3 = int((i2 + (i2-i1))%6)
            swing_angles[i,:] = [sixring[i3], zn2, zn1, c3]
            i += 1
swing_angles = swing_angles.astype('int')

print("Found %i swing angles:" %(len(swing_angles)))
for swing_angle in swing_angles:
    print(swing_angle)

zn_neighs1 = dict((i,set([])) for i in zn_ids)

for row in c3_zns:
    zn_neighs1[int(row[1])].add(int(row[2]))
    zn_neighs1[int(row[2])].add(int(row[1]))


def calculate_dihedral(pos, cell, indices, angle_shift = False, oriented = False):
    n_angles = indices.shape[0]
    n_atoms = indices.shape[1]
    n_dim = 3
    n_time = pos.shape[0]
       
    atom_vec = np.zeros((n_time, n_angles, n_dim, 3))
    plane_vec = np.zeros((n_time, n_angles, n_dim, 2))
    angle = np.zeros((n_time, n_angles))

    for t in range(n_time):
        cell_t = Cell(cell[t,:,:])
        for i in range(n_angles):
            for j in range(3):
                delta = np.array(pos[t,indices[i,j+1], :]-pos[t,indices[i,j], :])
                cell_t.mic(delta)
                atom_vec[t,i,:,j] = delta

        for i in range(n_angles):
            # calculate the plane normals
            for j in range(2):
                plane_vec[t,i,:,j] = np.cross(atom_vec[t,i,:,j], atom_vec[t,i,:,j+1])
            angle[t,i] = np.arccos((plane_vec[t,i,:,0]*plane_vec[t,i,:,1]).sum(axis=0)/np.sqrt((plane_vec[t,i,:,0]**2).sum(axis=0)*(plane_vec[t,i,:,1]**2).sum(axis=0)))/(np.pi/180)
            if oriented:
                # determine the orientation of the cross product of both planes wrt the mutual axis
                sign = np.sign((np.cross(plane_vec[t,i,:,0], plane_vec[t,i,:,1])*atom_vec[t,i,:,1]).sum(axis=0))
                angle[t,i] *= sign
                if angle_shift and angle[t,i] < 0: angle[t,i] += 360
    return angle

fns = ['traj_%d.h5' %i for i in range(1,21)]

angles1 = []
angles2 = []
angles3 = []
angles4 = []

for fn in fns:
    print('Processing %s...' %fn)
    with h5py.File(fn, 'r') as f:
        pos = np.array(f['trajectory/pos'])
        rvecs = np.array(f['trajectory/cell'])
        time = np.array(f['trajectory/time'])

        angles = calculate_dihedral(pos, rvecs, swing_angles, oriented=True, angle_shift=True)

    for angle in angles.T:
        if np.mean(angle) > 0 and np.mean(angle) < 108: angles1 += [angle]
        elif np.mean(angle) < 179: angles2 += [180-angle]
        elif np.mean(angle) < 251: angles3 += [-(180-angle)]
        else: angles4 += [angle]

np.savetxt('dihedrals_1.txt', np.array(angles1).T)
np.savetxt('dihedrals_2.txt', np.array(angles2).T)
np.savetxt('dihedrals_3.txt', np.array(angles3).T)
np.savetxt('dihedrals_4.txt', np.array(angles4).T)
