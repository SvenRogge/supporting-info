import numpy as np
import h5py
from yaff import System, Cell
import matplotlib.pyplot as pt
from molmod.units import picosecond, femtosecond, angstrom

class CP2K_xyz(object):

    def __init__(self, fn_prefix, directory='.'):
        
        from molmod.io.xyz import XYZFile

        xyz = XYZFile('{}/{}-pos-1.xyz'.format(directory, fn_prefix))
        self.geometries = xyz.geometries
        self.numbers = xyz.numbers
        self.symbols = xyz.symbols

        cell = np.loadtxt('{}/{}-1.cell'.format(directory, fn_prefix))[:, 2:-1]*angstrom
        self.cell = cell.reshape(-1, 3, 3)

sixrings = [[223, 220, 226, 217, 224, 218],
            [217, 224, 219, 223, 220, 227],
            [216, 227, 220, 222, 219, 224],
            [216, 224, 218, 222, 220, 226],
            [217, 227, 221, 223, 219, 225],
            [226, 216, 225, 218, 222, 221],
            [223, 218, 225, 217, 226, 221],
            [219, 225, 216, 227, 221, 222]]

system = System.from_file('ZIF8.chk')

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

temp = 300
xyz = CP2K_xyz('{}/ZIF8_{}K'.format(temp, temp))
pos = xyz.geometries[100:2200, :, :]
rvecs = xyz.cell[100:2200, :, :]
time = np.loadtxt('{}/ZIF8_{}K-1.cell'.format(temp, temp))[100:2200, 1]*femtosecond

print(rvecs.shape)
angles = calculate_dihedral(pos, rvecs, swing_angles, oriented = True, angle_shift=True)
angles_PBE = [61.2602,170.7308,189.2692,298.7398]
angles_closed = [61.83,171.30,188.7,298.17]
angles_open = [28.4,137.9,222.1,331.6]

np.savetxt('CP2K_angles_{}K.dat'.format(temp), angles)


colors = ['#377eb8','#4daf4a','#984ea3','#ff7f00', 'y']

pt.clf()

for idx, angle in enumerate(angles.T):
    if np.abs(np.mean(angle) - 54) < 15: color = colors[0]
    elif np.abs(np.mean(angle) - 163) < 15: color = colors[1]
    elif np.abs(np.mean(angle) - 196) < 15: color = colors[2]
    elif np.abs(np.mean(angle) - 306) < 15: color = colors[3]
    else: color = '0.5'
    pt.plot(time/picosecond, angle, color=color, alpha=0.4, label = '%s-%s-%s-%s' %(swing_angles[idx,0],swing_angles[idx,1],swing_angles[idx,2],swing_angles[idx,3]))

for angle in angles_PBE:
    pt.plot(time/picosecond, np.ones(len(time))*angle, color='k', linestyle=':', label='PBE')

for angle in angles_closed:
    pt.plot(time/picosecond, np.ones(len(time))*angle, color='#e41a1c', linestyle=':', label='closed')

for angle in angles_open:
    pt.plot(time/picosecond, np.ones(len(time))*angle, color='#a65628', linestyle=':', label='open')

pt.xlim([0, 10])
pt.xlabel('Time [ps]')
pt.ylabel('Dihedral angle [deg]')

lgd = pt.legend(bbox_to_anchor=(1.05, 1), loc=2, ncol=3)

pt.savefig('dihedrals_timeseries.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight', file='pdf')
pt.show()

pt.clf()

xaxis = np.arange(0,362,2)

angles1 = []
angles2 = []
angles3 = []
angles4 = []

max_angles = np.zeros(4)
mean_angles = np.zeros(4)
frac_open = np.zeros(4)

for angle in angles.T:
    if np.abs(np.mean(angle) - 54) < 15: angles1 += list(angle)
    elif np.abs(np.mean(angle) - 163) < 10: angles2 += list(angle)
    elif np.abs(np.mean(angle) - 196) < 10: angles3 += list(angle)
    elif np.abs(np.mean(angle) - 306) < 15: angles4 += list(angle)
    else: AssertionError('Mean angle does not correspond with known value')

def plot_hist(angles, xaxis, counter):
    hist, bin_edges = np.histogram(angles, xaxis, density = True)
    pt.plot( 0.5*(bin_edges[1:] + bin_edges[:-1]), hist, color = colors[counter])
    open_or_closed = ""
    for angle in np.array(angles):
        dist_open = np.min(np.abs(angles_open - angle))
        dist_closed = np.min(np.abs(angles_closed - angle))
        if dist_open <= dist_closed: open_or_closed += 'o'
        else: open_or_closed += 'c'
    n_open = np.sum([1 for c in open_or_closed if c == 'o'])
    frac_open[counter] = 1.0*n_open/len(angles)
    mean_angles[counter] = np.mean(angles)
    max_angles[counter] = 0.5*(bin_edges[np.argmax(hist)]+bin_edges[np.argmax(hist)+1])

plot_hist(angles1, xaxis, 0)
plot_hist(angles2, xaxis, 1)
plot_hist(angles3, xaxis, 2)
plot_hist(angles4, xaxis, 3)

pt.axvline(x=angles_PBE[0], color='k', linestyle='--', label='PBE')
for angle in angles_PBE[1:]:
    pt.axvline(x=angle, color='k', linestyle='--')

pt.axvline(x=angles_closed[0], color='#e41a1c', linestyle='--', label='closed')
for angle in angles_closed[1:]:
    pt.axvline(x=angle, color='#e41a1c', linestyle='--')

pt.axvline(x=angles_open[0], color='#a65628', linestyle='--', label='open')
for angle in angles_open[1:]:
    pt.axvline(x=angle, color='#a65628', linestyle='--')

pt.xlabel('Dihedral angle [deg]')
pt.ylabel('Histogram [1/deg]')
pt.xlim([0,360])

pt.legend(ncol=3)

pt.savefig('dihedrals_hist.pdf', bbox_inches='tight', file='pdf')
pt.show()

with open('max_hist.txt','w') as f:
    f.write('Location of the maxima in the histogram:\n')
    for value in max_angles:
        f.write('%.3f\t' %(value))
    f.write('\n\nAverage angles:\n')
    for value in mean_angles:
        f.write('%.3f\t' %(value))

with open('frac_open.txt', 'w') as f:
    f.write('Fraction of histogram closer to open structure:\n')
    for value in frac_open:
        f.write('%.2f\t' %(100*value))
