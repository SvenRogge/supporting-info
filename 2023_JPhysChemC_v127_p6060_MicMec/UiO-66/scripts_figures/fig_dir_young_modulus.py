#!/usr/bin/env python
# File name: fig_dir_young_modulus_3D_3x3x3_confs.py
# Author: Joachim Vandewalle
# Date: 01-10-2022

"""Create a 3D directional Young modulus plot."""

import json

import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt


with open("output_configurations.json") as jfile:
    jdict = json.load(jfile)
    

def _young_modulus_from_compliance_tensor(compliance_tensor):
    # Create the mesh in spherical coordinates and compute corresponding E.
    theta = np.linspace(0, np.pi, 100)
    phi = np.linspace(0, 2*np.pi, 100)
    PHI, THETA = np.meshgrid(phi, theta)
    U = [np.cos(PHI)*np.sin(THETA), np.sin(PHI)*np.sin(THETA), np.cos(THETA)]

    E = np.zeros(THETA.shape)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    E += U[i]*U[j]*U[k]*U[l]*compliance_tensor[i,j,k,l]
    return PHI, THETA, np.absolute(1/E)
    

def plot_dir_young(compliance_tensor, conf):

    PHI, THETA, E = _young_modulus_from_compliance_tensor(compliance_tensor)

    # Express the mesh in the cartesian system.
    X, Y, Z = E*np.cos(PHI)*np.sin(THETA), E*np.sin(PHI)*np.sin(THETA), E*np.cos(THETA)


    def set_axes_equal(ax: plt.Axes):
        """Set three-dimensional plot axes to equal scale.

        Make the axes of a three-dimensional plot have equal scale so that spheres appear as spheres and cubes as cubes.  
        Required since ``ax.axis("equal")`` and ``ax.set_aspect("equal")`` don't work on three-dimensional plots.
        """
        limits = np.array([
            ax.get_xlim3d(),
            ax.get_ylim3d(),
            ax.get_zlim3d(),
        ])
        origin = np.mean(limits, axis=1)
        radius = 0.5*np.max(np.abs(limits[:,1] - limits[:,0]))
        _set_axes_radius(ax, origin, radius)

    def _set_axes_radius(ax, origin, radius):
        x, y, z = origin
        # Alter the limits manually if these automatic limits are not to your liking.
        ax.set_xlim3d([-30, 30])
        ax.set_ylim3d([-30, 30])
        ax.set_zlim3d([-30, 30])


    # Plot the surface.
    fig = plt.figure()

    ax = fig.add_subplot(projection="3d")

    ax.set_box_aspect([1,1,1])
    ax.plot_surface(X, Y, Z)
    set_axes_equal(ax)
    ax.set_xlabel(r"$\mathrm{E_x \; [GPa]}$")
    ax.set_ylabel(r"$\mathrm{E_y \; [GPa]}$")
    ax.set_zlabel(r"$\mathrm{E_z \; [GPa]}$")

    ax = plt.gca()
    # Delete the tick labels of one or more axes.
    #ax.xaxis.set_ticklabels([])
    #ax.yaxis.set_ticklabels([])
    #ax.zaxis.set_ticklabels([])

    # Reduce the number of tick labels for better visibility.
    every_nth = 2
    for n, label in enumerate(ax.xaxis.get_ticklabels()):
        if n % every_nth != 0:
            label.set_visible(False)
    for n, label in enumerate(ax.yaxis.get_ticklabels()):
        if n % every_nth != 0:
            label.set_visible(False)
    for n, label in enumerate(ax.zaxis.get_ticklabels()):
        if n % every_nth != 0:
            label.set_visible(False)

    # Delete the tick lines of one or more axes.
    #for line in ax.xaxis.get_ticklines():
    #    line.set_visible(False)
    #for line in ax.yaxis.get_ticklines():
    #    line.set_visible(False)
    #for line in ax.zaxis.get_ticklines():
    #    line.set_visible(False)

    plt.savefig(f"fig_dir_young_modulus_3x3x3_conf{conf}.png")
    plt.show()


V = {
    0: (0,0),
    1: (1,1), 
    2: (2,2), 
    3: (1,2),
    4: (0,2),
    5: (0,1)
}


def voigt_inv(matrix, mode=None):
    """Map a (6 x 6) Voigt notation matrix to a (3 x 3 x 3 x 3) tensor.
    
    Parameters
    ----------
    matrix : numpy.ndarray, shape=(6, 6)
        The Voigt notation matrix to be mapped to a tensor.
    mode : {"compliance", "elasticity"}, optional
        Declare whether the input matrix is a compliance matrix or an elasticity matrix.

    Returns
    -------
    tensor : numpy.ndarray, shape=(3, 3, 3, 3)
        The resulting tensor.

    Notes
    -----
    Voigt notation differs depending on whether the tensor is a compliance tensor or an elasticity tensor,
    hence the (optional) keyword ``mode``.   
    """
    if matrix.shape == (6, 6):
        tensor = np.zeros((3,3,3,3))
        if (mode is None) or (mode == "compliance"):
            for index, _ in np.ndenumerate(tensor):
                ij = tuple(sorted(index[0:2]))
                kl = tuple(sorted(index[2:4]))
                for key in V.keys():
                    if V[key] == ij:
                        V_ij = key
                    if V[key] == kl:
                        V_kl = key
                tensor[index] = matrix[(V_ij, V_kl)]
                if V_ij >= 3:
                    tensor[index] *= 0.5
                if V_kl >= 3:
                    tensor[index] *= 0.5
        elif (mode == "elasticity"):
            for index, _ in np.ndenumerate(tensor):
                ij = tuple(sorted(index[0:2]))
                kl = tuple(sorted(index[2:4]))
                for key in V.keys():
                    if V[key] == ij:
                        V_ij = key
                    if V[key] == kl:
                        V_kl = key          
                tensor[index] = matrix[(V_ij, V_kl)]
        else:
            raise ValueError("Method `voigt_inv` did not receive valid input for keyword `mode`.")
        return tensor
    elif matrix.shape == (6,):
        tensor = np.zeros((3,3))
        if (mode is None) or (mode == "strain"):
            for idx, _ in np.ndenumerate(tensor):
                ij = tuple(sorted(idx))
                for key in V.keys():
                    if V[key] == ij:
                        V_ij = key
                tensor[idx] = matrix[V_ij]
                if V_ij >= 3:
                    tensor[idx] *= 0.5
        elif (mode == "stress"):
            for idx, _ in np.ndenumerate(tensor):
                ij = tuple(sorted(idx))
                for key in V.keys():
                    if V[key] == ij:
                        V_ij = key
                tensor[idx] = matrix[V_ij]
        else:
            raise ValueError("Method `voigt_inv` did not receive valid input for keyword `mode`.")
        return tensor
    else:
        raise ValueError("Method `voigt_inv` did not receive valid input for `matrix`.")


elasticity_matrix = np.array(jdict["C_atomic_conf0"])
plot_dir_young(voigt_inv(np.linalg.inv(elasticity_matrix), mode="compliance"), 6)


