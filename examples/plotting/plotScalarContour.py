import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from chemtools import Molecule 

def plotScalarFunctionContourPlot(inFile, step_size=0.3, title= 'Contour Plot'):
    """
    This script will plot a vector quality of a molecule as a gradient plot
    
    CONDA ENVIRONMENT
        This function must be executed in a conda environment containing: 
            chemtools
            
                link to configuring chemtools can be found here
                
                https://chemtools.org/usr_doc_installation.html
                
            Numpy
            MatPlotLib 
    PARAMETERS
    
        inFile, str
            This represents our input Gaussian Checkpoint *.fchk file
        
        mol, chemtools object
            load inFile as chemtools Molecule module
                
        coord, chemtools API Call
            cartesian coordinates of inFile
            
        rot, np.array()
            Rotation Matrix used to rotate affine plane and normalize for plotting of scalarFunc
        
    """
    # Obtain Cartesian Coordinates
    mol = Molecule.from_file(inFile)
    coord = mol.coordinates
    coord1 = coord[0] 
    coord2 = coord[1]
    coord3 = coord[2]

    # Orthonormalization
    v1, v2 = coord1 - coord3, coord2 - coord3
    normal = np.cross(v1, v2)
    normal /= np.linalg.norm(normal)
    print("Orthonormal : ")
    a, b, c = normal
  
    # This represents our Rotate and Scale Operation,
    cos_theta = c
    sin_theta = np.sqrt(a**2.0 + b**2.0)
    u_1 = b / np.sqrt(a**2.0 + b**2.0)
    u_2 = -a / np.sqrt(a**2.0 + b**2.0)
    rot = np.array([
        [cos_theta + u_1**2.0 * (1 - cos_theta), u_1 * u_2 * (1 - cos_theta), u_2 * sin_theta],
        [u_1 * u_2 * (1 - cos_theta), cos_theta + u_2**2.0 * (1 - cos_theta), -u_1 * sin_theta],
        [-u_2 * sin_theta, u_1 * sin_theta, cos_theta]
    ])
    
    
    
    print("Rotation Operation: ")
    print(rot)

    rot_coords = np.dot(rot, (coord - coord[0]).T).T
    print("Rotation Coordinates: ")
    print(rot_coords)
    l_bnd = np.min(rot_coords, axis=0) - 1
    u_bnd = np.max(rot_coords, axis=0) + 1

    x_grid = np.arange(l_bnd[0], u_bnd[0] + step_size, step_size)
    y_grid = np.arange(l_bnd[1], u_bnd[1] + step_size, step_size)
    grid_2d = np.array(np.meshgrid(x_grid, y_grid)).T.reshape(-1,2)
    
    # fill grid space with zeros to later populate with scalarFunc
    grid_zeros = np.hstack((grid_2d, np.zeros((grid_2d.shape[0], 1), dtype=np.float)))
    grid_plane = np.einsum("ij,kj->ki", rot.T, grid_zeros)
    grid_plane += coord[0]
    print("Grid Plane:")
    print(grid_plane)

    # setup x,y for contour plot 
    x = grid_2d[:,0]
    y = grid_2d[:,1]
    xy = np.meshgrid(grid_2d[:,0],grid_2d[:,1])
    
    scalarFunc = mol.compute_density(grid_plane)

    # establish number of level curves
    levels = np.array([0.001 * n * n for n in range(1000)])
    
    # reshape level curves array to resemble length of x and y , with Fortran style iteration
    scalarFuncPlot = scalarFunc.reshape((len(x_grid), len(y_grid)), order="F")
    plt.figure(figsize=(10,10))
    plt.contour(x_grid,y_grid,scalarFuncPlot, levels )
    plt.title("Contour Plot of Electrostatic Density")

    plt.show()
    #END FUNCTION

# Example using 2,6-Dichloropyridine
inFile = 'dichloropyridine26_q+0.fchk'

plotScalarFunctionContourPlot(inFile, step_size=0.3, title= 'Contour Plot')  
