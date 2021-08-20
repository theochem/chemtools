def plotVectorFieldOverPlane(inFile,step_size=0.32, title='Plot of Vector Field'):
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
    
    coord1 : np.ndarray()
        First point on plane 
        
    coord2 : np.ndarray() 
        Second Point on plane
    
    coord3 : np.ndarray() 
        Third Point on plane
    
    normal : np.cross()
        orthonormalization of coordinates
        
    vectorFunc : callabel() 
        Vector Field that inputs 3D point to 3D point 
        In this example, vectorFunc is mol.compute_gradient()
    
    
    stepSize : float
        Step Size of 2D Affine Grid
        
    title : str
        Title of plot
    
    color : str
        Color of Gradient Arrows  
    
    scale : int
        Scale of Gradient Arrows
    """
    # DEPENDENCIES
    import numpy as np

    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d

    from chemtools import Molecule
    #END DEPENDENCIES
    
    mol = Molecule.from_file(inFile)
    coord = mol.coordinates
    coord1 = coord[0] 
    coord2 = coord[1]
    coord3 = coord[2]
    
    v1, v2 = coord2 - coord1, coord3 - coord1
    normal = np.cross(v1, v2)
    normal /= np.linalg.norm(normal)

    a, b, c = normal
    
      
    a,b,c = normal
    cos_theta = c
    sin_theta = np.sqrt(a**2.0 + b**2.0)
    u_1 = b / np.sqrt(a**2.0 + b**2.0)
    u_2 = -a / np.sqrt(a**2.0 + b**2.0)
    rot = np.array([
            
            [cos_theta + u_1**2.0 * (1 - cos_theta), u_1 * u_2 * (1 - cos_theta), u_2 * sin_theta],
            [u_1 * u_2 * (1 - cos_theta), cos_theta + u_2**2.0 * (1 - cos_theta), -u_1 * sin_theta],
            [-u_2 * sin_theta, u_1 * sin_theta, cos_theta]
            
        ])
    rot_coords = np.dot(rot, (coord - coord[0]).T).T


    l_bnd = np.min(rot_coords, axis=0) - 1
    u_bnd = np.max(rot_coords, axis=0) + 1
    print(l_bnd)
    
    x_grid = np.arange(l_bnd[0], u_bnd[0] + step_size, step_size)
    y_grid = np.arange(l_bnd[1], u_bnd[1] + step_size, step_size)
    grid_2d = np.array(np.meshgrid(x_grid, y_grid)).T.reshape(-1,2)

    # Add zero z-axis and rotate it and translate it to the plane
    grid_zeros = np.hstack((grid_2d, np.zeros((grid_2d.shape[0], 1), dtype=np.float)))
    grid_plane = np.einsum("ij,kj->ki", rot.T, grid_zeros)
    grid_plane += coord[0]
    
    vectorFunc = mol.compute_gradient(grid_plane)
    gradients = vectorFunc
    print(gradients)

    proj_gradients = gradients - np.dot(gradients, normal)[:, np.newaxis] * normal
    rot_proj_gradients = np.dot(rot, proj_gradients.T).T
    
    # Rotate Projected Plane.
    rot_proj_gradients /= np.linalg.norm(rot_proj_gradients, axis=1).reshape((-1, 1))
    plt.figure(figsize=(10,10))
    plt.quiver(grid_2d[:, 0], grid_2d[:, 1],
         rot_proj_gradients[:, 0], rot_proj_gradients[:, 1],color='black', scale=45)
    plt.plot(rot_coords[:, 0], rot_coords[:, 1], "ro", label="Molecular Coordinates")
    plt.legend()
    
    plt.title(title)
    
    plt.show() 

    #END FUNCTION

inFile = 'dichloropyridine26_q+0.fchk'

plotVectorFieldOverPlane(inFile,step_size=0.32, title='Plot of Vector Field') 
