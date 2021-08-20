# Google Summer of Code 2021: Exploring Alternative Visualizations for OpenChem QCDevs
Nathan Wood, Year 4 Undergraduate, University of Florida

This project was prepared under the advisory of Dr. Farnaz Heidar-Zadeh, Dr. Paul Ayers, and Dr. Esteban Vohringer-Martinez 

## Note: Submitted Work
All work pertaining to objectives has been linked within the Results section.


## Introduction


## Goals

### Visualizing Electrostatic Potential Isosurfaces
  We wish to examine alternatives to Visualizing Molecular Dynamics (VMD) to visualize isosurfaces of molecules computed by Chemtools. Examples of VMD Alternatives include Avogadro/Avogadro2, ChimeraX, IQMol, JMol/JSMol, and PyMOL. 

### 3D Annotation of Molecules
  We wish to devise a means of taking IOData output, such as electrostatic potential charges, and applying these values as text annotations in a molecule rendered in a three dimensional environment.                        

### 2D Annotation of Molecules
  We wish to implement a means of not only visualizing a given molecule originating from a Gaussian Checkpoint file (fchk), but also annotating a 2-dimensional representation of this molecule using data from IOData.

### Plotting of Vector and Scalar Qualities of Molecules 
  We wish to develop a means of visualizing the computation of both scalar and vector qualities of molecules.

### Molecular Graphs of Bond and Ring Critical Points 
  We wish to implement a means of visualizing the bond and ring critical points of a given molecule.

## Results
### Visualizing Electrostatic Potential Isosurfaces
A script, consolidated as function `print_chimerax_isosurfaces`, has
been written, and accomplishes the following actions within the ChimeraX
Environment:

-   Open our surfaces file, generically named `*_rho.cube`

-   Open our map file, which allows us to color the surfaces,
    generically named `*_esp.cube`

-   Set to isosurface levels to visualize

-   Set the lighting, shadows, surface material appearance (ie dull or
    shiny), and representation (ie sufraces or meshes)

-   Set a color palette, both custom or built into ChimeraX

-   Set the minimum `scalemin` and maximum `scalemax` value range for
    colorizng the surfaces, or allow ChimeraX to determine the most
    suitable range using the string `'compute'` for both the `scalemin`
    and `scalemax` parameters within the script

## Post Mortem
### IOData Error in 2D Annotations 

### Discrepancies between VMD and ChimeraX 

### Molecular Graphs 
## Moving Forward
### Administrative

### Converting IOData Output into Chemical JSON 

### Molecular Graphs


## Special Thanks

## Works Cited


