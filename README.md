\maketitle
\vspace{1cm}
This report has been prepared for Google Summer of Code 2021, under the
advisory of Dr.Paul Ayers, Dr. Farnaz Heidar-Zadeh, Dr. Esteban
Vohringer-Martinez\

![image](openchem.png){width="4cm"}

\vspace{2cm}
**Abstract**\
OpenChem wishes to implement a means of visualizing Chemtools data. Such
means include exploring viable alternatives to Visualizing Molecular
Dynamics as a visualization platfrom. Further visualization goals
include the plotting of gradient vector fields and annotation of
molecules. Additionally, ChimeraX, developed by the University of
California, San Francisco, was chosen as a general purpose alternative
for VMD, with regards to visualizing isosurfaces as well as molecular
annotations.

Introduction
============

Visualization is important to make sense of computational chemistry
data. Visualization methods relevant to OpenChem include the
visualization of isosurfaces resulting from electrostatic potential
calculations, the annotations of molecules, as well as plotting the
physical and chemical qualities of specific molecules. In order to
accomplish this, the gap between the computational tools and currently
existing visualization tools must be bridged.

Project
=======

Visualizing Electrostatic Potential Isosurfaces in ChimeraX
-----------------------------------------------------------

We wish to examine alternatives to Visualizing Molecular Dynamics (VMD)
to visualize isosurfaces of molecules computed by Chemtools. Examples of
VMD Alternatives include Avogadro/Avogadro2, ChimeraX, IQMol,
JMol/JSMol, and PyMOL.

3D Annotation of Molecules Using IOData Output
----------------------------------------------

We wish to devise a means of taking IOData output, such as electrostatic
potential charges, and applying these values as text annotations in a
molecule rendered in a three dimensional environment.

2D Annotation of Molecules Using IOData Output
----------------------------------------------

' We wish to implement a means of not only visualizing a given molecule
originating from a Gaussian Checkpoint file (fchk), but also annotating
a 2-dimensional representation of this molecule using data from IOData.

Gradient Vector Fields and Scalar Contour Plots
-----------------------------------------------

We wish to develop a means of visualizing the computation of both scalar
and vector qualities of molecules.

Molecular Graphs
----------------

We wish to implement a means of visualizing the bond and ring critical
points of a given molecule.

Results
=======

Visualizing Electrostatic Potential Isosurfaces
-----------------------------------------------

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

3D Annotation of Molecules
--------------------------

IOData can output to Protein Data Bank (PDB) files, a text file that
enumerates monomer and a atomic data, such as element, XYZ cartesian
coordinates, and B Factor, which has historically been used to quantify
flexibility and freedom of motion in computational protein-ligand
docking. Originally, ChimeraX's Attribute Files, which allow a tab
separated table enumerating characteristics per unit, was conceived as a
means of annotating molecules as it was not considered to be a \"hack\"
and could allow multiple different annotations could be performed
without overwriting our PDB. However, this was abandoned due to the
complexity of assigning IOData values to an format ChimeraX finds
accessible.\
IOData possesses a feature known as `extra`, which allows writing in B
Factors when an intermediate PDB file is opened using IOData function
`load_one`, which can then be written to a finalized PDB file to be
loaded into ChimeraX. Once loaded, ChimeraX can visualize the B Factor
using its labeling utilities.

![3D Annotation of 2,6-Dichloropyridine With Electrostatic Potential
Charges](3dannotate.png){width="3in"}

2D Annotation of Molecules
--------------------------

IOData is used to dump an XYZ and PDB files from the original Gaussian
Checkpoint (fchk) file. XYZ files are capable of being converted to
accurate SMILES strings via the Chemml utilities, however, the order of
the atoms being enumerated differs from that of IOData, making it
non-trivial to assign charges resulting from IOData's `atcharges` call.
On the other hand, PDB Files maintain the atomic order as IOData,
however, RDKit cannot visualize any pi bonds from the PDB. This is due
to the fact that PDB Files store pi bond data based on the distance
between atoms in cartesian space, whereas RDKit prefers bonds to be
explicitly detailed.

To address both issues, we can first correct the geometry of the PDB by
using RDKit's `Compute2DCoords` and then force the PDB to assume the
visualization stile of the afforementioned SMILES string using
`AssignBondOrdersFromTemplate`. Afterwards, the charges from IOData can
be assigned iteratively.

![2D Annotation of 2,6-Dichloropyridine, Using Electrostatic Potential
Charges](test11.png){width="3in"}

Plotting
--------

### Gradient Vector Fields

A Jupyter Notebook and function `plotVectorFieldOverPlane` were made
that is capable of plotting gradient vector fields. This was performed
by:

-   Using 3D rotation matrix function `_rotate_3d()` to rotate an affine
    plane and normalize to plot gradients against

-   Assigning our `vectorFunc` as as the vector quality we want to
    explore, in this case `mol.compute_gradient()`

-   Projecting our vector quality against the rotated affine plane

-   Visualizing the vector quality as a gradient vector field \"Quiver
    Plot\"

![Gradient Vector Field of Electron Density,
2,6-Dichloropyridine](dichlorgradient.png){width="4in"}

To rotate a scalar quality, in this case electron density, we want to

-   Rotate an affine plane and normalize to plot gradients against using
    the 2D rotate and scale operation, using the helper function
    `_rotate_2d()`

-   Assign our `scalarFunc` which is a scalar quality of a molecule, in
    this case, electrostatic density

-   Visualizing the scalar quality as a contour plot, with our
    `scalarFunc` representing the level curves

![Contour Plot of Electrostatic Densities,
2,6-Dichloropyridine](contour.png){width="4in"}

Post Mortem
===========

IOData Errors in 2D Annotations
-------------------------------

When working with other molecules, such as Acrolein, IOData will throw
the following error

                Exception ignored in: <function LineIterator.__del__ at 0x7fd2980ec440>
                Traceback (most recent call last):
                  File "/home/woodn/.conda/envs/chemml_env/lib/python3.7/site-packages/iodata/utils.py", line 80, in __del__
                    self.f.close()
                AttributeError: 'LineIterator' object has no attribute 'f'
            

Despite this error, IOData still dumps the necessary XYZ and PDB Files,
and the RDKit visualization is complete.

![Acrolein, Annotated Using IOData Derived Electrostatic Potential
Charges](acroleinerror.png){width="3in"}

Discrepancies between ChimeraX and VMD
--------------------------------------

### Isosurface Levels

There is a discrepancy between appropriate values to display isosurface
levels between both ChimeraX and VMD. For the example
2,6-dichloropyridine, the previous VMD example recommends setting
`isoSurf` to .003, however, in ChimeraX, this will produce a noticeable
uncolored void that can be corrected when setting `isoSurf` to .005.

![Example of Volume discrepancy, .003 (left) and .005 (right) in
ChimeraX Rainbow Palette](vmdfail.png "fig:"){width="2in"} ![Example of
Volume discrepancy, .003 (left) and .005 (right) in ChimeraX Rainbow
Palette](vmdcorrect.png "fig:"){width="1.925in"}

Moving Forward
==============

Administrative and Documentation
--------------------------------

### Cleanup 

There is a need to cleanup old files, commented out or redundant code
portions, and intermediate files.

### Consolidation

Code should be consolidated into functions for ease of use and
modularity.

### Correspondence and Accountability

Future GSOC internships should make further use of correspondence
methods such as Email, Github threads, Slack, and Zoom.

Conversion of IOData Data into Chemical JSON
--------------------------------------------

In order to improve interoperability with Avogadro2, it will be
necessary to be able to output IOData and Chemtools data into a
\"Chemical JSON\" a form of JSON file that outlines chemical information
such as bonds, atoms, as well as charge and molecular/atomic orbital
data. While the formatting of basic data such as bonds and atoms is
nontrivial, considerable work will have to be done to correctly format
quantum chemistry data.

Extended ChimeraX Functionality
-------------------------------

ChimeraX is a versatile tool capable of being controlled via python
scripting with different levels of granularity. ChimeraX Documentation
can be read from the program itself, or from the developer's website. It
is worthwhile to explore the ability of ChimeraX

Molecular Graphs of Bond Critical Points
----------------------------------------

Special Thanks
==============

I would like to extend great thanks to the following individuals for
assisting and granting me flexibility:

-   Farnaz Heidar Zadeh, Professor Queens University

-   Paul Ayers, Professor University of Toronto

-   Esteban Vohringer Martinez, Professor Universidad de Concepciòn

-   Victor Canogil

-   Toon Vestralen

-   Geoff Hutchinson, Avogadro2 Lead

-   Ali Tehrani

-   Elaine Meng, PhD, ChimeraX Development Team, University of
    California San Francisco

-   Anthony Schaefer, PhD, ChimeraX Development Team, University of
    California San Francisco

Appendix A: Github Fork
=======================
