# Motion-Tree-project

This repository is aimed at reproducing motion tree's from molecular simulations. This is my research project at Birkbeck University. 

The motion tree python file, "motion_tree_full_nodes.py"
requires a ./data folder, where the MD simulation needs to be saved, and also set as the current working directory to run.

Currently, the file expects MD trajectory data, the example includes 53ligand.gro and 53alignedligand.xtc, available on the university cluster servers at: 
/d/projects/djosep08/3hpq/gromacs


data that has already had pre-selection for calpha atoms. The steps using gromacs are seen at the end of the readme.
It is possible to include all atoms for a potential AAMT, this can be implemented but not currently an option, 
supplying a full trajectory (of all atoms instead of just Calpha atoms), will still result in a calpha MT, due to the atom selection for the universe settings.

Required modules:
numpy
MDAnalysis
matplotlib.pyplot
scipy.cluster.hierarchy
calpha (found under the modules folder of the directory)

Settings for the 'graph settings' of the python file are available, this includes:
-removal or inclusion of the residual numbers on the x axis,
-size constraints,
-titles.





 ai was used for:
-the progress bars, 
-outputting cluster information for chimerax & numbering of the topmost nodes under 'topmost nodes' in calpha module.



To do:
-Upload the rmsd/rmsf functions, and the 'slicing' python files. although these were used for analysis
-AAMT can be implemented with changing of the selection of atoms, to instead account for all heavy atoms, instead of just the Calpha files.



Order of understanding: 


Paper 2013 = "Hierarchical Description and Extensive Classification of Protein Structural Changes by Motion Tree".
Includes computation of a matrix D, looking at the difference between 2 proteins.


Paper 2016 = "Motion Tree Delineates Hierarchical Structure  of Protein Dynamics Observed in Molecular  Dynamics Simulation", an extension, with a new matrix D




File reminders:
    - .xtc and .trr are both trajectory files, with xtc being more compressed but .trr inclues things like forces.
    - .tpr contains input for system, topology, coordinates/solvent box.
    - .gro/.pdb contain the reference structure, as an alternative to .tpr.
    - explore gmx-trjconv to clean trajectory file and then select only Calpha/backbone atoms. 

 group 0, group 1, group 2
solu(protein), solv(water), system
Files:
-f (trajectory, .xtc/.trr)
-s (structure+mass(db) .tpr, .gro,)
-n (index, .ndx optional)
-o output file


tips: https://epcced.github.io/20220421_GROMACS_introduction/04-post-analysis/index.html
https://manual.gromacs.org
https://www.youtube.com/watch?v=VLvB1vyltu8
https://www.blopig.com/blog/2025/08/taming-the-trajectory-beast-a-simpler-way-to-sample-your-md-simulations/
https://morphit-pro.cmp.uea.ac.uk/MorphItPro/faces/faces/about.xhtml?faces-redirect=true (server for visual changes)https://nglviewer.org/mdsrv/examples.html (view animations on the web)

Commands used on gmx:
1: (remove PBC)gmx trjconv -s step_10.tpr -f step5_10.xtc -o sys_md_nojump.xtc -pbc nojump -n index.ndx

For a MT, step 2 and 3 shoould be ignored due to potential dilution or alteration of results.
2: (Centering) gmx trjconv -s step_10.tpr -f sys_md_nojump.xtc -o sys_md_centered.xtc -center -pbc mol -ur compact -n index.ndx
3: (alighning) gmx trjconv -s step_10.tpr -f sys_md_centered.xtc -o sys_md_aligned.xtc -fit rot+trans -n index.ndx






https://www.researchgate.net/post/How_to_do_pca_analysis_of_c-alpha_atom_of_the_protein

Selecting C-alpha:
gmx convert-tpr

(maybe create .xtc of only c-alpha as well)

Creating index file:
gmx select -s protein.tpr (group 3 is the Calphas.)
>gmx make_ndx -f step5_10.gro -o ca.ndx
>1 & 3 (selected both protein and c-alpha)
> q (leave editor)

rmsd calc
>gmx rms -s step5_10.tpr -f step5_10.xtc -n ca.ndx -o rma_ca.xvg

gmx trjconv -s system.tpr -f trajectory.xtc -n index.ndx -o ca_trajectory.xtc

protein only .gro file
>gmx trjconv -s step5_10.tpr -f step5_10.gro -o protein.gro \-n index.ndx
# take 2 with removing rotation/translation
>gmx trjconv -s step5_10.tpr -f step5_10.xtc -n index.ndx -o fittedca.xtc -fit rot+trans