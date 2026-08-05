# Motion-Tree-project
This repository is aimed at reproducing motion tree's from molecular simulations. This is my research project at Birkbeck University. 


To do:

- Create simulation data from https://www.pnas.org/doi/10.1073/pnas.0608432104 molucule.
- ***Breakdown of the maths and submit document to Mark. 
- https://github.com/khinsen/MMTK use toolkit.
- pick a pipeline to use!



Order of understanding: 
Paper 2013 = "Hierarchical Description and Extensive Classification of Protein Structural Changes by Motion Tree".
Includes computation of a matrix D, looking at the difference between 2 proteins.


Paper 2016 = "Motion Tree Delineates Hierarchical Structure  of Protein Dynamics Observed in Molecular  Dynamics Simulation", an extension, with a new matrix D


Questions to ask:
- what solvent should proteins be studied in? water? the solvent presented in the puplication? 

source /usr/local/gromacs/bin/GMXRC # location of gmx



April To do:
    - .xtc and .trr are both trajectory files, with xtc being more compressed but .trr inclues things like forces.
    - .tpr contains input for system, topology, coordinates/solvent box.
    - .gro/.pdb contain the reference structure, as an alternative to .tpr.
    - explore gmx-trjconv to clean trajectory file and then select only Calpha/backbone atoms. 


tips: https://epcced.github.io/20220421_GROMACS_introduction/04-post-analysis/index.html
https://manual.gromacs.org
https://www.youtube.com/watch?v=VLvB1vyltu8
https://www.blopig.com/blog/2025/08/taming-the-trajectory-beast-a-simpler-way-to-sample-your-md-simulations/
https://morphit-pro.cmp.uea.ac.uk/MorphItPro/faces/faces/about.xhtml?faces-redirect=true (server for visual changes)https://nglviewer.org/mdsrv/examples.html (view animations on the web)

Commands used on gmx:
1: (remove PBC)gmx trjconv -s step_10.tpr -f step5_10.xtc -o sys_md_nojump.xtc -pbc nojump -n index.ndx
2: (Centering) gmx trjconv -s step_10.tpr -f sys_md_nojump.xtc -o sys_md_centered.xtc -center -pbc mol -ur compact -n index.ndx
3: (alighning) gmx trjconv -s step_10.tpr -f sys_md_centered.xtc -o sys_md_aligned.xtc -fit rot+trans -n index.xtc


scp:
scp djosep08@ssh.cryst.bbk.ac.uk:/d/user6/djosep08/Projects/project_gromacs/gromacs/rmsd.xvg .
scp -r step34 djosep08@ssh.cryst.bbk.ac.uk:Projects


 group 0, group 1, group 2
solu(protein), solv(water), system
Files:
-f (trajectory, .xtc/.trr)
-s (structure+mass(db) .tpr, .gro,)
-n (index, .ndx optional)
-o output file

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



vmd: 
was successful, simulation obtained.
