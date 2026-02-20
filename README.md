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

python openmm_run.py -i step5_production.inp -p step3_input.psf -c step3_input.crd -t toppar.str --platform CPU -opdb prod_output.pdb -odcd prod_traj.dcd