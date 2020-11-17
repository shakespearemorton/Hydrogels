This program was created to explore how characteristics of a rough surface influence adsorption of monomers. Maximization of adsorbed monomers will increase the toughness of the adhesive, providing a strong and non-toxic alternative to surgical glues. 

To run the program, open runMe.sh and change L, M, or R. This will change the spherical harmonic used (L and M) along with the internal radius of the nanoparticle (R).The program is set by default to scale the particle size so that the maximum radius is always 10. By scaling the radius, it is ensured that an increase in surface area comes from increased roughness, and not just having a larger particle. This can be changed by setting P = 0, and then altering S. 

When the program is run, the matlab file, Characterize will give a text file containing:
Surface Area, Average-Mean-Curvature, RMS, and Sb
Sb is a shape factor that has been found to correlate to adsorption of polymers to the surface of nanoparticles. Specifically it is ( Average Roughness / Average Wavelength ) ^2. 

A particle is also generated for LAMMPS called system.txt using the python file doMinimize. There is also an option to create a VMD file in the XYZ format that has been hashed out, but is available. The particle will be made up out of sub-particles, the number of which is related to the Surface Area. To change what factor of the surface area it is related to, alter nParticles in doMinimize.

LAMMPS is then submitted to the HPC cluster for 40 hours to run a compression simulation with the particle between two hydrogels. There is an error that occurs frequently where the gel does not fully compress around the particle. The cause of this error has not been found. If this error occurs the simulation will still run to completion, so please check for reasonable results, and glance at the .atom file produced. 

Post_Process will run after the 40 hour cut-off has been hit, giving
1, Max Stress, Max Strain, Toughness, Number of monomers attached before deformation, effective interaction potential, and bridging parameter. 

All of these programs can be altered to analyze multiple particles by adding a loop, or having the program reference an excel sheet containing particle details. 

If any files become corrupted, unaltered files can be found in the Master folder. Please do not alter the files here unless absolutely necessary.