command = 'gmsh geofile.geo -2 -format msh -o airfoiltest0110.msh';
command2 = 'gmsh -refine airfoiltest.msh';
system(command);