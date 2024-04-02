command = 'gmsh geofile.geo -2 -format msh2 -o airfoiltest0110.msh';
% command = 'gmsh geofile.geo -2 -format msh -o airfoiltest0110.msh';
command2 = 'gmsh -refine airfoiltest0110.msh';
system(command);
disp('Meshing1 done');
% system(command2);
% disp('Meshing2 done');