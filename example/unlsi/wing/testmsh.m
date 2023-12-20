%------------------------------------------------------------------------%
%------ Gmsh to Matlab script: Import mesh to matlab---------------------%
%------------------------------------------------------------------------%
clc
close all
clear 
%-----------------------------------------------------------------------%
% dlmread(filename,delimiter,[R1 C1 R2 C2]) reads only the range 
% bounded by row offsets R1 and R2 and column offsets C1 and C2.
%-----------------------------------------------------------------------%
file    =  ('wing_WingGeom_Struct0.msh');
[srf,vtx,fc,bc,simp,edg,mat_ind,phys_names] = gmsh_read_mesh(file);
%---- visualize in matlab ---------------------
figure(1)
trimesh(srf(fc==2,:),vtx(:,1),vtx(:,2),vtx(:,3))
xlabel('X','fontsize',14)
ylabel('Y','fontsize',14)
title('GMsh to MATLAB import','fontsize',14)
fh = figure(1);
set(fh, 'color', 'white'); 
%-------------------------------------------------------------------------
