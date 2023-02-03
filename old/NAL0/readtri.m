function [tri, Verts, surfID] = readtri(filename)
fp = fopen(filename,'r');
buff = fscanf(fp,'%d %d',2);
nVerts = buff(1);
nTri = buff(2);

Verts = fscanf(fp,'%f %f %f',[3,nVerts])';
tri = fscanf(fp,'%d %d %d',[3,nTri])';
surfID = fscanf(fp,'%d',[1,nTri])';
x = Verts(:,1);
y = Verts(:,2);
z = Verts(:,3);
fclose all;
trimesh(tri,x,y,z)
axis equal
end