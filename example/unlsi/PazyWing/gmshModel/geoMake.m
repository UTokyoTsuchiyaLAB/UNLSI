clear;
N=100;
airfoil = CST_airfoil([-0.1294 -0.0036 -0.0666], [0.206 0.2728 0.2292],0,N);
lc   = 3;                        % mesh size
% formatSpecPoint = 'Point(%4.0f)={%4.2f,%4.2f,}'

point = airfoil(1:N,:);
% FILEPATH: /c:/Users/aimit/git/UAVdesignTools/geoMake.m
% Create a matrix L of size N-by-2 filled with zeros.
L = zeros([N,2]);
for i = 1:N-1
    L(i,:) = [i,i+1];
end
L(N,:) = [N,1];
% loop = linspace(1,N,N);

% GMsh geometry file is stored as .txt file
%https://kawakamik.com/post/gmsh/
fileID = fopen('geofile.geo','w');
fprintf(fileID,'// Gmsh project created on Fri Dec 20 20:36:50 2019 \n SetFactory("OpenCASCADE");\n');
for i=1:N
fprintf(fileID,'Point(%d)={%d,%d,0,3};\n',i,point(i,:));%Poit(番号)= {x,y,z,メッシュサイズ}
end
for i = 1:N
fprintf(fileID,'Line(%d) = {%d,%d};\n',i,L(i,:));
end
% fprintf(fileID,'Curve Loop(1) ={'+linspace(1,N,N)+'};\n');
fprintf(fileID,'Curve Loop(1)={');
for i = 1:N-1
    fprintf(fileID,'-%d,',i);
end
fprintf(fileID,'%d};\n',N);

fprintf(fileID,'Plane Surface(1) = {1};\n');
fprintf(fileID,'Physical Surface("rib1") = {1};');

% for i=1:4
% fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',B(i,:));
% end
% for i=1:4
% fprintf(fileID,'Point{%d} In Surface{6};\n',B(i,1));
% end
fclose(fileID);
% type geofile.geo
% uiopen('C:\ROOT GEOMETRY\RootSys Surface\geofile.txt',1)

GmshRun;
