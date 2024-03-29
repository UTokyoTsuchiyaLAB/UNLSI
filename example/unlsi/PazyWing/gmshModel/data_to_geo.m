%-----------------------------------------------------------------------%
%----- data to .geo file --------------------------------------------%
%-----------------------------------------------------------------------%

clc
close all
clear

%-----------------------------------------------------------------------%
% define border
lc   = 3;                        % mesh size
Xc   = [ -25,25,-25,25];         % X - coordinates of corners
Zc   = [ 0,0,-50,-50];           % Z - coordinates of corners

A(1:4,1:5) = 0;
for i =1:4
   A(i,1)=i ; 
   A(i,2)=Xc(i);   
   A(i,3)=0;                      % Y-coordinate (X-Z plane)
   A(i,4)=Zc(i);
   A(i,5)=lc;
end

%----------------------------------------------------------------------%
%- lines that join the border points

Prev = [ 1;1;3;4];   % previous corner node
Next = [ 2;3;4;2];   % next corner node

L(1:4,1:3) = 0;
for i = 1:4
    L(i,1) = i;
    L(i,2) = Prev(i);
    L(i,3) = Next(i);
end
  
%-----------------------------------------------------------------------%
% internal points where finer mesh is needed

lc1   = 0.3;                     % mesh size
X    = [ -2,2,1,-1];             % X - coordinates of internal box
Z    = [ -1,-1,-3,-3];           % Z - coordinates of internal box

B(1:4,1:5) = 0;
for i =1:4
   B(i,1)=i+4 ; 
   B(i,2)=X(i);   
   B(i,3)=0;                      % Y-coordinate (X-Z plane)
   B(i,4)=Z(i);
   B(i,5)=lc1;
end

% GMsh geometry file is stored as .txt file
fileID = fopen('geofile.txt','w');
for i=1:4
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',A(i,:));
end
for i = 1:4
fprintf(fileID,'Line(%d) = {%d,%d};\n',L(i,:));
end
fprintf(fileID,'Line Loop(5) = {1,-4,-3,-2};\n');
fprintf(fileID,'Plane Surface(6) = {5};\n');
for i=1:4
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',B(i,:));
end
for i=1:4
fprintf(fileID,'Point{%d} In Surface{6};\n',B(i,1));
end
fclose(fileID);
type geofile.txt
uiopen('C:\ROOT GEOMETRY\RootSys Surface\geofile.txt',1)

%-----------------------------------------------------------------