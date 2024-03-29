clear;
N=20;
N_rib = 3;
airfoil = CST_airfoil([-0.1294 -0.0036 -0.0666], [0.206 0.2728 0.2292],0,N);
lc   = 3;                        % mesh size

point = airfoil(1:N,:);


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

%%%%%%%%%%点列の作成
for k=0:N_rib-1
    y = 0.1*k;
    for i=1:N
    fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',i+N*k,point(i,1),y,point(i,2),lc);%Poit(番号)= {x,y,z,メッシュサイズ}
    end
end
%center aluminun plate
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+1,0.4455,0,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+2,0.5545,0,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+3,0.5545,0.1*k,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+4,0.4455,0.1*k,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}

%%%%%%%%%%線分の作成
for k=0:N_rib-1
    for i = 1:N
    fprintf(fileID,'Line(%d) = {%d,%d};\n',i+N*k,L(i,:));
    end
    L=N+L;
end
%center aluminun plate
fprintf(fileID,'Line(%d) = {%d,%d};\n',N*N_rib+1,N*N_rib+1,N*N_rib+2);
fprintf(fileID,'Line(%d) = {%d,%d};\n',N*N_rib+2,N*N_rib+2,N*N_rib+3);
fprintf(fileID,'Line(%d) = {%d,%d};\n',N*N_rib+3,N*N_rib+3,N*N_rib+4);
fprintf(fileID,'Line(%d) = {%d,%d};\n',N*N_rib+4,N*N_rib+4,N*N_rib+1);

%%%%%%%%%%線分のループ
for k=0:N_rib-1
    fprintf(fileID,'Curve Loop(%d)={',k+1);
    for i = 1:N-1
        fprintf(fileID,'-%d,',i+k*N);
    end
    fprintf(fileID,'%d};\n',N*(k+1));
end
%center aluminun plate
fprintf(fileID,'Curve Loop(%d)={%d,%d,%d,%d};\n',N_rib+1,N*N_rib+1,N*N_rib+2,N*N_rib+3,N*N_rib+4);

%%%%%%%%%%面の作成
for k=0:N_rib-1
    fprintf(fileID,'Plane Surface(%d) = {%d};\n',k+1,k+1);
end
%center aluminun plate
fprintf(fileID,'Plane Surface(%d) = {%d};\n',N_rib+1,N_rib+1);

%%%%%%%%%%面の物理的な境界条件
fprintf(fileID,'Physical Surface("rib") = {');
for k=1:N_rib-1
    fprintf(fileID,'%d,',k);
end
fprintf(fileID,'%d};\n',N_rib);
%center aluminun plate
fprintf(fileID,'Physical Surface("aluminun_plate") = {%d};\n',N_rib+1);

% for k=0:N_rib-1
%     fprintf(fileID,'Physical Surface("rib") = {%d};\n',k+1);
% end

fclose(fileID);

GmshRun;
