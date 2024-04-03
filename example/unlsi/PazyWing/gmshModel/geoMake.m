clear;
N=100;
N_rib = 15;
chord = 0.1;
span = 0.55;
airfoil = CST_airfoil([-0.1294 -0.0036 -0.0666], [0.206 0.2728 0.2292],0,N);
lc   = 3;                        % mesh size
point = airfoil(1:N,:)*chord;
yRib = linspace(0.005,span,N_rib);

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
    y = yRib(k+1);
    for i=1:N
    fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',i+N*k,point(i,1),y,point(i,2),lc);%Poit(番号)= {x,y,z,メッシュサイズ}
    end
end
%center aluminun plate
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+1,0.02,0,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+2,0.08,0,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+3,0.08,span,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+4,0.02,span,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}

%tip rod
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+5,-0.1,0.55,-0.005,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+6,-0.1,0.555,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+7,-0.1,0.55,0.005,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+8,-0.1,0.545,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+9,0.2,0.55,-0.005,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+10,0.2,0.555,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+11,0.2,0.55,0.005,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+12,0.2,0.545,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}

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

%tip rod
for k = 0:1
    for i = 1:3
        fprintf(fileID,'Line(%d) = {%d,%d};\n',N*N_rib+4+i+4*k,N*N_rib+4+i+4*k,N*N_rib+i+5+4*k);
    end
    fprintf(fileID,'Line(%d) = {%d,%d};\n',N*N_rib+8+4*k,N*N_rib+8+4*k,N*N_rib+5+4*k);
end

for i=1:4
    fprintf(fileID,'Line(%d) = {%d,%d};\n',N*N_rib+12+i,N*N_rib+4+i,N*N_rib+8+i);
end

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

%tip rod
for k=1:2
    fprintf(fileID,'Curve Loop(%d)={%d,%d,%d,%d};\n',N_rib+k+1,N*N_rib+4*k+1,N*N_rib+4*k+2,N*N_rib+4*k+3,N*N_rib+4*k+4);
end
for k=1:4
    fprintf(fileID,'Curve Loop(%d)={%d,%d,%d,%d};\n',N_rib+3+k,N*N_rib+4+k,N*N_rib+5+k,N*N_rib+10+k,N*N_rib+9+k);
end


fprintf(fileID,'Mesh.CharacteristicLengthMax = 0.02; // 最大メッシュサイズを設定\n Mesh.CharacteristicLengthMin = 0.01; // 最小メッシュサイズを設定\n');
fprintf(fileID,'Mesh.Algorithm = 1; // Delaunay法に基づく三角形メッシュ生成を強制\n');
%%%%%%%%%%面の作成
for k=0:N_rib-1
    fprintf(fileID,'Plane Surface(%d) = {%d};\n',k+1,k+1);
end
%center aluminun plate
fprintf(fileID,'Plane Surface(%d) = {%d};\n',N_rib+1,N_rib+1);

%tip rod
for k=1:6
    fprintf(fileID,'Plane Surface(%d) = {%d};\n',N_rib+k+1,N_rib+k+1);
end

% fprintf(fileID,'Mesh.CharacteristicLengthMax = 0.2; // 最大メッシュサイズを設定 \n Mesh.CharacteristicLengthMin = 0.1; // 最小メッシュサイズを設定\n');
% fprintf(fileID,'Mesh.Algorithm = 2; // Delaunay法に基づく三角形メッシュ生成を強制\n');
fprintf(fileID,'Mesh.RecombineAll = 0;\n'); % 四角形要素に変換

%%%%%%%%%%面の物理的な境界条件
fprintf(fileID,'Physical Surface("rib") = {');
for k=1:N_rib-1
    fprintf(fileID,'%d,',k);
end
fprintf(fileID,'%d};\n',N_rib);
%center aluminun plate
fprintf(fileID,'Physical Surface("aluminun_plate") = {%d};\n',N_rib+1);

%tip rod
surfGroup = linspace(N_rib+2,N_rib+7,6);
fprintf(fileID,'Physical Surface("tip_rod") = {%d,%d,%d,%d,%d,%d};\n',surfGroup);

fclose(fileID);

command = 'gmsh geofile.geo -2 -format msh2 -o PazyWing.msh';
system(command);
system('gmsh PazyWing.msh');
% GmshRun;
