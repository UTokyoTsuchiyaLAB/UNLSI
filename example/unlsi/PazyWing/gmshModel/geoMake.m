clear;
N=20;
N_rib = 3;
chord = 0.1;
span = 0.55;
airfoil = CST_airfoil([-0.1294 -0.0036 -0.0666], [0.206 0.2728 0.2292],0,N);
lc   = 1;                        % mesh size
point = airfoil(1:N,:)*chord;


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
% yList = linspace(0,span,N_rib*10);
% xList = linspace(0.4455*chord,0.5545*chord,N);
% for i=1:N_rib*10
%     fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+i,0.4455*chord,yList(i),0,lc);%Point(番号)= {x,y,z,メッシュサイズ}
% end

% for i=1:N
%     fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',(N+10)*N_rib+i,xList(i),span,0,lc);
% end

% for i=1:N_rib*10
%     fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',(N+10)*N_rib+N+i,0.5545*chord,yList(N_rib*10-i+1),0,lc);%Point(番号)= {x,y,z,メッシュサイズ}
% end
% for i=1:N
%     fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',(N+20)*N_rib+N+i,xList(N-i+1),0,0,lc);
% end

fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+1,0.4455*chord,0,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+2,0.5545*chord,0,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+3,0.5545*chord,0.1*k,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}
fprintf(fileID,'Point(%d)={%d,%d,%d,%d};\n',N*N_rib+4,0.4455*chord,0.1*k,0,lc);%Point(番号)= {x,y,z,メッシュサイズ}

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

fprintf(fileID,'Mesh.CharacteristicLengthMax = 0.2; // 最大メッシュサイズを設定\n Mesh.CharacteristicLengthMin = 0.1; // 最小メッシュサイズを設定\n');
fprintf(fileID,'Mesh.Algorithm = 1; // Delaunay法に基づく三角形メッシュ生成を強制\n');
%%%%%%%%%%面の作成
for k=0:N_rib-1
    fprintf(fileID,'Plane Surface(%d) = {%d};\n',k+1,k+1);
end
%center aluminun plate
fprintf(fileID,'Plane Surface(%d) = {%d};\n',N_rib+1,N_rib+1);

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

% for k=0:N_rib-1
%     fprintf(fileID,'Physical Surface("rib") = {%d};\n',k+1);
% end

fclose(fileID);

command = 'gmsh geofile.geo -2 -format msh2 -o airfoiltest0110.msh';
system(command);
% GmshRun;
