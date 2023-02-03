function [verts,tri,surfID,surfThn,surfE,surfRho] = NAL02tri(vehicle,HR)
%wgs出力
O = [0 0 0];
surfData = surface2net(1,vehicle.surface.nose{1},vehicle.surface.nose{2},vehicle.surface.nose{3},O,[1 0],HR);
surfData = surface2net(1,vehicle.surface.wingkink{1},vehicle.surface.wingkink{2},vehicle.surface.wingkink{3},O,[1 0],HR,surfData);
surfData = surface2net(1,vehicle.surface.wingout{1},vehicle.surface.wingout{2},vehicle.surface.wingout{3},O,[1 0],HR,surfData);
surfData = surface2net(1,vehicle.surface.body1{1},vehicle.surface.body1{2},vehicle.surface.body1{3},O,[1 0],HR,surfData);
surfData = surface2net(1,vehicle.surface.body2{1},vehicle.surface.body2{2},vehicle.surface.body2{3},O,[1 0],HR,surfData);
surfData = surface2net(1,vehicle.surface.body3{1},vehicle.surface.body3{2},vehicle.surface.body3{3},O,[0 0],HR,surfData);
surfData = surface2net(1,vehicle.surface.vert{1},vehicle.surface.vert{2},vehicle.surface.vert{3},O,[1 0],HR,surfData);
surfData = surface2net(1,vehicle.surface.wingtip{1},vehicle.surface.wingtip{2},vehicle.surface.wingtip{3},O,[0 0],HR,surfData);
surfData = surface2net(1,vehicle.surface.verttip{1},vehicle.surface.verttip{2},vehicle.surface.verttip{3},O,[0 0],HR,surfData);
surfData = surface2net(1,vehicle.surface.base{1},vehicle.surface.base{2},vehicle.surface.base{3},O,[1 0],HR,surfData);
surfData = surface2net(18,vehicle.wake.wingkink{1},vehicle.wake.wingkink{2},vehicle.wake.wingkink{3},O,[1 0],HR,surfData);
surfData = surface2net(18,vehicle.wake.wingout{1},vehicle.wake.wingout{2},vehicle.wake.wingout{3},O,[1 0],HR,surfData);
surfData = surface2net(18,vehicle.wake.vert{1},vehicle.wake.vert{2},vehicle.wake.vert{3},O,[1 0],HR,surfData);
%翼構造
%前縁から何番目にスパーを通すか
n_foil = 12;
sparNo = [4 8 10];
for i= 1:size(sparNo,2)
    for j =1:3
        wing.str{i,j}=[];
        for iter = 1:size(vehicle.surface.wing)
            wing.str{i,j} = [wing.str{i,j},[vehicle.surface.wing{j}(n_foil+sparNo(i),:);vehicle.surface.wing{j}(n_foil-sparNo(i),:)]];
        end
    end
end
iter = i+1;
ribStep = 1;
for i = 1:ribStep:size(vehicle.surface.wing{1},2)
    for j =1:3
        wing.str{iter,j} = [vehicle.surface.wing{j}(1:n_foil,i),flipud(vehicle.surface.wing{j}(n_foil:end,i))];
    end
    iter = iter+1;
end
%垂直尾翼リブ
for i = 1:size(vehicle.surface.vert{1},2)
    wing.str{iter,1} = [vehicle.surface.vert{1}(1:4,i)';flipud(vehicle.surface.vert{1}(4:end,i))'];
    wing.str{iter,3} = [vehicle.surface.vert{3}(1:4,i)';flipud(vehicle.surface.vert{3}(4:end,i))'];
    wing.str{iter,2} = [vehicle.surface.vert{2}(1:4,i)';flipud(vehicle.surface.vert{2}(4:end,i))'];
    iter = iter+1;
end
%}
n_str = size(wing.str,1);
for i = 1:n_str
    surfData = surface2net(-1,wing.str{i,1},wing.str{i,2},wing.str{i,3},O,[0,1],HR,surfData);
end
%}
n_wake = 3;
[verts, tri, surfID]=surf2tri(surfData,n_wake,n_str,1);

for i = 1:numel(surfID)
    if surfID(i,1) ==10
        surfID(i,1) = 101;
    end
    if surfID(i,1) ==60
        surfID(i,1) = 102;
    end
end
E_material =  73500000000;rho_material = 2770;stress_yield=470*10^6; %A2000系
surfE = ones(size(surfID,1),1).*E_material;
surfRho = ones(size(surfID,1),1).*rho_material;
surfThn = ones(size(surfID,1),1).*0.0001;
end

function [surfData] = surface2net(netnum,dataX,dataY,dataZ,O,flip,HR,surfData)
	dataX = dataX.*HR-O(1);
	dataY = dataY.*HR-O(2);
	dataZ = dataZ.*HR-O(3);
	
	if flip(1) ~=0
		dataX=flipud(dataX);
		dataY=flipud(dataY);
		dataZ=flipud(dataZ);
	end
	if flip(2) ~=0
		dataX=fliplr(dataX);
		dataY=fliplr(dataY);
		dataZ=fliplr(dataZ);
    end
    if nargin==8
        no = size(surfData,2)+1;
    else
		no = 1;
    end
    surfData{no}.netnum = netnum;
    surfData{no}.cols = size(dataX,2);
    surfData{no}.rows = size(dataX,1);
	
    iter = 1;
    for j_iter = 1:size(dataX,2)
        for i_iter =1:size(dataX,1)
            surfData{no}.pnt(iter,:) = [dataX(i_iter,j_iter),dataY(i_iter,j_iter),dataZ(i_iter,j_iter)];
            iter = iter+1;
        end
    end
end

function [verts, tri, GRP]=surf2tri(surfData,n_wake,n_str,isSym)
    k=1;
    for nm_sf = 1:size(surfData,2)
        nsf = surfData{nm_sf}.cols;
        msf = surfData{nm_sf}.rows;
        iter = 1;
        for i = 1:nsf
            for j = 1:msf
                vertexbuff = surfData{nm_sf}.pnt(iter,:);
                surfaceX{nm_sf}(j,i) = vertexbuff(1);
                surfaceY{nm_sf}(j,i) = vertexbuff(2);
                surfaceZ{nm_sf}(j,i) = vertexbuff(3);
                surfaceID{nm_sf}(j,i) = k;
                k=k+1;
                if isSym == 1
                    invID{nm_sf}(j,i) = k;
                    k=k+1;
                end
                iter=iter+1;
            end
        end
    end


    %IDをセットする
    sfX = size(surfaceX,2);
    GRP = [];
    verts = [];
    tri=[];
    %正転
    for i = 1:sfX-n_wake-n_str
        for j = 1:size(surfaceX{i},2)-1
            for k = 1:size(surfaceX{i},1)-1
                [tri GRP] = addtri([i,j,k],[i,j+1,k],[i,j,k+1],[i,j,k+1],[i,j+1,k],[i,j+1,k+1],surfaceX,surfaceY,surfaceZ,surfaceID,tri,GRP,0,0,n_wake,0,n_str);
            end
        end
    end
    %反転
    if isSym ==1
        for i = 1:sfX-n_wake-n_str
            for j = 1:size(surfaceX{i},2)-1
                for k = 1:size(surfaceX{i},1)-1
                    [tri GRP] = addtri([i,j+1,k],[i,j,k],[i,j,k+1],[i,j+1,k],[i,j,k+1],[i,j+1,k+1],surfaceX,surfaceY,surfaceZ,invID,tri,GRP,1,0,n_wake,0,n_str);
                end
            end
        end
    end
    %wake
    if not(n_wake==0)
        for i = sfX-n_wake-n_str+1:sfX-n_str
            for j = 1:size(surfaceX{i},2)-1
                for k = 1:size(surfaceX{i},1)-1
                    [tri GRP] = addtri([i,j,k],[i,j+1,k],[i,j,k+1],[i,j,k+1],[i,j+1,k],[i,j+1,k+1],surfaceX,surfaceY,surfaceZ,surfaceID,tri,GRP,0,1,n_wake,0,n_str);
                end
            end
        end
        %反転
        if isSym ==1
            for i = sfX-n_wake-n_str+1:sfX-n_str
                for j = 1:size(surfaceX{i},2)-1
                    for k = 1:size(surfaceX{i},1)-1
                        [tri GRP] = addtri([i,j+1,k],[i,j,k],[i,j,k+1],[i,j+1,k],[i,j,k+1],[i,j+1,k+1],surfaceX,surfaceY,surfaceZ,invID,tri,GRP,1,1,n_wake,0,n_str);
                    end
                end
            end
        end
    end
    %structure
    if not(n_str==0)
        for i = sfX-n_str+1:sfX
            for j = 1:size(surfaceX{i},2)-1
                for k = 1:size(surfaceX{i},1)-1
                    [tri GRP] = addtri([i,j,k],[i,j+1,k],[i,j,k+1],[i,j,k+1],[i,j+1,k],[i,j+1,k+1],surfaceX,surfaceY,surfaceZ,surfaceID,tri,GRP,0,1,n_wake,1,n_str);
                end
            end
        end
        %反転
        if isSym ==1
            for i = sfX-n_str+1:sfX
                for j = 1:size(surfaceX{i},2)-1
                    for k = 1:size(surfaceX{i},1)-1
                        [tri GRP] = addtri([i,j+1,k],[i,j,k],[i,j,k+1],[i,j+1,k],[i,j,k+1],[i,j+1,k+1],surfaceX,surfaceY,surfaceZ,invID,tri,GRP,1,1,n_wake,1,n_str);
                    end
                end
            end
        end
    end
    %verts
    for i=1:sfX-n_wake
        for j = 1:size(surfaceX{i},2)
            for k = 1:size(surfaceX{i},1)
                verts = [verts;surfaceX{i}(k,j),surfaceY{i}(k,j),surfaceZ{i}(k,j)];
                if isSym == 1
                    verts = [verts;surfaceX{i}(k,j),-surfaceY{i}(k,j),surfaceZ{i}(k,j)];
                end
            end
        end
    end
    if not(n_wake == 0)
        for i=sfX-n_wake+1:sfX
            for j = 1:size(surfaceX{i},2)
                for k = 1:size(surfaceX{i},1)
                    verts = [verts;surfaceX{i}(k,j),surfaceY{i}(k,j),surfaceZ{i}(k,j)];
                    if isSym == 1
                        verts = [verts;surfaceX{i}(k,j),-surfaceY{i}(k,j),surfaceZ{i}(k,j)];
                    end
                end
            end
        end
    end
    %同じ点をマージ
    tlist = 0;
    for i = 1:size(verts,1)
        Qi = verts(i,:);
        for j = i:size(verts,1)
            Qj = verts(j,:);
            if all(norm(Qi-Qj)<sqrt(eps)/1000) && not(i==j) && not(any(tlist == i))
                tri(tri==j) = i;
                tlist = [tlist;j];
            end
        end
    end
end   

function [tri, GRP] = addtri(nP1,nP2,nP3,nQ1,nQ2,nQ3,surfaceX,surfaceY,surfaceZ,surfaceID,tri,GRP,ismirror,iswake,n_wake,isstr,n_str)
P1 = [surfaceX{nP1(1)}(nP1(3),nP1(2)),surfaceY{nP1(1)}(nP1(3),nP1(2)),surfaceZ{nP1(1)}(nP1(3),nP1(2))];
P2 = [surfaceX{nP2(1)}(nP2(3),nP2(2)),surfaceY{nP2(1)}(nP2(3),nP2(2)),surfaceZ{nP2(1)}(nP2(3),nP2(2))];
P3 = [surfaceX{nP3(1)}(nP3(3),nP3(2)),surfaceY{nP3(1)}(nP3(3),nP3(2)),surfaceZ{nP3(1)}(nP3(3),nP3(2))];
if norm(P1-P2)<sqrt(eps)/1000 || norm(P2-P3)<sqrt(eps)/1000 || norm(P3-P1)<sqrt(eps)/1000
    Pis2D = 1;
else
    Pis2D = 0;
end

Q1 = [surfaceX{nQ1(1)}(nQ1(3),nQ1(2)),surfaceY{nQ1(1)}(nQ1(3),nQ1(2)),surfaceZ{nQ1(1)}(nQ1(3),nQ1(2))];
Q2 = [surfaceX{nQ2(1)}(nQ2(3),nQ2(2)),surfaceY{nQ2(1)}(nQ2(3),nQ2(2)),surfaceZ{nQ2(1)}(nQ2(3),nQ2(2))];
Q3 = [surfaceX{nQ3(1)}(nQ3(3),nQ3(2)),surfaceY{nQ3(1)}(nQ3(3),nQ3(2)),surfaceZ{nQ3(1)}(nQ3(3),nQ3(2))];
if norm(Q1-Q2)<sqrt(eps)/1000 || norm(Q2-Q3)<sqrt(eps)/1000 || norm(Q3-Q1)<sqrt(eps)/1000
    Qis2D = 1;
else
    Qis2D = 0;
end
pmat = [P1;P2;P3];
qmat = [Q1;Q2;Q3];
%[buff PSI] = sort(pmat);
%[buff QSI] = sort(qmat);
%pmat = pmat(PSI(:,2),:);
%qmat = qmat(QSI(:,2),:);
for i = 1:3
    if all(pmat(1,:)==qmat(i,:))
        pmat(1,:) = [];
        qmat(i,:) = [];
        break;
    end
end
if size(pmat,1) == 2
    for i = 1:2
        if all(pmat(1,:)==qmat(i,:))
            pmat(1,:) = [];
            qmat(i,:) = [];
            break;
        end
    end
end
if size(pmat,1) == 1
    if all(pmat(1,:)==qmat(1,:))
        pmat(1,:) = [];
        qmat(1,:) = [];
    end
end
if isempty(pmat)
    isSame = 1;
else
    isSame = 0;
end

i = nP1(1);%必ず同一
sfX = size(surfaceID,2);
if isstr == 1
   addGRP = 301; 
elseif ismirror == 1 && iswake == 1
    addGRP = i-(sfX-n_wake-n_str)+250;
elseif ismirror == 0 && iswake == 1
    addGRP = i-(sfX-n_wake-n_str)+200;
elseif ismirror == 0 && iswake == 0
    addGRP = i;
    %addGRP = 1;
elseif ismirror == 1 && iswake == 0
    addGRP = 50+i;
    %addGRP=2;
end
    
if isSame == 1 && Pis2D == 0 && Qis2D==1
    tri =[tri;surfaceID{nP1(1)}(nP1(3),nP1(2)),surfaceID{nP2(1)}(nP2(3),nP2(2)),surfaceID{nP3(1)}(nP3(3),nP3(2))];
    GRP = [GRP;addGRP];
elseif isSame == 1 && Pis2D == 1 && Qis2D ==0
    tri =[tri;surfaceID{nQ1(1)}(nQ1(3),nQ1(2)),surfaceID{nQ2(1)}(nQ2(3),nQ2(2)),surfaceID{nQ3(1)}(nQ3(3),nQ3(2))];
    GRP = [GRP;addGRP];
elseif isSame == 1 && Pis2D == 0 && Qis2D ==0
    tri =[tri;surfaceID{nP1(1)}(nP1(3),nP1(2)),surfaceID{nP2(1)}(nP2(3),nP2(2)),surfaceID{nP3(1)}(nP3(3),nP3(2))];
    GRP = [GRP;addGRP];
elseif isSame == 0 && Pis2D == 0 && Qis2D ==1
    tri =[tri;surfaceID{nP1(1)}(nP1(3),nP1(2)),surfaceID{nP2(1)}(nP2(3),nP2(2)),surfaceID{nP3(1)}(nP3(3),nP3(2))];
    GRP = [GRP;addGRP];
elseif isSame == 0 && Pis2D == 1 && Qis2D ==0
    tri =[tri;surfaceID{nQ1(1)}(nQ1(3),nQ1(2)),surfaceID{nQ2(1)}(nQ2(3),nQ2(2)),surfaceID{nQ3(1)}(nQ3(3),nQ3(2))];
    GRP = [GRP;addGRP];
elseif isSame == 0 && Pis2D == 0 && Qis2D ==0
    tri =[tri;surfaceID{nP1(1)}(nP1(3),nP1(2)),surfaceID{nP2(1)}(nP2(3),nP2(2)),surfaceID{nP3(1)}(nP3(3),nP3(2))];
    GRP = [GRP;addGRP];
    tri =[tri;surfaceID{nQ1(1)}(nQ1(3),nQ1(2)),surfaceID{nQ2(1)}(nQ2(3),nQ2(2)),surfaceID{nQ3(1)}(nQ3(3),nQ3(2))];
    GRP = [GRP;addGRP];
end
end
