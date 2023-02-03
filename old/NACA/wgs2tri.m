function [verts,tri,GRP,surfaceX,surfaceY,surfaceZ]=wgs2tri(wgsname,triname,n_str,n_wake,isSym)
    %clear all
    %n_wake = 1;
    %isSym = 1;
    %wgsname = 'HSBJg.wgs';
    %triname = 'HSBJ.tri';
    %------------------------
    %wgsì«Ç›çûÇ›
    %-----------------------
    fp=fopen(wgsname,'rt');
    buff = fgetl(fp);
    nm_sf = 1;
    k = 1;
    while 1
        buff = fgetl(fp);
        if buff ==-1
            break;
        end
        matbuff = fscanf(fp,'%d',[3,1]);
        nsf = matbuff(2);
        msf = matbuff(3);
        buff = fgetl(fp);
        for i = 1:nsf
            for j = 1:msf
                vertexbuff = fscanf(fp,'%f',[3,1])';
                surfaceX{nm_sf}(j,i) = vertexbuff(1);
                surfaceY{nm_sf}(j,i) = vertexbuff(2);
                surfaceZ{nm_sf}(j,i) = vertexbuff(3);
                surfaceID{nm_sf}(j,i) = k;
                k=k+1;
                if isSym == 1
                    invID{nm_sf}(j,i) = k;
                    k=k+1;
                end
            end
        end
        buff = fgetl(fp);
        nm_sf =nm_sf+1;
    end
    fclose(fp);

%IDÇÉZÉbÉgÇ∑ÇÈ
    sfX = size(surfaceX,2);
    GRP = [];
    verts = [];
    tri=[];
    %ê≥ì]
    for i = 1:sfX-n_wake-n_str
        for j = 1:size(surfaceX{i},2)-1
            for k = 1:size(surfaceX{i},1)-1
                [tri GRP] = addtri([i,j,k],[i,j+1,k],[i,j,k+1],[i,j,k+1],[i,j+1,k],[i,j+1,k+1],surfaceX,surfaceY,surfaceZ,surfaceID,tri,GRP,0,0,n_wake,0,n_str);
            end
        end
    end
    %îΩì]
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
        %îΩì]
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
        %îΩì]
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
    %ìØÇ∂ì_ÇÉ}Å[ÉW
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

i = nP1(1);%ïKÇ∏ìØàÍ
sfX = size(surfaceID,2);
if isstr == 1
   addGRP = i+(sfX-n_str)+300; 
elseif ismirror == 1 && iswake == 1
    addGRP = i-(sfX-n_wake-n_str)+250;
elseif ismirror == 0 && iswake == 1
    addGRP = i-(sfX-n_wake-n_str)+200;
elseif ismirror == 0 && iswake == 0
    addGRP = i;
    %addGRP = 1;
elseif ismirror == 1 && iswake == 0
    addGRP = sfX-n_wake+i;
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