function frame = UNLSI_makeFrame(verts,tri,surfID,LLTwakeID,onlySuperSonic)
%surfIDはパネルの要素を記述する
%1~100   : ボディパネル(空力+構造)
%101~200 : ベース面パネル(圧力分布統計推算のみ)
%201~300 : 後流パネル
%301~400 : 構造パネル(構造)
%REFS 
%[XREF YREF ZREF;
% SREF BREF CREF]
if nargin == 3
    LLTwakeID = [];
    onlySuperSonic = 0;
end
if nargin == 4
    onlySuperSonic = 0;
end

%データの格納
frame.verts = verts;
frame.tri = tri;
frame.surfID = surfID;

%パネルの種類分け
frame.isBody = and(frame.surfID>=1,frame.surfID<=100);
frame.isBase = and(frame.surfID>=101,frame.surfID<=200);
frame.isWake = and(frame.surfID>=201,frame.surfID<=300);
frame.isStr = and(frame.surfID>=301,frame.surfID<=400);
frame.nPanel = size(frame.tri,1);
frame.nbPanel = sum(frame.isBody);
frame.tri2ubody = zeros(frame.nPanel,1);
iter = 1;
for i = 1:frame.nPanel
    if frame.isBody(i,1)==1 
        frame.tri2ubody(i,1) = iter;
        iter = iter +1;
    end
end

%法線ベクトルとパネル面積の計算
frame.nTri = size(frame.tri,1);
frame.nVerts = size(frame.verts,1);

for i = 1:size(frame.tri,1)
    [buff1, buff2, frame.n(i,:)] = vartex_iter(frame.verts(frame.tri(i,1),:),frame.verts(frame.tri(i,2),:),frame.verts(frame.tri(i,3),:));
end
if onlySuperSonic == 0
    %
    %Connectivityの作成
    %Conectivity Check
    for i = 1:frame.nPanel
        if frame.isWake(i,1)==0
            nb = [];
            ea = [];
            [buff1 buff2] = neighbourcheck(i,[frame.tri(i,1),frame.tri(i,2)],frame);
            nb = [nb buff1];
            ea = [ea buff2];
            [buff1 buff2] = neighbourcheck(i,[frame.tri(i,2),frame.tri(i,3)],frame);
            nb = [nb buff1];
            ea = [ea buff2];
            [buff1 buff2] = neighbourcheck(i,[frame.tri(i,3),frame.tri(i,1)],frame);
            nb = [nb buff1];
            ea = [ea buff2];
            neighbour{i} = nb;
            edgeangle{i} = ea;
            nbw = [];
            eaw = [];
            [nbbuff, eabuff] = neighbourwakecheck(i,[frame.tri(i,1),frame.tri(i,2)],frame);
            nbw = [nbw nbbuff];
            eaw = [eaw eabuff];
            [nbbuff, eabuff] = neighbourwakecheck(i,[frame.tri(i,2),frame.tri(i,3)],frame);
            nbw = [nbw nbbuff];
            eaw = [eaw eabuff];
            [nbbuff, eabuff] = neighbourwakecheck(i,[frame.tri(i,3),frame.tri(i,1)],frame);
            nbw = [nbw nbbuff];
            eaw = [eaw eabuff];
            wakeneibour{i} = nbw;
            wakeea{i} = eaw;
            if not(isempty(nbw))
                attachwake(i,1) = 1;
            else
                attachwake(i,1) = 0;
            end
        else
            attachwake(i,1) = 0;
        end
    end
    %Cluster探索。
    nCluster = round(frame.nbPanel/5)+20;
    for i = 1:frame.nPanel
        if frame.isBody(i,1)==1
            toAdd = i;
            j = 1;
            try
                for iter = 1:nCluster
                    for k = 1:numel(neighbour{toAdd(j)})
                       if edgeangle{toAdd(j)}(k)<=60*pi/180 &&acos(dot(frame.n(i,:),frame.n(neighbour{toAdd(j)}(k),:)))<=60*pi/180 && not(any(toAdd==neighbour{toAdd(j)}(k)))
                       %if edgeangle{toAdd(j)}(k)<=60*pi/180 && not(any(toAdd==neighbour{toAdd(j)}(k)))
                           toAdd = [toAdd,neighbour{toAdd(j)}(k)];
                       end
                    end
                    j = j+1;
                    if numel(toAdd)>nCluster
                        break
                    end
                end
                frame.cluster{i} = [i,toAdd(2:nCluster)];
            catch
                frame.cluster{i} = [i,toAdd(2:end)];
            end
        end
    end


    %attachWakeしているパネルの中からwakeedgeを探し出す
    frame.wakelines = [];
    if not(sum(frame.isWake)==0)
        nodeID = [];
        bpanelID = [];
        wpanelID = [];
        for i = 1:size(frame.tri,1)
            if frame.isBody(i,1)==1
                if attachwake(i,1) == 1
                    boolbuff = neighbourcheckID(i,[frame.tri(i,1),frame.tri(i,2)],wakeneibour{i},frame);
                    if boolbuff == 1
                        nodeID = [nodeID; frame.tri(i,1),frame.tri(i,2)];
                        bpanelID = [bpanelID;i];
                        wpanelID = [wpanelID;wakeneibour{i}];
                    end
                    boolbuff = neighbourcheckID(i,[frame.tri(i,2),frame.tri(i,3)],wakeneibour{i},frame);
                    if boolbuff == 1
                        nodeID = [nodeID; frame.tri(i,2),frame.tri(i,3)];
                        bpanelID = [bpanelID;i];
                        wpanelID = [wpanelID;wakeneibour{i}];
                    end
                    boolbuff = neighbourcheckID(i,[frame.tri(i,3),frame.tri(i,1)],wakeneibour{i},frame);
                    if boolbuff == 1
                        nodeID = [nodeID; frame.tri(i,3),frame.tri(i,1)];
                        bpanelID = [bpanelID;i];
                        wpanelID = [wpanelID;wakeneibour{i}];
                    end
                end
            end
        end
        %wakeのマージ
        j=1;
        i=2;
        while(1)
           if or(nodeID(1,1) == nodeID(i,1) && nodeID(1,2) == nodeID(i,2),nodeID(1,2) == nodeID(i,1) && nodeID(1,1) == nodeID(i,2))
               frame.wakelines{j}.nodeID = [nodeID(1,1) nodeID(1,2)];
               frame.wakelines{j}.bpanelID = [bpanelID(1),bpanelID(i)];
               frame.wakelines{j}.wpanelID = wpanelID(1);
               nodeID(i,:) = [];
               bpanelID(i,:) = [];
               wpanelID(i,:) = [];
               nodeID(1,:) = [];
               bpanelID(1,:) = [];
               wpanelID(1,:) = [];
               i = 2;
               j = j+1;
           else
               i = i+1;
           end
           if isempty(nodeID)
               break;
           end
        end

        %upperとlowerを分ける
        for i = 1:size(frame.wakelines,2)
            %wakelineのID
            frame.wakelines{i}.ID = frame.surfID(frame.wakelines{i}.wpanelID);
            if acos(dot(frame.n(frame.wakelines{i}.wpanelID,:),frame.n(frame.wakelines{i}.bpanelID(1),:)))<acos(dot(frame.n(frame.wakelines{i}.wpanelID,:),frame.n(frame.wakelines{i}.bpanelID(2),:)))
                frame.wakelines{i}.upperID = frame.wakelines{i}.bpanelID(1);
                frame.wakelines{i}.lowerID = frame.wakelines{i}.bpanelID(2);
            else
                frame.wakelines{i}.upperID = frame.wakelines{i}.bpanelID(2);
                frame.wakelines{i}.lowerID = frame.wakelines{i}.bpanelID(1);
            end       
        end
    end
    %}
    %sparseLLTを使う場合、その設定
    if not(isempty(LLTwakeID))
        iter = 1;
        LLTwakeNo{iter} = [];
        for i = 1:size(LLTwakeID,1)
            if iter == LLTwakeID(i,1)
                LLTwakeNo{iter} = [LLTwakeNo{iter},LLTwakeID(i,2)];
            else
                iter = LLTwakeID(i,1);
                LLTwakeNo{iter} = [];
                LLTwakeNo{iter} = [LLTwakeNo{iter},LLTwakeID(i,2)];
            end
        end
        %wake計算するsurfIDを特定
        wakeCalcID = [];
        for i = 1:size(frame.wakelines,2)
            if not(any(wakeCalcID==frame.surfID(frame.wakelines{i}.upperID,1)))
                wakeCalcID=[wakeCalcID,frame.surfID(frame.wakelines{i}.upperID,1)];
            end
            if not(any(wakeCalcID==frame.surfID(frame.wakelines{i}.lowerID,1)))
                wakeCalcID=[wakeCalcID,frame.surfID(frame.wakelines{i}.lowerID,1)];
            end
        end
        for j = 1:size(LLTwakeNo,2)
            y{j} = []; yVertsID{j} = [];muBodyID{j} = [];
            for i = 1:size(frame.wakelines,2)
                if any(LLTwakeNo{j} == frame.wakelines{i}.ID)
                    y{j} = [y{j},(frame.verts(frame.wakelines{i}.nodeID(1),2)+frame.verts(frame.wakelines{i}.nodeID(2),2))/2];
                    yVertsID{j} = [yVertsID{j},[frame.wakelines{i}.nodeID(1);frame.wakelines{i}.nodeID(2)]];
                    muBodyID{j} =[muBodyID{j}, [frame.tri2ubody(frame.wakelines{i}.upperID);frame.tri2ubody(frame.wakelines{i}.lowerID)]];
                end
            end
            [~,index] = sort(y{j});y{j} = y{j}(index);yVertsID{j} = yVertsID{j}(:,index);muBodyID{j} = muBodyID{j}(:,index);
        end
        frame.wakeCalcID = wakeCalcID;
        frame.yVertsID = yVertsID;frame.muBodyID = muBodyID;
    end
end
end

function [darea, dvolume, n] = vartex_iter(P0,P1,P2)
    %座標を与えてSTLの繰り返し部分を追加書き込みする
    x10 = P1(1)-P0(1);y10 = P1(2)-P0(2);z10 = P1(3)-P0(3);
    x20 = P2(1)-P0(1);y20 = P2(2)-P0(2);z20 = P2(3)-P0(3);
    nx=y10*z20-z10*y20;
    ny=z10*x20-x10*z20;
    nz=x10*y20-y10*x20;
    N = norm([nx,ny,nz]);
    darea = sqrt(nx*nx+ny*ny+nz*nz)/2;
    dvolume = (nx*P0(1)+ny*P0(2)+nz*P0(3))/6;
    n = [nx,ny,nz]./N;
end

function [ID, edgeangle] = neighbourcheck(nodeID,node,frame)
ID = [];
edgeangle = [];
for i = 1:size(frame.tri,1)
    if frame.isBody(i,1)==1
        if not(nodeID == i)
            if any(node(1) == frame.tri(i,:))
                if any(node(2) == frame.tri(i,:))
                    ID = [ID,i];
                    edgeangle = [edgeangle,acos(dot(frame.n(nodeID,:),frame.n(ID(end),:)))];
                end
            end
        end
    end
end
end

function bool = neighbourcheckID(nodeID,node,checkID,frame)
bool = 0;
i = checkID;
if not(nodeID == i)
    if any(node(1) == frame.tri(i,:))
        if any(node(2) == frame.tri(i,:))
            bool = 1;
            return;
        end
    end
end
end

function [ID, edgeangle] = neighbourwakecheck(nodeID,node,frame)
ID = [];
edgeangle = [];
for i = 1:size(frame.tri,1)
    if frame.isWake(i,1)==1
        if not(nodeID == i)
            if any(node(1) == frame.tri(i,:))
                if any(node(2) == frame.tri(i,:))
                    ID = [ID,i];
                    edgeangle = [edgeangle,acos(dot(frame.n(nodeID,:),frame.n(ID(end),:)))];
                end
            end
        end
    end
end
end