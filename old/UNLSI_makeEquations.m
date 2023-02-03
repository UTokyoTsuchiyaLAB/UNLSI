function frame = UNLSI_makeEquations(verts,REFS,isSymm,frame,onlySuperSonic)

%データの格納
frame.verts = verts;
frame.REFS = REFS;
frame.isSymm = isSymm;

%法線ベクトルとパネル面積の計算(再計算)
frame.nVerts = size(frame.verts,1);
frame.areasum = 0;
frame.volume = 0;
frame.pntX = [frame.verts(frame.tri(:,1),1),frame.verts(frame.tri(:,2),1),frame.verts(frame.tri(:,3),1)];
frame.pntY = [frame.verts(frame.tri(:,1),2),frame.verts(frame.tri(:,2),2),frame.verts(frame.tri(:,3),2)];
frame.pntZ = [frame.verts(frame.tri(:,1),3),frame.verts(frame.tri(:,2),3),frame.verts(frame.tri(:,3),3)];
frame.center = [mean(frame.pntX,2),mean(frame.pntY,2),mean(frame.pntZ,2)];
for i = 1:size(frame.tri,1)
    [frame.area(i,1), dvolume, frame.n(i,:)] = vartex_iter(frame.verts(frame.tri(i,1),:),frame.verts(frame.tri(i,2),:),frame.verts(frame.tri(i,3),:));
    if frame.isBody(i,1) == 1
        frame.areasum = frame.areasum + frame.area(i,1);
        frame.volume = frame.volume + dvolume;
    end
end
if frame.isSymm == 1
    frame.areasum = frame.areasum*2;
    frame.volume = frame.volume*2;
end
if nargin<5
%
%パネル微分行列の作成
frame.VxMat = zeros(frame.nPanel,frame.nPanel);
frame.VyMat = zeros(frame.nPanel,frame.nPanel);
frame.VzMat = zeros(frame.nPanel,frame.nPanel);

for i = 1:frame.nPanel
    if frame.isBody(i,1) == 1
        CPmat =frame.center(frame.cluster{i},1:3);
        pnt = frame.center(i,:);
        m = frame.verts(frame.tri(i,1),:)'-pnt(:);
        m = m./norm(m);
        l = cross(m,frame.n(i,:)');
        Minv = [l,m,frame.n(i,:)'];
        lmnMat = (Minv\(CPmat-repmat(pnt,[size(CPmat,1),1]))')';
        bb = [lmnMat(1:end,1),lmnMat(1:end,2),lmnMat(1:end,3),ones(size(lmnMat,1),1)];
        Bmat=pinv(bb);
        Vnmat = Minv(:,[1,2])*[1,0,0,0;0,1,0,0]*Bmat;
        frame.VxMat(i,frame.cluster{i}(1:end)) = Vnmat(1,:);
        frame.VyMat(i,frame.cluster{i}(1:end)) = Vnmat(2,:);
        frame.VzMat(i,frame.cluster{i}(1:end)) = Vnmat(3,:);
    end
end

%sparseLLTを行う場合、LLTmatを作成
if isfield(frame,'yVertsID')
    for i = 1:size(frame.yVertsID,2)
        y = ((verts(frame.yVertsID{i}(1,:),2)+verts(frame.yVertsID{i}(2,:),2))./2)';
        frame.bLLT(i) = (y(end)+mean(diff(y))/2)*2;
        frame.LLTMat{i} = makeLLTMat(y,frame.bLLT(i),frame.isSymm);
    end
end



frame.infA = zeros(frame.nbPanel);
frame.infB = zeros(frame.nbPanel);
POI.X = zeros(frame.nbPanel,1);
POI.Y = zeros(frame.nbPanel,1);
POI.Z = zeros(frame.nbPanel,1);
Amat = zeros(1,frame.nbPanel);
c.X = zeros(1,frame.nbPanel);
c.Y = zeros(1,frame.nbPanel);
c.Z = zeros(1,frame.nbPanel);
n.X = zeros(1,frame.nbPanel);
n.Y = zeros(1,frame.nbPanel);
n.Z = zeros(1,frame.nbPanel);
POI.X(:,1) = frame.center(frame.isBody==1,1);
POI.Y(:,1) = frame.center(frame.isBody==1,2);
POI.Z(:,1) = frame.center(frame.isBody==1,3);
Amat(1,:) = frame.area(frame.isBody==1,1)';
c.X(1,:) = frame.center(frame.isBody==1,1)';
c.Y(1,:) = frame.center(frame.isBody==1,2)';
c.Z(1,:) = frame.center(frame.isBody==1,3)';
n.X(1,:) = frame.n(frame.isBody==1,1)';
n.Y(1,:) = frame.n(frame.isBody==1,2)';
n.Z(1,:) = frame.n(frame.isBody==1,3)';
N1.X(1,:) = frame.pntX(frame.isBody==1,1)';
N1.Y(1,:) = frame.pntY(frame.isBody==1,1)';
N1.Z(1,:) = frame.pntZ(frame.isBody==1,1)';
N2.X(1,:) = frame.pntX(frame.isBody==1,2)';
N2.Y(1,:) = frame.pntY(frame.isBody==1,2)';
N2.Z(1,:) = frame.pntZ(frame.isBody==1,2)';
N3.X(1,:) = frame.pntX(frame.isBody==1,3)';
N3.Y(1,:) = frame.pntY(frame.isBody==1,3)';
N3.Z(1,:) = frame.pntZ(frame.isBody==1,3)';
%行列に拡張tf
POI.X = repmat(POI.X,[1,size(POI.X,1)]);
POI.Y = repmat(POI.Y,[1,size(POI.Y,1)]);
POI.Z = repmat(POI.Z,[1,size(POI.Z,1)]);
c.X = repmat(c.X,[frame.nbPanel,1]);
c.Y = repmat(c.Y,[frame.nbPanel,1]);
c.Z = repmat(c.Z,[frame.nbPanel,1]);
n.X = repmat(n.X,[frame.nbPanel,1]);
n.Y = repmat(n.Y,[frame.nbPanel,1]);
n.Z = repmat(n.Z,[frame.nbPanel,1]);
N1.X = repmat(N1.X,[frame.nbPanel,1]);
N1.Y = repmat(N1.Y,[frame.nbPanel,1]);
N1.Z = repmat(N1.Z,[frame.nbPanel,1]);
N2.X = repmat(N2.X,[frame.nbPanel,1]);
N2.Y = repmat(N2.Y,[frame.nbPanel,1]);
N2.Z = repmat(N2.Z,[frame.nbPanel,1]);
N3.X = repmat(N3.X,[frame.nbPanel,1]);
N3.Y = repmat(N3.Y,[frame.nbPanel,1]);
N3.Z = repmat(N3.Z,[frame.nbPanel,1]);

n12.X = (N1.X+N2.X)./2;
n12.Y = (N1.Y+N2.Y)./2;
n12.Z = (N1.Z+N2.Z)./2;
m = getUnitVector(c,n12);
l = matrix_cross(m,n);
pjk.X = POI.X-c.X;
pjk.Y = POI.Y-c.Y;
pjk.Z = POI.Z-c.Z;
PN = matrix_dot(pjk,n);

%1回目
a.X = POI.X-N1.X;
a.Y = POI.Y-N1.Y;
a.Z = POI.Z-N1.Z;
b.X = POI.X-N2.X;
b.Y = POI.Y-N2.Y;
b.Z = POI.Z-N2.Z;
s.X = N2.X-N1.X;
s.Y = N2.Y-N1.Y;
s.Z = N2.Z-N1.Z;
Al = matrix_dot(n,matrix_cross(s,a));

PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
PB = PA-Al.*matrix_dot(s,m);
num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
phiV = atan2(num,denom);
A = matrix_norm(a);
B = matrix_norm(b);
S = matrix_norm(s);
GL = 1./S.*log(((A+B+S)./(A+B-S)));
srcV = Al.*GL-PN.*phiV;
VortexA = phiV;
VortexB = srcV;
%2回目
a.X = POI.X-N2.X;
a.Y = POI.Y-N2.Y;
a.Z = POI.Z-N2.Z;
b.X = POI.X-N3.X;
b.Y = POI.Y-N3.Y;
b.Z = POI.Z-N3.Z;
s.X = N3.X-N2.X;
s.Y = N3.Y-N2.Y;
s.Z = N3.Z-N2.Z;

Al = matrix_dot(n,matrix_cross(s,a));

PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
PB = PA-Al.*matrix_dot(s,m);
num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
phiV = atan2(num,denom);
A = matrix_norm(a);
B = matrix_norm(b);
S = matrix_norm(s);
GL = 1./S.*log(((A+B+S)./(A+B-S)));
srcV = Al.*GL-PN.*phiV;
VortexA = VortexA+phiV;
VortexB = VortexB+srcV;
%3回目
a.X = POI.X-N3.X;
a.Y = POI.Y-N3.Y;
a.Z = POI.Z-N3.Z;
b.X = POI.X-N1.X;
b.Y = POI.Y-N1.Y;
b.Z = POI.Z-N1.Z;
s.X = N1.X-N3.X;
s.Y = N1.Y-N3.Y;
s.Z = N1.Z-N3.Z;

Al = matrix_dot(n,matrix_cross(s,a));

PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
PB = PA-Al.*matrix_dot(s,m);
num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
phiV = atan2(num,denom);
A = matrix_norm(a);
B = matrix_norm(b);
S = matrix_norm(s);
GL = 1./S.*log(((A+B+S)./(A+B-S)));
srcV = Al.*GL-PN.*phiV;
VortexA = VortexA+phiV;
VortexB = VortexB+srcV;

if frame.isSymm == 1
    POIm = POI;
    POIm.Y = -POI.Y;
    n12.X = (N1.X+N2.X)./2;
    n12.Y = (N1.Y+N2.Y)./2;
    n12.Z = (N1.Z+N2.Z)./2;
    m = getUnitVector(c,n12);
    l = matrix_cross(m,n);
    pjk.X = POIm.X-c.X;
    pjk.Y = POIm.Y-c.Y;
    pjk.Z = POIm.Z-c.Z;
    PN = matrix_dot(pjk,n);

    %1回目
    a.X = POIm.X-N1.X;
    a.Y = POIm.Y-N1.Y;
    a.Z = POIm.Z-N1.Z;
    b.X = POIm.X-N2.X;
    b.Y = POIm.Y-N2.Y;
    b.Z = POIm.Z-N2.Z;
    s.X = N2.X-N1.X;
    s.Y = N2.Y-N1.Y;
    s.Z = N2.Z-N1.Z;
    Al = matrix_dot(n,matrix_cross(s,a));

    PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
    PB = PA-Al.*matrix_dot(s,m);
    num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
    denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
    phiV = atan2(num,denom);
    A = matrix_norm(a);
    B = matrix_norm(b);
    S = matrix_norm(s);
    GL = 1./S.*log(((A+B+S)./(A+B-S)));
    srcV = Al.*GL-PN.*phiV;
    VortexA = VortexA+phiV;
    VortexB = VortexB+srcV;
    %2回目
    a.X = POIm.X-N2.X;
    a.Y = POIm.Y-N2.Y;
    a.Z = POIm.Z-N2.Z;
    b.X = POIm.X-N3.X;
    b.Y = POIm.Y-N3.Y;
    b.Z = POIm.Z-N3.Z;
    s.X = N3.X-N2.X;
    s.Y = N3.Y-N2.Y;
    s.Z = N3.Z-N2.Z;

    Al = matrix_dot(n,matrix_cross(s,a));

    PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
    PB = PA-Al.*matrix_dot(s,m);
    num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
    denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
    phiV = atan2(num,denom);
    A = matrix_norm(a);
    B = matrix_norm(b);
    S = matrix_norm(s);
    GL = 1./S.*log(((A+B+S)./(A+B-S)));
    srcV = Al.*GL-PN.*phiV;
    VortexA = VortexA+phiV;
    VortexB = VortexB+srcV;
    %3回目
    a.X = POIm.X-N3.X;
    a.Y = POIm.Y-N3.Y;
    a.Z = POIm.Z-N3.Z;
    b.X = POIm.X-N1.X;
    b.Y = POIm.Y-N1.Y;
    b.Z = POIm.Z-N1.Z;
    s.X = N1.X-N3.X;
    s.Y = N1.Y-N3.Y;
    s.Z = N1.Z-N3.Z;

    Al = matrix_dot(n,matrix_cross(s,a));

    PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
    PB = PA-Al.*matrix_dot(s,m);
    num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
    denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
    phiV = atan2(num,denom);
    A = matrix_norm(a);
    B = matrix_norm(b);
    S = matrix_norm(s);
    GL = 1./S.*log(((A+B+S)./(A+B-S)));
    srcV = Al.*GL-PN.*phiV;
    VortexA = VortexA+phiV;
    VortexB = VortexB+srcV;
end
frame.infA = VortexA; 
frame.infB = VortexB; 
eyeR = eye(size(n.X));
frame.infA(eyeR==1) = -2.*pi;

%wakeの影響
xwake = 100.*frame.REFS(2,3);
aw = 0;
for j = 1:numel(frame.wakelines)
    POIw.X(1:frame.nbPanel,1) = frame.center(frame.isBody==1,1);
    POIw.Y(1:frame.nbPanel,1) = frame.center(frame.isBody==1,2);
    POIw.Z(1:frame.nbPanel,1) = frame.center(frame.isBody==1,3);
    wakepos(1,:) = frame.verts(frame.wakelines{j}.nodeID(1),:);
    wakepos(2,:) = frame.verts(frame.wakelines{j}.nodeID(2),:);
    wakepos(3,:) = frame.verts(frame.wakelines{j}.nodeID(2),:)+[xwake*cosd(aw),0,xwake*sind(aw)];
    wakepos(4,:) = frame.verts(frame.wakelines{j}.nodeID(1),:)+[xwake*cosd(aw),0,xwake*sind(aw)];
    ulvec = frame.center(frame.wakelines{j}.lowerID,:)-frame.center(frame.wakelines{j}.upperID,:);
    [b1, b2, nbuff] = vartex_iter(wakepos(1,:),wakepos(2,:),wakepos(3,:));
    uldist = dot(ulvec,nbuff);
    if uldist >= 0
        Nw1.X = repmat(wakepos(1,1),[frame.nbPanel,1]);
        Nw1.Y = repmat(wakepos(1,2),[frame.nbPanel,1]);
        Nw1.Z = repmat(wakepos(1,3),[frame.nbPanel,1]);
        Nw2.X = repmat(wakepos(2,1),[frame.nbPanel,1]);
        Nw2.Y = repmat(wakepos(2,2),[frame.nbPanel,1]);
        Nw2.Z = repmat(wakepos(2,3),[frame.nbPanel,1]);
        Nw3.X = repmat(wakepos(3,1),[frame.nbPanel,1]);
        Nw3.Y = repmat(wakepos(3,2),[frame.nbPanel,1]);
        Nw3.Z = repmat(wakepos(3,3),[frame.nbPanel,1]);
        Nw4.X = repmat(wakepos(4,1),[frame.nbPanel,1]);
        Nw4.Y = repmat(wakepos(4,2),[frame.nbPanel,1]);
        Nw4.Z = repmat(wakepos(4,3),[frame.nbPanel,1]);
    else
        Nw1.X = repmat(wakepos(2,1),[frame.nbPanel,1]);
        Nw1.Y = repmat(wakepos(2,2),[frame.nbPanel,1]);
        Nw1.Z = repmat(wakepos(2,3),[frame.nbPanel,1]);
        Nw2.X = repmat(wakepos(1,1),[frame.nbPanel,1]);
        Nw2.Y = repmat(wakepos(1,2),[frame.nbPanel,1]);
        Nw2.Z = repmat(wakepos(1,3),[frame.nbPanel,1]);
        Nw3.X = repmat(wakepos(4,1),[frame.nbPanel,1]);
        Nw3.Y = repmat(wakepos(4,2),[frame.nbPanel,1]);
        Nw3.Z = repmat(wakepos(4,3),[frame.nbPanel,1]);
        Nw4.X = repmat(wakepos(3,1),[frame.nbPanel,1]);
        Nw4.Y = repmat(wakepos(3,2),[frame.nbPanel,1]);
        Nw4.Z = repmat(wakepos(3,3),[frame.nbPanel,1]);
    end

    [b1, b2, nw] = vartex_iter([Nw1.X(1),Nw1.Y(1),Nw1.Z(1)],[Nw2.X(1),Nw2.Y(1),Nw2.Z(1)],[Nw3.X(1),Nw3.Y(1),Nw3.Z(1)]);
    cw = mean([[Nw1.X(1),Nw1.Y(1),Nw1.Z(1)];[Nw2.X(1),Nw2.Y(1),Nw2.Z(1)];[Nw3.X(1),Nw3.Y(1),Nw3.Z(1)];[Nw4.X(1),Nw4.Y(1),Nw4.Z(1)]],1);

    n.X = repmat(nw(1),[frame.nbPanel,1]);
    n.Y = repmat(nw(2),[frame.nbPanel,1]);
    n.Z = repmat(nw(3),[frame.nbPanel,1]);
    Amat = repmat(b1*2,[frame.nbPanel,1]);
    c.X = repmat(cw(1),[frame.nbPanel,1]);
    c.Y = repmat(cw(2),[frame.nbPanel,1]);
    c.Z = repmat(cw(3),[frame.nbPanel,1]);
    n12.X = (Nw3.X+Nw4.X)./2;
    n12.Y = (Nw3.Y+Nw4.Y)./2;
    n12.Z = (Nw3.Z+Nw4.Z)./2;
    m = getUnitVector(c,n12);
    l = matrix_cross(m,n);
    pjk.X = POIw.X-c.X;
    pjk.Y = POIw.Y-c.Y;
    pjk.Z = POIw.Z-c.Z;
    PN = matrix_dot(pjk,n);
    %VortexA = PN.*Amat./matrix_norm(pjk).^3;
    %
   %1回目
    a.X = POIw.X-Nw1.X;
    a.Y = POIw.Y-Nw1.Y;
    a.Z = POIw.Z-Nw1.Z;
    b.X = POIw.X-Nw2.X;
    b.Y = POIw.Y-Nw2.Y;
    b.Z = POIw.Z-Nw2.Z;
    s.X = Nw2.X-Nw1.X;
    s.Y = Nw2.Y-Nw1.Y;
    s.Z = Nw2.Z-Nw1.Z;
    Al = matrix_dot(n,matrix_cross(s,a));

    PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
    PB = PA-Al.*matrix_dot(s,m);
    num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
    denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
    %sv = sigmoid(denom);
    %phiV = atan2(num,denom).*sv+(1-sv).*pi./2;
    phiV = atan2(num,denom);
    VortexA = phiV;
    %2回目
    a.X = POIw.X-Nw2.X;
    a.Y = POIw.Y-Nw2.Y;
    a.Z = POIw.Z-Nw2.Z;
    b.X = POIw.X-Nw3.X;
    b.Y = POIw.Y-Nw3.Y;
    b.Z = POIw.Z-Nw3.Z;
    s.X = Nw3.X-Nw2.X;
    s.Y = Nw3.Y-Nw2.Y;
    s.Z = Nw3.Z-Nw2.Z;
    Al = matrix_dot(n,matrix_cross(s,a));

    PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
    PB = PA-Al.*matrix_dot(s,m);
    num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
    denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
    %sv = sigmoid(denom);
    %phiV = atan2(num,denom).*sv+(1-sv).*pi./2;
    phiV = atan2(num,denom);
    VortexA = VortexA+phiV;
    %3回目
    a.X = POIw.X-Nw3.X;
    a.Y = POIw.Y-Nw3.Y;
    a.Z = POIw.Z-Nw3.Z;
    b.X = POIw.X-Nw4.X;
    b.Y = POIw.Y-Nw4.Y;
    b.Z = POIw.Z-Nw4.Z;
    s.X = Nw4.X-Nw3.X;
    s.Y = Nw4.Y-Nw3.Y;
    s.Z = Nw4.Z-Nw3.Z;
    Al = matrix_dot(n,matrix_cross(s,a));
    PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
    PB = PA-Al.*matrix_dot(s,m);
    num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
    denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
    %sv = sigmoid(denom);
    %phiV = atan2(num,denom).*sv+(1-sv).*pi./2;
    phiV = atan2(num,denom);
    VortexA = VortexA+phiV;
    %4回目
    a.X = POIw.X-Nw4.X;
    a.Y = POIw.Y-Nw4.Y;
    a.Z = POIw.Z-Nw4.Z;
    b.X = POIw.X-Nw1.X;
    b.Y = POIw.Y-Nw1.Y;
    b.Z = POIw.Z-Nw1.Z;
    s.X = Nw1.X-Nw4.X;
    s.Y = Nw1.Y-Nw4.Y;
    s.Z = Nw1.Z-Nw4.Z;
    Al = matrix_dot(n,matrix_cross(s,a));
    PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
    PB = PA-Al.*matrix_dot(s,m);
    num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
    denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
    %sv = sigmoid(denom);
    %phiV = atan2(num,denom).*sv+(1-sv).*pi./2;
    phiV = atan2(num,denom);
    VortexA = VortexA+phiV;
    %}
    influence = VortexA;
    interpID(1) = frame.tri2ubody(frame.wakelines{j}.upperID);
    interpID(2) = frame.tri2ubody(frame.wakelines{j}.lowerID);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    frame.infA(:,interpID(1)) = frame.infA(:,interpID(1)) - influence;
    frame.infA(:,interpID(2)) = frame.infA(:,interpID(2)) + influence;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if frame.isSymm == 1
        POIw.Y = -POIw.Y;
        
        n12.X = (Nw3.X+Nw4.X)./2;
        n12.Y = (Nw3.Y+Nw4.Y)./2;
        n12.Z = (Nw3.Z+Nw4.Z)./2;
        m = getUnitVector(c,n12);
        l = matrix_cross(m,n);
        pjk.X = POIw.X-c.X;
        pjk.Y = POIw.Y-c.Y;
        pjk.Z = POIw.Z-c.Z;
        PN = matrix_dot(pjk,n);
        %VortexA = PN.*Amat./matrix_norm(pjk).^3;
        %
       %1回目
        a.X = POIw.X-Nw1.X;
        a.Y = POIw.Y-Nw1.Y;
        a.Z = POIw.Z-Nw1.Z;
        b.X = POIw.X-Nw2.X;
        b.Y = POIw.Y-Nw2.Y;
        b.Z = POIw.Z-Nw2.Z;
        s.X = Nw2.X-Nw1.X;
        s.Y = Nw2.Y-Nw1.Y;
        s.Z = Nw2.Z-Nw1.Z;
        Al = matrix_dot(n,matrix_cross(s,a));

        PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
        PB = PA-Al.*matrix_dot(s,m);
        num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
        denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
        %sv = sigmoid(denom);
        %phiV = atan2(num,denom).*sv+(1-sv).*pi./2;
        phiV = atan2(num,denom);
        VortexA = phiV;
        %2回目
        a.X = POIw.X-Nw2.X;
        a.Y = POIw.Y-Nw2.Y;
        a.Z = POIw.Z-Nw2.Z;
        b.X = POIw.X-Nw3.X;
        b.Y = POIw.Y-Nw3.Y;
        b.Z = POIw.Z-Nw3.Z;
        s.X = Nw3.X-Nw2.X;
        s.Y = Nw3.Y-Nw2.Y;
        s.Z = Nw3.Z-Nw2.Z;
        Al = matrix_dot(n,matrix_cross(s,a));

        PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
        PB = PA-Al.*matrix_dot(s,m);
        num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
        denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
        %sv = sigmoid(denom);
        %phiV = atan2(num,denom).*sv+(1-sv).*pi./2;
        phiV = atan2(num,denom);
        VortexA = VortexA+phiV;
        %3回目
        a.X = POIw.X-Nw3.X;
        a.Y = POIw.Y-Nw3.Y;
        a.Z = POIw.Z-Nw3.Z;
        b.X = POIw.X-Nw4.X;
        b.Y = POIw.Y-Nw4.Y;
        b.Z = POIw.Z-Nw4.Z;
        s.X = Nw4.X-Nw3.X;
        s.Y = Nw4.Y-Nw3.Y;
        s.Z = Nw4.Z-Nw3.Z;
        Al = matrix_dot(n,matrix_cross(s,a));
        PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
        PB = PA-Al.*matrix_dot(s,m);
        num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
        denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
        %sv = sigmoid(denom);
        %phiV = atan2(num,denom).*sv+(1-sv).*pi./2;
        phiV = atan2(num,denom);
        VortexA = VortexA+phiV;
        %4回目
        a.X = POIw.X-Nw4.X;
        a.Y = POIw.Y-Nw4.Y;
        a.Z = POIw.Z-Nw4.Z;
        b.X = POIw.X-Nw1.X;
        b.Y = POIw.Y-Nw1.Y;
        b.Z = POIw.Z-Nw1.Z;
        s.X = Nw1.X-Nw4.X;
        s.Y = Nw1.Y-Nw4.Y;
        s.Z = Nw1.Z-Nw4.Z;
        Al = matrix_dot(n,matrix_cross(s,a));
        PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
        PB = PA-Al.*matrix_dot(s,m);
        num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
        denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
        %sv = sigmoid(denom);
        %phiV = atan2(num,denom).*sv+(1-sv).*pi./2;
        phiV = atan2(num,denom);
        VortexA = VortexA+phiV;
        %}
        influence = VortexA;
        interpID(1) = frame.tri2ubody(frame.wakelines{j}.upperID);
        interpID(2) = frame.tri2ubody(frame.wakelines{j}.lowerID);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    frame.infA(:,interpID(1)) = frame.infA(:,interpID(1)) - influence;
    frame.infA(:,interpID(2)) = frame.infA(:,interpID(2)) + influence;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   end
end
%}
if cond(frame.infA)>1000
    error('bad condition of influence matrix A');
end
end
end


function m = getUnitVector(p1,p2)
m.X = p2.X-p1.X;
m.Y = p2.Y-p1.Y;
m.Z = p2.Z-p1.Z;
no = sqrt(m.X.^2+m.Y.^2+m.Z.^2);
index = find(no);
m.X(index) = m.X(index) ./ no(index);
m.Y(index) = m.Y(index) ./ no(index);
m.Z(index) = m.Z(index) ./ no(index);
end

function no = matrix_norm(m)
no = sqrt(m.X.^2+m.Y.^2+m.Z.^2);
end

function md = matrix_dot(A,B)
md = A.X.*B.X+A.Y.*B.Y+A.Z.*B.Z;
end

function mc = matrix_cross(A,B)
mc.X = A.Y.*B.Z-A.Z.*B.Y;
mc.Y = A.Z.*B.X-A.X.*B.Z;
mc.Z = A.X.*B.Y-A.Y.*B.X;
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

function R = RbfMatMake(xd,rbfMode,r0)
    if rbfMode ==4
        phi = @(r,r0)phi4(r,r0);
    elseif rbfMode == 1
        phi = @(r,r0)phi1(r,r0);
    elseif rbfMode == 2
        phi = @(r,r0)phi2(r,r0);
    else
        error('未実装')
    end
    
    X = xd';
    H = sum(xd.^2,2);
    H = repmat(H,[1,size(X,2)]);
    r = sqrt(H'-2.*X'*X+H);
    R = phi(r,r0);

end

function  phi = phi1(r,r0)
    phi = sqrt(r.*r+r0.*r0);
end

function  phi = phi2(r,r0)
    phi = 1./sqrt(r.*r+r0.*r0);
end

function phi = phi4(r,r0)
    phi = exp(-0.5.*r.^2/r0.^2);
end

function rdot = RBFdot(pnt,xd,rbfMode,r0)
rdot = zeros(numel(pnt),size(xd,1));
for i = 1:size(xd,1)
    if rbfMode ==1
        rdot(:,i) = phi1_deriv(pnt(:),xd(i,:)',r0);
    elseif rbfMode == 4
        rdot(:,i) = phi4_deriv(pnt(:),xd(i,:)',r0);
    end
end
end

function deriv = phi1_deriv(xi,xj,r0)
    r = norm(xi-xj);
    deriv = (xi-xj)./sqrt(r*r+r0*r0);
end

function deriv = phi4_deriv(xi,xj,r0)
    r = norm(xi-xj);
    deriv = -(xi-xj)./r0^2.*exp(-0.5.*r.^2/r0.^2);
end

function res = sigmoid(x)
    a = 1000;
    res=1./(1+exp(-a.*x));
end

function optMat = makeLLTMat(y,b,isSymm)
    theta = acos(y./(b/2));
    N = size(theta,2)*2;
    optMat = zeros(N+size(theta,2));
    if isSymm == 1
        for i = 1:N
            optMat(i,i) = pi/2*(2*i-1)^4;
            optMat(i,N+1:N+size(theta,2)) = sin((2*i-1).*theta);
        end
        for j = 1:size(theta,2)
            optMat(N+j,1:N) = sin((2*(1:N)'-1)*theta(j));
        end
    else
        for i = 1:N
            optMat(i,i) = pi/2*i^4;
            optMat(i,N+1:N+size(theta,2)) = sin(i.*theta);
        end
        for j = 1:size(theta,2)
            optMat(N+j,1:N) = sin((1:N)'*theta(j));
        end
    end
end