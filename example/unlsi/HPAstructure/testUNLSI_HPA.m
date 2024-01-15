%
%%%%%%%%ここから
clear;

%%%%%%%%%
%パイプを等価な剛性を持つspar2本に変換
[t,b] = InatiaPipe2biRec(80,1,80); %mm単位
%%%%%%%%%

[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "FWpeller2.vspgeom", 0);

wing = UNLSI(p',con',id',wedata,1);
wing = wing.setREFS(25.44,33.75,0.766);
wing = wing.setRotationCenter([0,0,0]);
wing = wing.setUNLSISettings("propCalcFlag",1,"laminarRatio",0.4,"propWakeLength",1,"Vinf",7.2);
wing = wing.setProp(1,4,1.6,1);
wing = wing.makeCluster();
wing = wing.makeEquation();
%}
%
[confem,vertsfem,idfem] = readFemMesh('FWpeller2_WingGeom_Struct0.msh');
wing = wing.setFemMesh(vertsfem,confem,idfem,[1]);
wing = wing.setFemMaterials([1,2,3,4],[0.0001,0.0011,0.0011,0.01],[5000000000,215744100000,215744100000,11500000]);
wing = wing.makeFemEquation();
%%%%%%%%ここまでは一度計算すればスキップできる
%}

%
wing = wing.setPropState(1,0.4,0.8,160); %プロペラの回転状況を指定
wing = wing.setDeflAngle(2,[0,1,0],0);%6番のIDを回転軸[0,1,0]、角度0degに設定（動翼として登録）
wing = wing.setDeflAngle(3,[0,0,1],0);%7番のIDを回転軸[0,1,0]、角度0degに設定（動翼として登録）
wing = wing.setDeflGroup(1,"elev",2,1); % ID2番をgain[1,1](軸が同じなら同じ方向）に動かすことをelevatorとして登録(動翼ID1番） 
wing = wing.setDeflGroup(2,"rud",3,1); % ID3番をgain[1,1](軸が同じなら同じ方向）に動かすことをrudderとして登録(動翼ID2番)
alpha = [-3,0,3];
beta = [0,0,0];

wing = wing.solveFlow(alpha,beta,0.001,300000); %パネル法を解く
wing.plotGeometry(1,wing.getCp(0),[-3,1.5]);%圧力係数のプロット
disp(wing.getAERODATA(alpha));
wing2 = wing;%変位後の空力解析用のインスタンスを準備
for iter = 1
    delta = wing.solveFem(wing2.getCp(0).*7.2.^2.*1.225.*0.5);
    modVerts = wing.calcModifiedVerts(delta{1});
    wing2 = wing2.setVerts(modVerts);
    wing2 = wing2.makeEquation(); %パネル法行列の作成
    wing2 = wing2.solveFlow(alpha,beta,0.001,300000);%パネル法を解く
    disp(wing2.getAERODATA(alpha));
    wing2.plotGeometry(2,wing2.getCp(0),[-3,1.5]);%圧力係数のプロット
end
