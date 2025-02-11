%
%%%%%%%%ここから
clear;
%%%%%%%%%
%パイプを等価な剛性を持つspar2本に変換
[tspar,b,Rarea] = InatiaPipe2biRec(80,1.5,766*0.12); %mm単位
%%%%%%%%%
Vinf = 10;%7.2;
Re = Vinf * 0.375 /1.512 *1e5;%機体レイノルズ数
alpha = 0; %迎え角
beta = 0; %横滑り角

[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom("flora.vspgeom", 0); %空力メッシュの作成

wing = UNLSI(p',con',id',wedata,1,0.005); %インスタンスの作成
wing = wing.setREFS(25.44,33.75,0.766,[0,0,0]); %基準値をセット
wing = wing.setUNLSISettings("laminarRatio",0.4,"Vinf",Vinf); %設定を変更 アクチュエータディスク（プロペラ計算）の有効化、層流比を0.4に、プロペラwakeの長さをプロペラ直径の4倍に、飛行速度を設定
%wing = wing.setWakeShape([1,0,0]);
wing = wing.makeCluster();
wing = wing.makeEquation();
%}
%
[confem,vertsfem,idfem] = readFemMesh('FWpeller2_WingGeom_Struct0.msh'); %FEMメッシュを読み込み
wing = wing.setFemMesh(vertsfem,confem,idfem,[1]); %空力メッシュID1番をFEMメッシュに対応付ける
[wing,weight] = wing.setFemMaterials([1,2,3,4,5,6,7,8],[0.0005,tspar/1000,tspar/1000,0.01,0.01,0.01,0.01,0.01],[5*1e9,9*1e9,9*1e9,5*1e9,5*1e9,5*1e9,5*1e9,5*1e9],[10,1600,1600,100,100,100,100,100],[100,100,100,100,100,100,100,100]);%物性値をセット 肉厚,ヤング率,密度,減衰パラメータ
disp("weight");
disp(weight);
wing = wing.makeFemEquation(); %FEM行列の作成
wing = wing.femModalAnalysis(15);%モード解析用の固有値解析

%%%%%%%%ここまでは一度計算すればスキップできる
%}
%
wing = wing.solveFlow(alpha,beta,0.001,Re); %パネル法を解く
wing.plotGeometry(1,wing.getCp(alpha),[-3,1.5]);%圧力係数のプロット
disp(wing.getAERODATA(alpha));
for i = 1:10
    delta{1} = wing.femSol2Delta(wing.femEigenVec(:,i));
    modVerts = wing.calcModifiedVerts(delta{1}.*10); %構造メッシュの変形に従って空力メッシュを変形
    wing2 = wing.setVerts(modVerts); %節点の移動のみ
    wing2.plotGeometry(i+2,wing.getCp(alpha),[-3,1.5]);
end

%
wing2 = wing; %空力弾性計算用に元のインスタンスを取っておく
dt = 0.05; %時間区切り
for i = 1:10
    disp(i);
    tic;
    if i == 1
        %[delta,deltadot]  = wing2.solveAeroelastic([0,dt],[],[],wing2.getCp(alpha,beta,0.001,Re).*Vinf.^2.*1.225.*0.5,1); %空力弾性計算　初期値0
        [z,zdot,delta,deltadot]  = wing.solveModalAeroelastic([0,dt],[],[],wing2.getCp(alpha).*Vinf.^2.*1.225.*0.5,1); %空力弾性計算　初期値0
    else
        %[delta,deltadot]  = wing2.solveAeroelastic([0,dt],delta,deltadot,wing2.getCp(alpha,beta,0.001,Re).*Vinf.^2.*1.225.*0.5,1);%空力弾性計算
        [z,zdot,delta,deltadot]  = wing.solveModalAeroelastic([0,dt],z,zdot,wing2.getCp(alpha).*Vinf.^2.*1.225.*0.5,1);%空力弾性計算
    end
    disp(full(z{1}));
    toc;
    modVerts = wing.calcModifiedVerts(delta{1}); %構造メッシュの変形に従って空力メッシュを変形
    wing2 = wing.setVerts(modVerts); %節点の移動のみ

    wing2 = wing2.makeEquation(); %パネル法行列の作成
    wing2 = wing2.solveFlow(alpha,0,0.001,Re);%パネル法を解く
    disp(wing2.getAERODATA(alpha));
    wing2.plotGeometry(2,wing2.getCp(alpha),[-3,1.5]);
    %wing2.plotWakeShape(2);
    M(i) = getframe(gcf); %動画用のフレームを取得
    videoMaker(M,"hpaCalculating",1/dt)
end
movie(M,1,1/dt);
%}
%
%空力係数の補間曲面の作成
wing = wing.makeSurrogateModel([-10,-5,0,5,10],[-10,-5,0,5,10],0.0001,300000);
[ppCoef,ppDyn] = wing.getSurrogateModel();
wing.plotSurrogateModel(ppCoef,ppDyn,linspace(-10,10,100),linspace(-10,10,100),0.001,1);
UREF = 7.2;
Mach = 0.001;
REFS = [wing.SREF,wing.BREF,wing.CREF];
save aeroRBF.mat ppCoef ppDyn REFS mass Inatia UREF Mach -mat;
%}