%
%%%%%%%%ここから
clear;
%%%%%%%%%
%パイプを等価な剛性を持つspar2本に変換
[tspar,b,Rarea] = InatiaPipe2biRec(80,2.5,766*0.12); %mm単位
%%%%%%%%%
Vinf = 7.2;%飛行速度
Re = Vinf * 0.88 /1.512 *1e5; %機体レイノルズ数
alpha = 0; %迎え角
beta = 0; %横滑り角

[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "FWpeller2.vspgeom", 0); %空力メッシュの作成

wing = UNLSI(p',con',id',wedata,1); %インスタンスの作成
wing = wing.setREFS(25.44,33.75,0.766,[0,0,0]); %基準値をセット
wing = wing.setUNLSISettings("propCalcFlag",1,"laminarRatio",0.4,"propWakeLength",1,"Vinf",Vinf); %設定を変更 アクチュエータディスク（プロペラ計算）の有効化、層流比を0.4に、プロペラwakeの長さをプロペラ直径の4倍に、飛行速度を設定
wing = wing.setProp(1,4,1.6,1); %空力メッシュID4晩を
wing = wing.makeCluster();
wing = wing.makeEquation();
%}
%
[confem,vertsfem,idfem] = readFemMesh('FWpeller2_WingGeom_Struct0.msh'); %FEMメッシュを読み込み
wing = wing.setFemMesh(vertsfem,confem,idfem,[1]); %空力メッシュID1番をFEMメッシュに対応付ける
[wing,weight] = wing.setFemMaterials([1,2,3,4],[0.0005,tspar/1000,0.01,tspar/1000],[5*1e9,12*1e9,5*1e9,12*1e9],[10,1600*Rarea,100,1600*Rarea],[100,100,100,100]);%物性値をセット 肉厚,ヤング率,密度,減衰パラメータ
disp("weight");
disp(weight);
wing = wing.makeFemEquation(); %FEM行列の作成
%%%%%%%%ここまでは一度計算すればスキップできる
%}
%
wing = wing.setPropState(1,0.4,0.8,160); %プロペラの回転状況を指定
wing = wing.solveFlow(alpha,beta,0.001,Re); %パネル法を解く
wing.plotGeometry(1,wing.getCp(0),[-3,1.5]);%圧力係数のプロット
disp(wing.getAERODATA(alpha));
%}

wing2 = wing;%変位後の空力解析用のインスタンスを準備
for iter = 1
    delta = wing.solveFem(wing.getCp(0).*Vinf.^2.*1.225.*0.5,1); %femを解く
    modVerts = wing.calcModifiedVerts(delta{1}); %結果を空力メッシュに投影する
    wing.plotFemMesh(3,delta{1});
    wing2 = wing2.setVerts(modVerts); %変位後の空力メッシュをセット
    wing2 = wing2.makeEquation(); %パネル法行列の作成
    wing2 = wing2.solveFlow(alpha,beta,0.001,Re);%パネル法を解く
    disp(wing2.getAERODATA(alpha));
    wing2.plotGeometry(2,wing2.getCp(0),[-3,1.5]);%圧力係数のプロット
end
%
