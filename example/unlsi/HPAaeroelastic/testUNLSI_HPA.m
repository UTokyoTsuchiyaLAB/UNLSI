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
wing = wing.setWakeShape([1,0,0]);
wing = wing.makeCluster();
wing = wing.makeEquation();
%}
%
[confem,vertsfem,idfem] = readFemMesh('FWpeller2_WingGeom_Struct0.msh'); %FEMメッシュを読み込み
wing = wing.setFemMesh(vertsfem,confem,idfem,[1]); %空力メッシュID1番をFEMメッシュに対応付ける
[wing,weight] = wing.setFemMaterials([1,2,3,4],[0.0005,tspar/1000,0.01,tspar/1000],[5*1e9,5*1e9,5*1e9,5*1e9],[10,1600*Rarea,100,1600*Rarea],[10,10,10,10]);%物性値をセット 肉厚,ヤング率,密度,減衰パラメータ
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
%
wing2 = wing; %空力弾性計算用に元のインスタンスを取っておく
dt = 0.05; %時間区切り
v = VideoWriter('hpaElasticAnalysis.mp4', 'MPEG-4');
% 時間区切りからフレームレートの計算と適用
framerate = 1/dt;
v.FrameRate = framerate;
% 保存する動画の画質。数字の大きいほうが高画質.[0~100]
v.Quality = 95;
% ビデオの書き込みを開始
open(v);
for i = 1:100
    disp(i);
    tic;
    if i == 1
        [delta,deltadot]  = wing2.solveAeroelastic([0,dt],[],[],wing2.getCp(alpha,beta,0.001,Re).*Vinf.^2.*1.225.*0.5,1); %空力弾性計算　初期値0
    else
        [delta,deltadot]  = wing2.solveAeroelastic([0,dt],delta,deltadot,wing2.getCp(alpha,beta,0.001,Re).*Vinf.^2.*1.225.*0.5,1);%空力弾性計算
    end
    toc;
    modVerts = wing.calcModifiedVerts(delta{1}); %構造メッシュの変形に従って空力メッシュを変形
    wing2 = wing2.setVerts(modVerts); %節点の移動のみ
    wing2 = wing2.marchWake(dt,alpha,beta,0.001,Re);

    wing2 = wing2.makeEquation(); %パネル法行列の作成
    wing2 = wing2.solveFlow(alpha,0,0.001,Re);%パネル法を解く
    disp(wing2.getAERODATA(alpha,0));
    wing2.plotGeometry(2,wing2.getCp(alpha,0,0.001,Re),[-3,1.5]);
    wing2.plotWakeShape(2);
    M(i) = getframe(gcf); %動画用のフレームを取得
    writeVideo(v, M(i));
end
close(v);
movie(M,1,1/dt);
%}
