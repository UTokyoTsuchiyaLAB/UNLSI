%
%%%%%%%%ここから
%
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "PazyWing.vspgeom", 0); %形状の読み込み
wing = UNLSI(p',con',id',wedata); %コンストラクタの実行
wing = wing.setREFS(80,20,4); %基準面積 基準長の設定
wing = wing.setRotationCenter([0,0,0]); %回転中心の設定
wing = wing.setUNLSISettings("nCalcDivide",2);%パネル法行列の作成における分割数を設定
wing = wing.makeCluster(); %速度分布を求めるためのパネルクラスターを作成
wing = wing.makeEquation(); %パネル法行列の作成
%}
%}
%%%%%%%%ここまでは一度計算すればスキップできる
%
alpha = 1;
Vinf =10;
Re = Vinf * 0.1 / 1.512 * 1e5;
wing = wing.solveFlow(alpha,0,0.001,Re);%パネル法を解く
wing.plotGeometry(1,wing.getCp(alpha,0,0.001,Re),[-3,1.5]);%圧力係数のプロット
disp(wing.getAERODATA(alpha,0));
[con,verts,femID] = readFemMesh('PazyWing.msh');
wing = wing.setFemMesh(verts,con,femID);%すべての空力メッシュIDとfemメッシュを関連付ける（第二引数省略）
[wing,weight] = wing.setFemMaterials([1,2],[0.003,0.0025],[1.31e9*1000,71.7e9*1000],[1000,2810],[1000,1000]);%物性値をセット 肉厚,ヤング率,密度,減衰パラメータ skin rib sparの順
disp("weight");
disp(weight);
wing = wing.makeFemEquation();
%}

%%%%%%%以下空力弾性計算
wing2 = wing;
dt = 0.01;
% VideoWriter オブジェクトを作成
v = VideoWriter('AeroAnalysis2.mp4', 'MPEG-4');
% 時間区切りからフレームレートの計算と適用
framerate = 1.0/dt;
v.FrameRate = framerate;
% 保存する動画の画質。数字の大きいほうが高画質.[0~100]
v.Quality = 100;
% ビデオの書き込みを開始
open(v);


for i = 1:10
    disp(i)
    tic;
    if i == 1
        [delta,deltadot]  = wing.solveAeroelastic([0,dt],[],[],wing2.getCp(alpha,0,0.001,Re).*Vinf.^2.*1.225.*0.5,1);%deltaとdeltadotに空行列を渡すと、初期値0,0からスタート
    else
        [delta,deltadot]  = wing.solveAeroelastic([0,dt],delta,deltadot,wing2.getCp(alpha,0,0.001,Re).*Vinf.^2.*1.225.*0.5,1);%初期値deltaとdeltadotから、tspan間での空力弾性応答を計算
    end
    toc;
    modVerts = wing.calcModifiedVerts(delta{1});
    wing2 = wing2.setVerts(modVerts);
    wing2 = wing2.makeEquation(); %パネル法行列の作成
    wing2 = wing2.solveFlow(alpha,0,0.001,Re);%パネル法を解く
    disp(wing2.getAERODATA(alpha,0));
    wing2.plotGeometry(2,wing2.getCp(alpha,0,0.001,Re),[-3,1.5]);
    M(i) = getframe;
    % フレームを作成
    frame = getframe(gcf);
    % フレームをビデオに書き込み
    writeVideo(v, frame);
end
movie(M,1,1/dt);
close(v);