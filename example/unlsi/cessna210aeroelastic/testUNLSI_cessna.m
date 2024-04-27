%
%%%%%%%%ここから
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "Cessna-210.vspgeom", 0);%形状の読み込み
[confem,vertsfem,idfem] = readFemMesh('Cessna-210_NormalWing_Struct0.msh');
cessna = UNLSI(p',con',id',wedata,1); %コンストラクタの実行
cessna = cessna.setREFS(175,36.75,4.91,[0,0,0]); %基準面積 基準長の設定
cessna= cessna.setUNLSISettings("Vinf",70.82);
cessna = cessna.makeCluster(); %速度分布を求めるためのパネルクラスターを作成
cessna = cessna.makeEquation(); %パネル法行列の作成
%}
%
cessna = cessna.setFemMesh(vertsfem,confem,idfem,[5]);%空力メッシュ5番とfemメッシュを対応付ける
[cessna,weight] = cessna.setFemMaterials([1,2,3],[0.001,0.003,0.003],[73500000000,73500000000,73500000000],[2700,2700,2700],[100,100,100]);
disp("weight");
disp(weight);
cessna = cessna.makeFemEquation();
%%%%%%%%ここまでは一度計算すればスキップできる
%}
%
alpha = 0;
cessna = cessna.setWakeShape([1,0,0]);
cessna = cessna.solveFlow(alpha,0,0.001,2.3*10^6); %パネル法を解く
cessna.plotGeometry(1,cessna.getCp(alpha),[-2,1]);%圧力係数のプロット
disp(cessna.getAERODATA(alpha))
cessna2 = cessna; %空力弾性計算用に元のインスタンスを取っておく
dt = 0.05; %時間区切り

for i = 1:100
    disp(i);
    tic;
    if i == 1
        [delta,deltadot]  = cessna2.solveAeroelastic([0,dt],[],[],cessna2.getCp(alpha,0,0.001,2.3*10^6).*80.^2.*1.225.*0.5,1); %空力弾性計算　初期値0
    else
        [delta,deltadot]  = cessna2.solveAeroelastic([0,dt],delta,deltadot,cessna2.getCp(alpha,0,0.001,2.3*10^6).*80.^2.*1.225.*0.5,1);%空力弾性計算
    end
    toc;
    modVerts = cessna.calcModifiedVerts(delta{1}); %構造メッシュの変形に従って空力メッシュを変形
    cessna2 = cessna2.setVerts(modVerts); %節点の移動のみ
    cessna2 = cessna2.marchWake(dt,alpha,0,0.001,2.3*10^6);
    cessna2 = cessna2.makeEquation(); %パネル法行列の作成
    cessna2 = cessna2.solveFlow(alpha,0,0.001,2.3*10^6);%パネル法を解く
    disp(cessna2.getAERODATA(alpha,0));
    cessna2.plotGeometry(2,cessna2.getCp(alpha,0,0.001,2.3*10^6),[-3,1.5]);
    cessna2.plotWakeShape(2);
    M(i) = getframe(gcf);
    videoMaker(M,"cessnaCalculating",1/dt)
end
movie(M,1,1/dt);