%
%%%%%%%%ここから
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "Cessna-210.vspgeom", 0);%形状の読み込み
[confem,vertsfem,idfem] = readFemMesh('Cessna-210_NormalWing_Struct0.msh');
cessna = UNLSI(p',con',id',wedata,1); %コンストラクタの実行
cessna = cessna.setREFS(175,36.75,4.91,[0,0,0]); %基準面積 基準長の設定
cessna = cessna.makeCluster(); %速度分布を求めるためのパネルクラスターを作成
cessna = cessna.makeEquation(); %パネル法行列の作成
%}
%
cessna = cessna.setFemMesh(vertsfem,confem,idfem,[5]);%空力メッシュ5番とfemメッシュを対応付ける
[cessna,weight] = cessna.setFemMaterials([1,2,3],[0.001,0.003,0.003],[73500000000,73500000000,73500000000],[2700,2700,2700]);
disp("weight");
disp(weight);
cessna = cessna.makeFemEquation();
%%%%%%%%ここまでは一度計算すればスキップできる
%}
%
alpha = 0;
cessna = cessna.solveFlow(alpha,0,0.001,2.3*10^6); %パネル法を解く
cessna.plotGeometry(1,cessna.getCp(alpha),[-2,1]);%圧力係数のプロット
disp(cessna.getAERODATA(alpha))
cessna2 = cessna;%変位後の空力解析用のインスタンスを準備
for iter = 1
    delta = cessna.solveFem(cessna2.getCp(alpha).*80.^2.*1.225.*0.5);
    modVerts = cessna.calcModifiedVerts(delta{1});
    cessna2 = cessna2.setVerts(modVerts);
    cessna2 = cessna2.makeEquation(); %パネル法行列の作成
    cessna2 = cessna2.solveFlow(alpha,0,0.001,2.3*10^6);%パネル法を解く
    disp(cessna2.getAERODATA(alpha));
    cessna2.plotGeometry(2,cessna2.getCp(alpha),[-2,1]);%圧力係数のプロット
end