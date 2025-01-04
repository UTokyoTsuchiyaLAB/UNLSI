clear;
orgVal = [3,3,3,3,3,3]; %設計変数の初期値
%[~,x] = modifyDesFile("org.des","org.des");

Machrange = [0.001]; %解析するマッハ数
alpharange = [5]; %解析する迎角
betarange = [0]; %解析する横滑り角
geomErr = 0.005;

lb = [0.5,     0.5,   0.5,0.5,0.5,0.5]; %設計変数の下限値
ub = [5,     5,   5,5 5,5]; %設計変数の上限値
cmin = [2.0, 0, 0, 0, 0, 0,  0]';%制約条件の下限値
cmax = [2.5,10,10,10,10,10,0.50]';%制約条件の上限値

ungradetest = UNGRADE(@(x)vspMeshGen(x,"wing","org.des"),@(x)vspGeomGen(x,"wing","org.des"),orgVal,lb,ub,1);%コンストラクタの実行
ungradetest.checkGeomGenWork(0.5);%設計変数が動いているかチェックする
%%%%%%%%%%%%%%%%%%ここで種々の設定をする%%%%%%%%%%%%%%%%%%%
ungradetest = ungradetest.setUNGRADESettings('femCouplingFlag',1);%弱連成最適化
[con,verts,femID] = readFemMesh('wing_WingGeom_Struct0.msh');
ungradetest = ungradetest.setFemMesh(verts,con,femID);%すべての空力メッシュIDとfemメッシュを関連付ける（第二引数省略）
[ungradetest,weight] = ungradetest.setFemMaterials([1,2,3],[0.002,0.003,0.01],[73500000000,73500000000,73500000000],[2700,2700,2700]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ungradetest = ungradetest.setFlowCondition(alpharange,betarange,Machrange,500000,0.5*1.225*150^2);
[Iorg,conorg,ungradetest] = ungradetest.evaluateObjFun(@objFun); %評価関数を計算する
modVerts = ungradetest.calcModifiedVerts(ungradetest.eqnSol.delta{1,1});
ungradetest2 = ungradetest.setVerts(modVerts);
ungradetest2.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %機体形状と圧力係数をプロットする


%%%%%%%%%%%%%%%%%初級者向け -- とりあえず実行%%%%%%%%%%%%%%%%%%%%%%%%%
%
M = getframe(gcf);
for i = 1:200
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax);%次の設計変数を計算する
    [con,verts,femID] = ungradetest.calcFemMeshfromVariables(nextVar);
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar,0);%次の設計変数を用いて機体形状を更新する
    ungradetest = ungradetest.setFemMesh(verts,con,femID);%すべての空力メッシュIDとfemメッシュを関連付ける（第二引数省略）
    [ungradetest,weight] = ungradetest.setFemMaterials([1,2,3],[0.002,0.003,0.01],[73500000000,73500000000,73500000000],[2700,2700,2700]);
    [Lorg,Iorg,conorg,ungradetest] = ungradetest.calcLagrangian(@objFun,cmin,cmax); 
    modVerts = ungradetest.calcModifiedVerts(ungradetest.eqnSol.delta{1,1});
    ungradetest2 = ungradetest.setVerts(modVerts);
    ungradetest2.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %↑の結果をプロットする
    M = [M,getframe(gcf)];
    videoMaker(M,"geometry",10);
    ungradetest.plotOptimizationState(2,"wing"); %最適化の状況をプロットする
end
%}


