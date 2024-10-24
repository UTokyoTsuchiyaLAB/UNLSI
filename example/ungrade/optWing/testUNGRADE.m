clear;
orgVal = modifyDesFile("init.des");%設計変数の初期値

lb = [-0.4,-0.4, 0.0,-0.1,-0.1,-5,0.15,-10,-0.4,-0.4,0.0,-0.1,-0.1,0.30]; %設計変数の下限値
ub = [ 0.1, 0.1, 0.4, 0.4, 0.4,25,0.45,  0, 0.1, 0.1,0.4, 0.4, 0.4,0.45]; %設計変数の上限値
cmin = [0.10,-6.0]'; %制約条件の下限値
cmax = [0.11,-2.0]'; %制約条件の上限値
geomErr = 0.001;

ungradetest = UNGRADE(@(x)vspMeshGen(x,"optWing2","org.des"),@(x)vspGeomGen(x,"optWing2","org.des"),orgVal,lb,ub,1); %コンストラクタの実行
ungradetest.checkGeomGenWork(0.5); %設計変数が動いているかチェックする
%%%%%%%%%%%%%%%%%%ここで種々の設定をする%%%%%%%%%%%%%%%%%%%
ungradetest = ungradetest.setUNGRADESettings("femCouplingFlag",1,'TrustRegion',0.1,"dynCoefFlag",1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[con,verts,femID] = readFemMesh('optWing2_WingGeom_Struct0.msh');
ungradetest = ungradetest.setFemMesh(verts,con,femID,1);%すべての空力メッシュIDとfemメッシュを関連付ける（第二引数省略）
[ungradetest,weight] = ungradetest.setFemMaterials([1,2],[0.001,0.003],[2.1,2.1].*1e9,[2700,2700]);
ungradetest = ungradetest.setFlowCondition(0,0,0.001,200000,0.5*1.225*20^2);
[Iorg,conorg,ungradetest] = ungradetest.evaluateObjFun(@objFun);%評価関数を計算する
modVerts = ungradetest.calcModifiedVerts(ungradetest.eqnSol.delta{1,1});
ungradetest2 = ungradetest.setVerts(modVerts);
ungradetest2.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %機体形状と圧力係数をプロットする

M = getframe(gcf);
for i = 1:240
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax);%次の設計変数を計算する
    [con,verts,femID] = ungradetest.calcFemMeshfromVariables(nextVar);
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar,0);%次の設計変数を用いて機体形状を更新する
    ungradetest = ungradetest.setFemMesh(verts,con,femID,1);%すべての空力メッシュIDとfemメッシュを関連付ける（第二引数省略）
    [ungradetest,weight] = ungradetest.setFemMaterials([1,2],[0.001,0.003],[2.1,2.1].*1e9,[2700,2700]);
    %%%%%%%%%ここから
    [Lorg,Iorg,conorg,ungradetest] = ungradetest.calcLagrangian(@objFun,cmin,cmax); 
    modVerts = ungradetest.calcModifiedVerts(ungradetest.eqnSol.delta{1,1});
    ungradetest2 = ungradetest.setVerts(modVerts);
    ungradetest2.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]); %↑の結果をプロットする
    M = [M,getframe(gcf)];
    videoMaker(M,"geometry",10);
    ungradetest.plotOptimizationState(2,"optWing2"); %最適化の状況をプロットする
    %%%%%%%%%ここはなくてもよい
end
save optWing2.mat
