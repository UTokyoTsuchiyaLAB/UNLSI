clear;
orgVal = [4,4,4,4,4];
%orgVal = modifyDesFile("org.des","org.des");

Machrange = [0.001,0.001,0.001];
alpharange = [2,4,6];
betarange = [0,0,0];

ub = [  4,   4,   4,   4,   4];
lb = [  1,   1,   1,   1,   1];
cmin = [0,0,0,0]';
cmax = [10,10,10,10]';

ungradetest = UNGRADE(@(x)vspMeshGen(x,"wing","org.des"),@(x)vspGeomGen(x,"wing","org.des"),orgVal,lb,ub,1,Machrange,alpharange,betarange);
ungradetest.checkGeomGenWork(0.5);
ungradetest = ungradetest.setCfParameter(500000,4,0.052*(10^-5),0,1);
ungradetest = ungradetest.setOptions('n_divide',3);
ungradetest = ungradetest.solveAnalysis(1,alpharange(1),betarange(1));
ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);

for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.2,"betaLM",0.5);%設計変数微分の計算方法の指定
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);
    ungradetest = ungradetest.solveAnalysis(1,alpharange(1),betarange(1));
    ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end
for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.1,"betaLM",0.5);%設計変数微分の計算方法の指定
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);
    ungradetest = ungradetest.solveAnalysis(1,alpharange(1),betarange(1));
    ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end
for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.05,"betaLM",0.0);%設計変数微分の計算方法の指定
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);
    ungradetest = ungradetest.solveAnalysis(1,alpharange(1),betarange(1));
    ungradetest.plotGeometry(1,ungradetest.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end
%}
