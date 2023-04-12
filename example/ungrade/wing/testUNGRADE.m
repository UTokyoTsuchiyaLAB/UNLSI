clear;
orgVal = [4,4,4,4,4];
%orgVal = modifyDesFile("org.des","org.des");


ub = [  4,   4,   4,   4,   4];
lb = [  1,   1,   1,   1,   1];
cmin = [0,0,0,0]';
cmax = [10,10,10,10]';
[orgVerts,orgCon,surfID,wakelineID,orgVal] = vspMeshGen(orgVal,"wing","org.des");
ungradetest = UNGRADE(@(x)vspMeshGen(x,"wing","org.des"),@(x)vspGeomGen(x,"wing","org.des"),orgVal,lb,ub,1,0.001,5,0);
ungradetest.checkGeomGenWork(0.5);
ungradetest = ungradetest.setCfParameter(500000,4,0.052*(10^-5),0,1);
ungradetest = ungradetest.setOptions('H0',eye(5));

for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.2,"betaLM",0.5);%設計変数微分の計算方法の指定
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);
    ungradetest.plotOptimizationState(2);
end
for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.1,"betaLM",0.5);%設計変数微分の計算方法の指定
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);
    ungradetest.plotOptimizationState(2);
end
for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.05,"betaLM",0.0);%設計変数微分の計算方法の指定
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);
    ungradetest.plotOptimizationState(2);
end
%}
