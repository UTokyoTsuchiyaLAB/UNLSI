clear;
orgVal = [0.04,0.4,0.12,0,0.04,0.4,0.12,5.5,0,4,2];

lb = [-0.05,0.2,0.1,0,-0.05,0.2,0.1,4.0,0,2,0];
ub = [ 0.05,0.5,0.14,6,0.05,0.5,0.14,8,20,8,6];
cmin = [40/100, 0.2 -0.3]';
cmax = [45/100, 1.0  0.3]';

ungradetest = UNGRADE(@(x)vspMeshGen(x,"Cessna-210","org.des"),@(x)vspGeomGen(x,"Cessna-210","org.des"),orgVal,lb,ub,1,[0.001,0.001,0.001],[0,5,10],[0,0,0]);
ungradetest.checkGeomGenWork(0.5);
ungradetest = ungradetest.setCfParameter(500000,4,0.052*(10^-5),0,1);
ungradetest = ungradetest.solveAnalysis(1,0,0);
ungradetest.plotGeometry(1,obj.Cp{1}(:,1),[-2,1]);

for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.2,"betaLM",0.5);%設計変数微分の計算方法の指定
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);
    ungradetest = ungradetest.solveAnalysis(1,0,0);
    ungradetest.plotGeometry(1,obj.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end
for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.1,"betaLM",0.5);%設計変数微分の計算方法の指定
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);
    ungradetest = ungradetest.solveAnalysis(1,0,0);
    ungradetest.plotGeometry(1,obj.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end
for i = 1:10
    [nextVar,ungradetest] = ungradetest.calcNextVariables(@objFun,cmin,cmax,"TrustRegion",0.05,"betaLM",0.0);%設計変数微分の計算方法の指定
    ungradetest= ungradetest.updateMeshGeomfromVariables(nextVar);
    ungradetest = ungradetest.solveAnalysis(1,0,0);
    ungradetest.plotGeometry(1,obj.Cp{1}(:,1),[-2,1]);
    ungradetest.plotOptimizationState(2);
end
%}
