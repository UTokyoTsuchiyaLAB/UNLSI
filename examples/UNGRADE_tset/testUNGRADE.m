%
%%%%%%%%ここから
clear;
orgVal = [4.0000,4.0000, 4.0000,4.0000,4.0000];
%orgVal = modifyDesFile("org.des","org.des");

ub = [  4,   4,   4,   4,   4];
lb = [0.5, 0.5, 0.5, 0.5, 0.5];
cmin = [0.4]';
cmax = [0.45]';
[orgVerts,orgCon,surfID,wakelineID,orgVal] = vspMeshGen(orgVal);
ungradetest = UNGRADE(@(x)vspSurfGen(x,"wing","org.des"),orgVal,lb,ub,orgVerts,orgCon,surfID,wakelineID,1);
ungradetest = ungradetest.setMeshGenFun(@(x)vspMeshGen(x,"wing","org.des"));
ungradetest = ungradetest.setOptCondition(0.0001,5,0,20,5,10,50,50,500000,4,0.052*(10^-5),0,1);
ungradetest = ungradetest.setHessianUpdate(eye(numel(orgVal)),0.1,"SSR1");
ungradetest = ungradetest.setREFS(80,20,4);
ungradetest = ungradetest.setRotationCenter([0,0,0]);
ungradetest = ungradetest.setCf(1,ungradetest.unlsiParam.Re,ungradetest.unlsiParam.Lch,ungradetest.unlsiParam.k,ungradetest.unlsiParam.LTratio,ungradetest.unlsiParam.coefficient);
ungradetest = ungradetest.makeCluster(ungradetest.unlsiParam.nCluster,ungradetest.unlsiParam.edgeAngleThreshold);

for i = 1:50
    ungradetest = ungradetest.updateVariables(@objFun,cmin,cmax);
end
%ungradetest = ungradetest.solve

%}