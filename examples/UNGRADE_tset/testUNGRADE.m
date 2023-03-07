%
%%%%%%%%ここから
clear;
orgVal = [4.0000,4.0000, 4.0000,4.0000,4.0000];
modVal = orgVal;
modVal(1) = orgVal(1)-0.8;
%orgVal = modifyDesFile("org.des","org.des");

ub = [  4,   4,   4,   4,   4];
lb = [0.5, 0.5, 0.5, 0.5, 0.5];
cmin = [ 0.4]';
cmax = [0.45]';
[~,orgVal] = modifyDesFile("org.des","org.des",orgVal);
CommandVSP = strcat("vsp wing.vsp3 -script makeOrgSurfandMesh.vspscript -des org.des"); %手順2のコマンド環境に合わせて各自書き換え
[~,~] = system(CommandVSP);
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "orgMesh.vspgeom", 0);
ungradetest = UNGRADE(@(x)vspSurfGen(x,"wing","org.des"),orgVal,lb,ub,p',con',id',wedata,1);
ungradetest = ungradetest.flowCondition(1,0.0001);
ungradetest = ungradetest.setREFS(80,20,4);
ungradetest = ungradetest.setRotationCenter([0,0,0]);
ungradetest = ungradetest.setCf(1,500000,0.2,0.052*(10^-5),0);

ungradetest = ungradetest.makeCluster(50,50);
ungradetest = ungradetest.makeEquation(20,5,10);

%ungradetest.checkSurfGenWork(1,1);
ungradetest = ungradetest.makeMeshGradient();
[modSurf1,modMesh1] = ungradetest.variables2Mesh(modVal,"raw");
[modSurf2,modMesh2] = ungradetest.variables2Mesh(modVal,"linear");
figure(1);clf;hold on;
ungradetest.viewMesh(modMesh1,1);
ungradetest.viewMesh(modMesh2,1);
figure(2);clf;hold on;
ungradetest.viewSurf(modSurf1,2);
ungradetest.viewSurf(modSurf2,2);


%}