%
%%%%%%%%ここから
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "wing.vspgeom", 0);
wing = UNLSI(p',con',id',wedata,1);
wing.checkMesh(sqrt(eps));
wing = wing.makeCluster(50,50);
wing = wing.makeEquation(20,5,3);
%%%%%%%%ここまでは一度計算すればスキップできる
%}

wing = wing.flowCondition(1,0.0001);
wing = wing.flowCondition(2,4);
wing = wing.setREFS(80,20,4);
wing = wing.setRotationCenter([0,0,0]);
wing = wing.setCf(1,500000,4,0.052*(10^-5),0);
wing = wing.setCf(2,500000,4,0.052*(10^-5),0);
wing = wing.solveFlow(1,[0,5,10],[0,0,0]);
wing = wing.solveFlow(2,[0,5,10],[0,0,0]);
disp(wing.AERODATA{1});
disp(wing.AERODATA{2});
wing.plotGeometry(1,wing.Cp{1}(:,2),[-2,1]);
wing.plotGeometry(2,wing.Cp{2}(:,2),[-0.5,0.5]);