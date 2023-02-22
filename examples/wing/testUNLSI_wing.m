%{
%%%%%%%%ここから
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "wing.vspgeom", 0);
wing = UNLSI(p',con',id',wedata,1);
wing.checkMesh(sqrt(eps));
wing = wing.makeCluster(50,50);
wing = wing.makeEquation(20,5,3);
%%%%%%%%ここまでは一度インスタンスを作成したらスキップできる
%}

wing = wing.flowCondition(1,0.0001);
wing = wing.flowCondition(2,4.0);
wing = wing.setREFS(72,18,4);
wing = wing.setRotationCenter([0,0,0]);
wing = wing.setCf(1,500000,0.2,0.052*(10^-5),0);
wing = wing.setCf(2,500000,0.2,0.052*(10^-5),0);
wing = wing.solveFlow(1,10,3);
wing.plotGeometry(1,wing.Cp,[-2,1]);