
%
%%%%%%%%ここから
clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "Anemoi.vspgeom", 0);
anemoi = UNLSI(p',con',id',wedata,1);
anemoi.checkMesh(sqrt(eps));
anemoi = anemoi.makeCluster(50,50);
anemoi = anemoi.makeEquation(20,5,10);
%%%%%%%%ここまでは一度インスタンスを作成したらスキップできる
%}


anemoi = anemoi.flowCondition(1,0.0001);
anemoi = anemoi.setREFS(0.8,4,0.2);
anemoi = anemoi.setRotationCenter([-0.07,0,0]);
anemoi = anemoi.setCf(1,50000,0.2,0.052*(10^-5),0,1.5);
anemoi = anemoi.solveFlow(1,0,0);
anemoi.plotGeometry(1,anemoi.Cp,[-2,1]);
