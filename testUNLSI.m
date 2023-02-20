%

clear;
[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "test.vspgeom", 0 );
%
test = UNLSI(p',con',id',wedata,1);
test.checkMesh(sqrt(eps));
test = test.flowCondition(1,0.001);
test = test.makeCluster(50,50);
test = test.makeEquation(20,5,10);
%}
test = test.setREFS(22.812,1.78,14.4);
test = test.setMomentCenter([14.4*0.6,0,0]);
test = test.setCf(1,500000,0.2,0.052*(10^-5),0);
test = test.solveFlow(1,5,0);
test.plotGeometry(1,test.Cp,[-0.3,0.3]);
%test.plotGeometry(1);