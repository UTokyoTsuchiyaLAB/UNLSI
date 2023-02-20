%

%clear;
%[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "HST_nominal_ellipse.facet", 0);
%
%test = UNLSI((p./1000)',con',id',wedata,1);
%test.checkMesh(sqrt(eps));
%test = test.flowCondition(2,4.0);
%test = test.makeCluster(50,50);
%test = test.makeEquation(20,5,10);
%
test = test.setREFS(22.812,1.78,14.4);
test = test.setMomentCenter([14.4*0.6,0,0]);
test = test.setCf(1,500000,0.2,0.052*(10^-5),0);
alist = linspace(-5,15,10);
for i = 1:numel(alist)
    test = test.solveFlow(1,alist(i),0);
end

%test.plotGeometry(1,test.Cp,[-0.5,0.5]);
%}
