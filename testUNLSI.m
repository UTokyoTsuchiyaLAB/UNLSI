clear;


[con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "test.vspgeom", 0 );
test = UNLSI(p',con',id',wedata);
test = test.addFlowCondition(2);
test = test.addFlowCondition(5);
test = test.makeCluster(10);
clf;
plot(test.flow{1}.pp.GridVectors{1},test.flow{1}.pp.Values);
hold on;
plot(test.flow{2}.pp.GridVectors{1},test.flow{2}.pp.Values);
hold off;