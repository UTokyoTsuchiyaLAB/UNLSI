clear;

verts = [0,0,0;
         1,0,0;
         0,1,0;
         0,0,1];
tri = [1,2,3;
       2,3,4;
       1,3,4];
surfID = [1;1;1;];
wakeline = [];
test = UNLSI(verts,tri,surfID,wakeline);
test = test.addFlowCondition(2);
test = test.addFlowCondition(5);
clf;
plot(test.flow{1}.pp.GridVectors{1},test.flow{1}.pp.Values);
hold on;
plot(test.flow{2}.pp.GridVectors{1},test.flow{2}.pp.Values);
hold off;