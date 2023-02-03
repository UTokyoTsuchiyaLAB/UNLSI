function frame = UNLSI_initFem(frame)
frame.g0 = 9.8;
frame.usedVerts = unique(frame.tri(frame.isWake == 0,:));
frame.usedVerts=frame.usedVerts(:)';
nbVerts = numel(frame.usedVerts);

for iter = 1:size(frame.tri,1)
    if frame.isWake(iter,1) == 0
        frame.IndexRow{iter} = zeros(18,18);
        frame.IndexCol{iter} = zeros(18,18);
        for i = 1:6
            for j = 1:3
                [r,c] = find(frame.usedVerts==frame.tri(iter,j));
                frame.IndexRow{iter}(3*(i-1)+j,:) = c+nbVerts*(i-1);
                frame.IndexCol{iter}(:,3*(i-1)+j) = c+nbVerts*(i-1);
            end
        end
    end
end
%‹«ŠEğŒİ’è
InvMatIndex = [];
for i = 1:size(frame.usedVerts,2)
   if  frame.verts(frame.usedVerts(i),2)<=sqrt(eps)*100
       MatIndex(i,1) = 0;
       MatIndex(1*nbVerts+i,1) = 0;
       MatIndex(2*nbVerts+i,1) = 0;
       MatIndex(3*nbVerts+i,1) = 0;
       MatIndex(4*nbVerts+i,1) = 0;
       MatIndex(5*nbVerts+i,1) = 0;
   else
       MatIndex(i,1) = 1;
       MatIndex(1*nbVerts+i,1) = 1;
       MatIndex(2*nbVerts+i,1) = 1;
       MatIndex(3*nbVerts+i,1) = 1;
       MatIndex(4*nbVerts+i,1) = 1;
       MatIndex(5*nbVerts+i,1) = 0;
       InvMatIndex=[InvMatIndex,i];
   end
end
frame.MatIndex = MatIndex;
frame.InvMatIndex=[InvMatIndex,1*nbVerts+InvMatIndex,2*nbVerts+InvMatIndex,3*nbVerts+InvMatIndex,4*nbVerts+InvMatIndex]';
end