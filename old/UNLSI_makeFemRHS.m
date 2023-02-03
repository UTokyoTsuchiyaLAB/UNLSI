function fem = UNLSI_makeFemRHS(ansf,frame)
fem.g0 = 9.8;
fem.Fp = zeros(6*size(frame.usedVerts,2),1);
for i = 1:size(frame.tri,1)
    if frame.isWake(i,1) == 0
        fem.Fg{i}(1,1:3) = -1./3.*ansf.Cp(i)*frame.area(i,1).*frame.n(i,:);fem.Fg{i}(1,4:6) = 0;
        fem.Fg{i}(2,1:3) = -1./3.*ansf.Cp(i)*frame.area(i,1).*frame.n(i,:);fem.Fg{i}(2,4:6) = 0;
        fem.Fg{i}(3,1:3) = -1./3.*ansf.Cp(i)*frame.area(i,1).*frame.n(i,:);fem.Fg{i}(3,4:6) = 0;

        fem.Fp(frame.IndexRow{i}(:,1),1) = fem.Fp(frame.IndexRow{i}(:,1),1)+fem.Fg{i}(:);
    end
end
fem.femRHSp = fem.Fp(frame.MatIndex==1,1);
end