function anss = UNLSI_calcFemResidual(frame,fem,dynPressure,delta)
    anss.dynPressure = dynPressure;
    nbVerts = numel(frame.usedVerts);
    if nargin > 3
        anss.delta_p = delta;
        anss.residual= frame.femLHS*delta-fem.femRHSp.*dynPressure; %à≥óÕÇ…ÇÊÇÈïœà ÇÃÇ›çló∂
        disp_buff(frame.InvMatIndex,1)=delta;
        anss.delta_Press = zeros(size(frame.verts));
        anss.delta_Press(frame.usedVerts,1)=disp_buff(1:nbVerts,1);
        anss.delta_Press(frame.usedVerts,2)=disp_buff(nbVerts+1:2*nbVerts,1);
        anss.delta_Press(frame.usedVerts,3)=disp_buff(2*nbVerts+1:3*nbVerts,1);
        anss.displacement = anss.delta_Press+frame.delta_SelfWeight;
    else
        anss.delta_p = frame.invFemLHS*(fem.femRHSp.*dynPressure);
        anss.residual= frame.femLHS*anss.delta_p-fem.femRHSp.*dynPressure;
        disp_buff(frame.InvMatIndex,1)=anss.delta_p;
        anss.delta_Press = zeros(size(frame.verts));
        anss.delta_Press(frame.usedVerts,1)=disp_buff(1:nbVerts,1);
        anss.delta_Press(frame.usedVerts,2)=disp_buff(nbVerts+1:2*nbVerts,1);
        anss.delta_Press(frame.usedVerts,3)=disp_buff(2*nbVerts+1:3*nbVerts,1);
        anss.displacement = anss.delta_Press+frame.delta_SelfWeight;
        
    end
    stress = reshape(frame.stressMat*anss.displacement(:),[3,size(frame.tri,1)])';
    anss.vonMises = sqrt(0.5*((stress(:,1)-stress(:,2)).^2+stress(:,1).^2+stress(:,2).^2+6.*stress(:,3).^2));
    anss.strain = reshape(frame.strainMat*anss.displacement(:),[3,size(frame.tri,1)])';
end