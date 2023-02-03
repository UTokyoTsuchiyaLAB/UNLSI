function dInf_dverts = UNLSI_calcInfCoeffDerivative(frame)
orgVerts = frame.verts;
dv = sqrt(eps);
dInf_dverts = [];
for xyz = 1:3
    for iter = 1:size(frame.verts,1)
        dInf_dverts{iter,xyz}.infA = sparse(frame.nbPanel,frame.nbPanel);
        dInf_dverts{iter,xyz}.infB = sparse(frame.nbPanel,frame.nbPanel);
       [infA1,infB1,cnt] = calcVertsCoeffs_row(frame,iter,orgVerts,dv,xyz);
       [infA2,infB2] = calcVertsCoeffs_col(frame,iter,orgVerts,dv,xyz,cnt);
       dInf_dverts{iter,xyz}.infA(cnt,:) = infA1;dInf_dverts{iter,xyz}.infB(cnt,:) = infB1;
       dInf_dverts{iter,xyz}.infA(:,cnt) = infA2;dInf_dverts{iter,xyz}.infB(:,cnt) = infB2;
       dInf_dverts{iter,xyz}.infA = (dInf_dverts{iter,xyz}.infA-frame.infA)./dv;
       dInf_dverts{iter,xyz}.infB = (dInf_dverts{iter,xyz}.infB-frame.infB)./dv;
    end
end
end

function [infA,infB,cnt] = calcVertsCoeffs_row(frame,iter,orgVerts,dv,xyz)
   cnt = [];
   for j = 1:frame.nbPanel
      if any(frame.tri(j,:)==iter)
          cnt = [cnt,j];
      end
   end
   if not(isempty(cnt))
        frame.verts = orgVerts;
        frame.verts(iter,xyz) = orgVerts(iter,xyz)+dv;
        frame.areasum = 0;
        frame.volume = 0;
        frame.pntX = [frame.verts(frame.tri(:,1),1),frame.verts(frame.tri(:,2),1),frame.verts(frame.tri(:,3),1)];
        frame.pntY = [frame.verts(frame.tri(:,1),2),frame.verts(frame.tri(:,2),2),frame.verts(frame.tri(:,3),2)];
        frame.pntZ = [frame.verts(frame.tri(:,1),3),frame.verts(frame.tri(:,2),3),frame.verts(frame.tri(:,3),3)];
        frame.center = [mean(frame.pntX,2),mean(frame.pntY,2),mean(frame.pntZ,2)];
        for i = cnt
            [b1, b2, frame.n(i,:)] = vartex_iter(frame.verts(frame.tri(i,1),:),frame.verts(frame.tri(i,2),:),frame.verts(frame.tri(i,3),:));
        end
        %%
        %Ç‹Ç∏ÇÕâ°
        for i = 1:numel(cnt)
            POI.X(i,1) = frame.center(i,1);
            POI.Y(i,1) = frame.center(i,2);
            POI.Z(i,1) = frame.center(i,3);
        end
        c.X(1,1:frame.nbPanel) = frame.center(1:frame.nbPanel,1)';
        c.Y(1,1:frame.nbPanel) = frame.center(1:frame.nbPanel,2)';
        c.Z(1,1:frame.nbPanel) = frame.center(1:frame.nbPanel,3)';
        n.X(1,1:frame.nbPanel) = frame.n(1:frame.nbPanel,1)';
        n.Y(1,1:frame.nbPanel) = frame.n(1:frame.nbPanel,2)';
        n.Z(1,1:frame.nbPanel) = frame.n(1:frame.nbPanel,3)';
        N1.X(1,1:frame.nbPanel) = frame.pntX(1:frame.nbPanel,1)';
        N1.Y(1,1:frame.nbPanel) = frame.pntY(1:frame.nbPanel,1)';
        N1.Z(1,1:frame.nbPanel) = frame.pntZ(1:frame.nbPanel,1)';
        N2.X(1,1:frame.nbPanel) = frame.pntX(1:frame.nbPanel,2)';
        N2.Y(1,1:frame.nbPanel) = frame.pntY(1:frame.nbPanel,2)';
        N2.Z(1,1:frame.nbPanel) = frame.pntZ(1:frame.nbPanel,2)';
        N3.X(1,1:frame.nbPanel) = frame.pntX(1:frame.nbPanel,3)';
        N3.Y(1,1:frame.nbPanel) = frame.pntY(1:frame.nbPanel,3)';
        N3.Z(1,1:frame.nbPanel) = frame.pntZ(1:frame.nbPanel,3)';
        POI.X = repmat(POI.X,[1,frame.nbPanel]);
        POI.Y = repmat(POI.Y,[1,frame.nbPanel]);
        POI.Z = repmat(POI.Z,[1,frame.nbPanel]);
        c.X = repmat(c.X,[numel(cnt),1]);
        c.Y = repmat(c.Y,[numel(cnt),1]);
        c.Z = repmat(c.Z,[numel(cnt),1]);
        n.X = repmat(n.X,[numel(cnt),1]);
        n.Y = repmat(n.Y,[numel(cnt),1]);
        n.Z = repmat(n.Z,[numel(cnt),1]);
        N1.X = repmat(N1.X,[numel(cnt),1]);
        N1.Y = repmat(N1.Y,[numel(cnt),1]);
        N1.Z = repmat(N1.Z,[numel(cnt),1]);
        N2.X = repmat(N2.X,[numel(cnt),1]);
        N2.Y = repmat(N2.Y,[numel(cnt),1]);
        N2.Z = repmat(N2.Z,[numel(cnt),1]);
        N3.X = repmat(N3.X,[numel(cnt),1]);
        N3.Y = repmat(N3.Y,[numel(cnt),1]);
        N3.Z = repmat(N3.Z,[numel(cnt),1]);

        n12.X = (N1.X+N2.X)./2;
        n12.Y = (N1.Y+N2.Y)./2;
        n12.Z = (N1.Z+N2.Z)./2;
        m = getUnitVector(c,n12);
        l = matrix_cross(m,n);
        pjk.X = POI.X-c.X;
        pjk.Y = POI.Y-c.Y;
        pjk.Z = POI.Z-c.Z;
        PN = matrix_dot(pjk,n);

        %1âÒñ⁄
        a.X = POI.X-N1.X;
        a.Y = POI.Y-N1.Y;
        a.Z = POI.Z-N1.Z;
        b.X = POI.X-N2.X;
        b.Y = POI.Y-N2.Y;
        b.Z = POI.Z-N2.Z;
        s.X = N2.X-N1.X;
        s.Y = N2.Y-N1.Y;
        s.Z = N2.Z-N1.Z;
        Al = matrix_dot(n,matrix_cross(s,a));

        PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
        PB = PA-Al.*matrix_dot(s,m);
        num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
        denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
        phiV = atan2(num,denom);
        A = matrix_norm(a);
        B = matrix_norm(b);
        S = matrix_norm(s);
        GL = 1./S.*log(((A+B+S)./(A+B-S+sqrt(eps))));
        srcV = Al.*GL-PN.*phiV;
        VortexA = phiV;
        VortexB = srcV;
        %2âÒñ⁄
        a.X = POI.X-N2.X;
        a.Y = POI.Y-N2.Y;
        a.Z = POI.Z-N2.Z;
        b.X = POI.X-N3.X;
        b.Y = POI.Y-N3.Y;
        b.Z = POI.Z-N3.Z;
        s.X = N3.X-N2.X;
        s.Y = N3.Y-N2.Y;
        s.Z = N3.Z-N2.Z;

        Al = matrix_dot(n,matrix_cross(s,a));

        PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
        PB = PA-Al.*matrix_dot(s,m);
        num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
        denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
        phiV = atan2(num,denom);
        A = matrix_norm(a);
        B = matrix_norm(b);
        S = matrix_norm(s);
        GL = 1./S.*log(((A+B+S)./(A+B-S+sqrt(eps))));
        srcV = Al.*GL-PN.*phiV;
        VortexA = VortexA+phiV;
        VortexB = VortexB+srcV;
        %3âÒñ⁄
        a.X = POI.X-N3.X;
        a.Y = POI.Y-N3.Y;
        a.Z = POI.Z-N3.Z;
        b.X = POI.X-N1.X;
        b.Y = POI.Y-N1.Y;
        b.Z = POI.Z-N1.Z;
        s.X = N1.X-N3.X;
        s.Y = N1.Y-N3.Y;
        s.Z = N1.Z-N3.Z;

        Al = matrix_dot(n,matrix_cross(s,a));

        PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
        PB = PA-Al.*matrix_dot(s,m);
        num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
        denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
        phiV = atan2(num,denom);
        A = matrix_norm(a);
        B = matrix_norm(b);
        S = matrix_norm(s);
        GL = 1./S.*log(((A+B+S)./(A+B-S+sqrt(eps))));
        srcV = Al.*GL-PN.*phiV;
        VortexA = VortexA+phiV;
        VortexB = VortexB+srcV;

        if frame.isSymm == 1
            c.Y = -c.Y;
            n.Y = -n.Y;
            N1.Y = -N1.Y;
            N2.Y = -N2.Y;
            N3.Y = -N3.Y;
            n12.X = (N1.X+N2.X)./2;
            n12.Y = (N1.Y+N2.Y)./2;
            n12.Z = (N1.Z+N2.Z)./2;
            m = getUnitVector(c,n12);
            l = matrix_cross(m,n);
            pjk.X = POI.X-c.X;
            pjk.Y = POI.Y-c.Y;
            pjk.Z = POI.Z-c.Z;
            PN = matrix_dot(pjk,n);

            %1âÒñ⁄
            a.X = POI.X-N1.X;
            a.Y = POI.Y-N1.Y;
            a.Z = POI.Z-N1.Z;
            b.X = POI.X-N2.X;
            b.Y = POI.Y-N2.Y;
            b.Z = POI.Z-N2.Z;
            s.X = N2.X-N1.X;
            s.Y = N2.Y-N1.Y;
            s.Z = N2.Z-N1.Z;
            Al = matrix_dot(n,matrix_cross(s,a));

            PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
            PB = PA-Al.*matrix_dot(s,m);
            num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
            denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
            phiV = atan2(num,denom);
            A = matrix_norm(a);
            B = matrix_norm(b);
            S = matrix_norm(s);
            GL = 1./S.*log(((A+B+S)./(A+B-S+sqrt(eps))));
            srcV = Al.*GL-PN.*phiV;
            VortexA = VortexA-phiV;
            VortexB = VortexB-srcV;
            %2âÒñ⁄
            a.X = POI.X-N2.X;
            a.Y = POI.Y-N2.Y;
            a.Z = POI.Z-N2.Z;
            b.X = POI.X-N3.X;
            b.Y = POI.Y-N3.Y;
            b.Z = POI.Z-N3.Z;
            s.X = N3.X-N2.X;
            s.Y = N3.Y-N2.Y;
            s.Z = N3.Z-N2.Z;

            Al = matrix_dot(n,matrix_cross(s,a));

            PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
            PB = PA-Al.*matrix_dot(s,m);
            num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
            denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
            phiV = atan2(num,denom);
            A = matrix_norm(a);
            B = matrix_norm(b);
            S = matrix_norm(s);
            GL = 1./S.*log(((A+B+S)./(A+B-S+sqrt(eps))));
            srcV = Al.*GL-PN.*phiV;
            VortexA = VortexA-phiV;
            VortexB = VortexB-srcV;
            %3âÒñ⁄
            a.X = POI.X-N3.X;
            a.Y = POI.Y-N3.Y;
            a.Z = POI.Z-N3.Z;
            b.X = POI.X-N1.X;
            b.Y = POI.Y-N1.Y;
            b.Z = POI.Z-N1.Z;
            s.X = N1.X-N3.X;
            s.Y = N1.Y-N3.Y;
            s.Z = N1.Z-N3.Z;

            Al = matrix_dot(n,matrix_cross(s,a));

            PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
            PB = PA-Al.*matrix_dot(s,m);
            num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
            denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
            phiV = atan2(num,denom);
            A = matrix_norm(a);
            B = matrix_norm(b);
            S = matrix_norm(s);
            GL = 1./S.*log(((A+B+S)./(A+B-S+sqrt(eps))));
            srcV = Al.*GL-PN.*phiV;
            VortexA = VortexA-phiV;
            VortexB = VortexB-srcV;
        end
        infA = VortexA./4./pi; 
        infB = VortexB./4./pi;
        eyeRbuff = eye(frame.nbPanel);
        eyeR=eyeRbuff(cnt,:);
        infA(find(eyeR)) = -0.5;
        
        for j = 1:size(frame.tri,1)
            if frame.isWake(j,1) == 1
                interpCoeff = zeros(1,size(frame.wakelines,2));
                interpID = zeros(2,size(frame.wakelines,2));
                for i = 1:size(frame.wakelines,2)
                    upLag = 1;
                    loLag = 1;
                    if frame.surfID(j)==frame.wakelines{i}.ID
                        iwakey = (frame.verts(frame.wakelines{i}.nodeID(1),2)+frame.verts(frame.wakelines{i}.nodeID(2),2))/2;
                        for k = 1:size(frame.wakelines,2)
                            if frame.surfID(j)==frame.wakelines{k}.ID && not(i==k)
                                    kwakey = (frame.verts(frame.wakelines{k}.nodeID(1),2)+frame.verts(frame.wakelines{k}.nodeID(2),2))/2;
                                    upLag = upLag*(frame.center(j,2)-kwakey);
                                    loLag = loLag*(iwakey-kwakey);
                            end
                        end
                        interpCoeff(i) = upLag/loLag;
                        interpID(:,i) = [frame.wakelines{i}.upperID;frame.wakelines{i}.lowerID];
                    end
                end
                Nw1.X = frame.verts(frame.tri(j,1),1);
                Nw1.Y = frame.verts(frame.tri(j,1),2);
                Nw1.Z = frame.verts(frame.tri(j,1),3);
                Nw2.X = frame.verts(frame.tri(j,2),1);
                Nw2.Y = frame.verts(frame.tri(j,2),2);
                Nw2.Z = frame.verts(frame.tri(j,2),3);
                Nw3.X = frame.verts(frame.tri(j,3),1);
                Nw3.Y = frame.verts(frame.tri(j,3),2);
                Nw3.Z = frame.verts(frame.tri(j,3),3);
                n.X = frame.n(j,1);
                n.Y = frame.n(j,2);
                n.Z = frame.n(j,3);
                c.X = frame.center(j,1);
                c.Y = frame.center(j,2);
                c.Z = frame.center(j,3);
                n12.X = (Nw1.X+Nw2.X)./2;
                n12.Y = (Nw1.Y+Nw2.Y)./2;
                n12.Z = (Nw1.Z+Nw2.Z)./2;
                m = getUnitVector(c,n12);
                l = matrix_cross(m,n);
                for i = cnt
                    if frame.isBody(i,1)==1
                        pjk.X = frame.center(i,1)-c.X;
                        pjk.Y = frame.center(i,2)-c.Y;
                        pjk.Z = frame.center(i,3)-c.Z;
                        PN = matrix_dot(pjk,n);
                       %1âÒñ⁄
                        a.X = frame.center(i,1)-Nw1.X;
                        a.Y = frame.center(i,2)-Nw1.Y;
                        a.Z = frame.center(i,3)-Nw1.Z;
                        b.X = frame.center(i,1)-Nw2.X;
                        b.Y = frame.center(i,2)-Nw2.Y;
                        b.Z = frame.center(i,3)-Nw2.Z;
                        s.X = Nw2.X-Nw1.X;
                        s.Y = Nw2.Y-Nw1.Y;
                        s.Z = Nw2.Z-Nw1.Z;
                        Al = matrix_dot(n,matrix_cross(s,a));

                        PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
                        PB = PA-Al.*matrix_dot(s,m);
                        num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
                        denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
                        phiV = atan2(num,denom);
                        VortexA = phiV;
                        %2âÒñ⁄
                        a.X = frame.center(i,1)-Nw2.X;
                        a.Y = frame.center(i,2)-Nw2.Y;
                        a.Z = frame.center(i,3)-Nw2.Z;
                        b.X = frame.center(i,1)-Nw3.X;
                        b.Y = frame.center(i,2)-Nw3.Y;
                        b.Z = frame.center(i,3)-Nw3.Z;
                        s.X = Nw3.X-Nw2.X;
                        s.Y = Nw3.Y-Nw2.Y;
                        s.Z = Nw3.Z-Nw2.Z;
                        Al = matrix_dot(n,matrix_cross(s,a));

                        PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
                        PB = PA-Al.*matrix_dot(s,m);
                        num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
                        denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
                        phiV = atan2(num,denom);
                        VortexA = VortexA+phiV;
                        %3âÒñ⁄
                        a.X = frame.center(i,1)-Nw3.X;
                        a.Y = frame.center(i,2)-Nw3.Y;
                        a.Z = frame.center(i,3)-Nw3.Z;
                        b.X = frame.center(i,1)-Nw1.X;
                        b.Y = frame.center(i,2)-Nw1.Y;
                        b.Z = frame.center(i,3)-Nw1.Z;
                        s.X = Nw1.X-Nw3.X;
                        s.Y = Nw1.Y-Nw3.Y;
                        s.Z = Nw1.Z-Nw3.Z;
                        Al = matrix_dot(n,matrix_cross(s,a));

                        PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
                        PB = PA-Al.*matrix_dot(s,m);
                        num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
                        denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
                        phiV = atan2(num,denom);
                        VortexA = VortexA+phiV;
                        influence = VortexA./pi./4;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        A = zeros(1,frame.nbPanel);B = zeros(1,frame.nbPanel);
                        A(1,interpID(1,find(interpID(1,:)))) =  + influence*interpCoeff(find(interpCoeff));
                        B(1,interpID(2,find(interpID(2,:)))) =  - influence*interpCoeff(find(interpCoeff));
                        infA(cnt==i,:) = infA(cnt==i,:)+A(1,:);
                        infA(cnt==i,:) = infA(cnt==i,:)+B(1,:);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
               end
               if frame.isSymm == 1
                    Nw1.X = frame.verts(frame.tri(j,1),1);
                    Nw1.Y = -frame.verts(frame.tri(j,1),2);
                    Nw1.Z = frame.verts(frame.tri(j,1),3);
                    Nw2.X = frame.verts(frame.tri(j,2),1);
                    Nw2.Y = -frame.verts(frame.tri(j,2),2);
                    Nw2.Z = frame.verts(frame.tri(j,2),3);
                    Nw3.X = frame.verts(frame.tri(j,3),1);
                    Nw3.Y = -frame.verts(frame.tri(j,3),2);
                    Nw3.Z = frame.verts(frame.tri(j,3),3);
                    n.X = frame.n(j,1);
                    n.Y = -frame.n(j,2);
                    n.Z = frame.n(j,3);
                    c.X = frame.center(j,1);
                    c.Y = -frame.center(j,2);
                    c.Z = frame.center(j,3);
                    n12.X = (Nw1.X+Nw2.X)./2;
                    n12.Y = (Nw1.Y+Nw2.Y)./2;
                    n12.Z = (Nw1.Z+Nw2.Z)./2;
                    m = getUnitVector(c,n12);
                    l = matrix_cross(m,n);
                    for i = cnt
                        if frame.isBody(i,1)==1
                            pjk.X = frame.center(i,1)-c.X;
                            pjk.Y = frame.center(i,2)-c.Y;
                            pjk.Z = frame.center(i,3)-c.Z;
                            PN = matrix_dot(pjk,n);
                            %1âÒñ⁄
                            a.X = frame.center(i,1)-Nw1.X;
                            a.Y = frame.center(i,2)-Nw1.Y;
                            a.Z = frame.center(i,3)-Nw1.Z;
                            b.X = frame.center(i,1)-Nw2.X;
                            b.Y = frame.center(i,2)-Nw2.Y;
                            b.Z = frame.center(i,3)-Nw2.Z;
                            s.X = Nw2.X-Nw1.X;
                            s.Y = Nw2.Y-Nw1.Y;
                            s.Z = Nw2.Z-Nw1.Z;
                            Al = matrix_dot(n,matrix_cross(s,a));

                            PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
                            PB = PA-Al.*matrix_dot(s,m);
                            num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
                            denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
                            phiV = atan2(num,denom);
                            VortexA = phiV;
                            %2âÒñ⁄
                            a.X = frame.center(i,1)-Nw2.X;
                            a.Y = frame.center(i,2)-Nw2.Y;
                            a.Z = frame.center(i,3)-Nw2.Z;
                            b.X = frame.center(i,1)-Nw3.X;
                            b.Y = frame.center(i,2)-Nw3.Y;
                            b.Z = frame.center(i,3)-Nw3.Z;
                            s.X = Nw3.X-Nw2.X;
                            s.Y = Nw3.Y-Nw2.Y;
                            s.Z = Nw3.Z-Nw2.Z;
                            Al = matrix_dot(n,matrix_cross(s,a));

                            PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
                            PB = PA-Al.*matrix_dot(s,m);
                            num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
                            denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
                            phiV = atan2(num,denom);
                            VortexA = VortexA+phiV;
                            %3âÒñ⁄
                            a.X = frame.center(i,1)-Nw3.X;
                            a.Y = frame.center(i,2)-Nw3.Y;
                            a.Z = frame.center(i,3)-Nw3.Z;
                            b.X = frame.center(i,1)-Nw1.X;
                            b.Y = frame.center(i,2)-Nw1.Y;
                            b.Z = frame.center(i,3)-Nw1.Z;
                            s.X = Nw1.X-Nw3.X;
                            s.Y = Nw1.Y-Nw3.Y;
                            s.Z = Nw1.Z-Nw3.Z;
                            Al = matrix_dot(n,matrix_cross(s,a));

                            PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
                            PB = PA-Al.*matrix_dot(s,m);
                            num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
                            denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
                            phiV = atan2(num,denom);
                            VortexA = VortexA+phiV;
                            influence = VortexA./pi./4;
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            A = zeros(1,frame.nbPanel);B = zeros(1,frame.nbPanel);
                            A(1,interpID(1,find(interpID(1,:)))) =  + influence*interpCoeff(find(interpCoeff));
                            B(1,interpID(2,find(interpID(2,:)))) =  - influence*interpCoeff(find(interpCoeff));
                            infA(cnt==i,:) = infA(cnt==i,:)+A(1,:);
                            infA(cnt==i,:) = infA(cnt==i,:)+B(1,:);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        end
                    end
               end
            end
        end
   else
       infA = [];
       infB = [];       
   end
end

function [infA,infB] = calcVertsCoeffs_col(frame,iter,orgVerts,dv,xyz,cnt)
   if not(isempty(cnt))
        frame.verts = orgVerts;
        frame.verts(iter,xyz) = orgVerts(iter,xyz)+dv;
        frame.areasum = 0;
        frame.volume = 0;
        frame.pntX = [frame.verts(frame.tri(:,1),1),frame.verts(frame.tri(:,2),1),frame.verts(frame.tri(:,3),1)];
        frame.pntY = [frame.verts(frame.tri(:,1),2),frame.verts(frame.tri(:,2),2),frame.verts(frame.tri(:,3),2)];
        frame.pntZ = [frame.verts(frame.tri(:,1),3),frame.verts(frame.tri(:,2),3),frame.verts(frame.tri(:,3),3)];
        frame.center = [mean(frame.pntX,2),mean(frame.pntY,2),mean(frame.pntZ,2)];
        for i = cnt
            [b1, b2, frame.n(i,:)] = vartex_iter(frame.verts(frame.tri(i,1),:),frame.verts(frame.tri(i,2),:),frame.verts(frame.tri(i,3),:));
        end
        %%
        POI.X(1:frame.nbPanel,1) = frame.center(1:frame.nbPanel,1);
        POI.Y(1:frame.nbPanel,1) = frame.center(1:frame.nbPanel,2);
        POI.Z(1:frame.nbPanel,1) = frame.center(1:frame.nbPanel,3);
        c.X(1,:) = frame.center(cnt,1)';
        c.Y(1,:) = frame.center(cnt,2)';
        c.Z(1,:) = frame.center(cnt,3)';
        n.X(1,:) = frame.n(cnt,1)';
        n.Y(1,:) = frame.n(cnt,2)';
        n.Z(1,:) = frame.n(cnt,3)';
        N1.X(1,:) = frame.pntX(cnt,1)';
        N1.Y(1,:) = frame.pntY(cnt,1)';
        N1.Z(1,:) = frame.pntZ(cnt,1)';
        N2.X(1,:) = frame.pntX(cnt,2)';
        N2.Y(1,:) = frame.pntY(cnt,2)';
        N2.Z(1,:) = frame.pntZ(cnt,2)';
        N3.X(1,:) = frame.pntX(cnt,3)';
        N3.Y(1,:) = frame.pntY(cnt,3)';
        N3.Z(1,:) = frame.pntZ(cnt,3)';
        POI.X = repmat(POI.X,[1,numel(cnt)]);
        POI.Y = repmat(POI.Y,[1,numel(cnt)]);
        POI.Z = repmat(POI.Z,[1,numel(cnt)]);
        c.X = repmat(c.X,[frame.nbPanel,1]);
        c.Y = repmat(c.Y,[frame.nbPanel,1]);
        c.Z = repmat(c.Z,[frame.nbPanel,1]);
        n.X = repmat(n.X,[frame.nbPanel,1]);
        n.Y = repmat(n.Y,[frame.nbPanel,1]);
        n.Z = repmat(n.Z,[frame.nbPanel,1]);
        N1.X = repmat(N1.X,[frame.nbPanel,1]);
        N1.Y = repmat(N1.Y,[frame.nbPanel,1]);
        N1.Z = repmat(N1.Z,[frame.nbPanel,1]);
        N2.X = repmat(N2.X,[frame.nbPanel,1]);
        N2.Y = repmat(N2.Y,[frame.nbPanel,1]);
        N2.Z = repmat(N2.Z,[frame.nbPanel,1]);
        N3.X = repmat(N3.X,[frame.nbPanel,1]);
        N3.Y = repmat(N3.Y,[frame.nbPanel,1]);
        N3.Z = repmat(N3.Z,[frame.nbPanel,1]);

        n12.X = (N1.X+N2.X)./2;
        n12.Y = (N1.Y+N2.Y)./2;
        n12.Z = (N1.Z+N2.Z)./2;
        m = getUnitVector(c,n12);
        l = matrix_cross(m,n);
        pjk.X = POI.X-c.X;
        pjk.Y = POI.Y-c.Y;
        pjk.Z = POI.Z-c.Z;
        PN = matrix_dot(pjk,n);

        %1âÒñ⁄
        a.X = POI.X-N1.X;
        a.Y = POI.Y-N1.Y;
        a.Z = POI.Z-N1.Z;
        b.X = POI.X-N2.X;
        b.Y = POI.Y-N2.Y;
        b.Z = POI.Z-N2.Z;
        s.X = N2.X-N1.X;
        s.Y = N2.Y-N1.Y;
        s.Z = N2.Z-N1.Z;
        Al = matrix_dot(n,matrix_cross(s,a));

        PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
        PB = PA-Al.*matrix_dot(s,m);
        num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
        denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
        phiV = atan2(num,denom);
        A = matrix_norm(a);
        B = matrix_norm(b);
        S = matrix_norm(s);
        GL = 1./S.*log(((A+B+S)./(A+B-S+sqrt(eps))));
        srcV = Al.*GL-PN.*phiV;
        VortexA = phiV;
        VortexB = srcV;
        %2âÒñ⁄
        a.X = POI.X-N2.X;
        a.Y = POI.Y-N2.Y;
        a.Z = POI.Z-N2.Z;
        b.X = POI.X-N3.X;
        b.Y = POI.Y-N3.Y;
        b.Z = POI.Z-N3.Z;
        s.X = N3.X-N2.X;
        s.Y = N3.Y-N2.Y;
        s.Z = N3.Z-N2.Z;

        Al = matrix_dot(n,matrix_cross(s,a));

        PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
        PB = PA-Al.*matrix_dot(s,m);
        num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
        denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
        phiV = atan2(num,denom);
        A = matrix_norm(a);
        B = matrix_norm(b);
        S = matrix_norm(s);
        GL = 1./S.*log(((A+B+S)./(A+B-S+sqrt(eps))));
        srcV = Al.*GL-PN.*phiV;
        VortexA = VortexA+phiV;
        VortexB = VortexB+srcV;
        %3âÒñ⁄
        a.X = POI.X-N3.X;
        a.Y = POI.Y-N3.Y;
        a.Z = POI.Z-N3.Z;
        b.X = POI.X-N1.X;
        b.Y = POI.Y-N1.Y;
        b.Z = POI.Z-N1.Z;
        s.X = N1.X-N3.X;
        s.Y = N1.Y-N3.Y;
        s.Z = N1.Z-N3.Z;

        Al = matrix_dot(n,matrix_cross(s,a));

        PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
        PB = PA-Al.*matrix_dot(s,m);
        num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
        denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
        phiV = atan2(num,denom);
        A = matrix_norm(a);
        B = matrix_norm(b);
        S = matrix_norm(s);
        GL = 1./S.*log(((A+B+S)./(A+B-S+sqrt(eps))));
        srcV = Al.*GL-PN.*phiV;
        VortexA = VortexA+phiV;
        VortexB = VortexB+srcV;

        if frame.isSymm == 1
            c.Y = -c.Y;
            n.Y = -n.Y;
            N1.Y = -N1.Y;
            N2.Y = -N2.Y;
            N3.Y = -N3.Y;
            n12.X = (N1.X+N2.X)./2;
            n12.Y = (N1.Y+N2.Y)./2;
            n12.Z = (N1.Z+N2.Z)./2;
            m = getUnitVector(c,n12);
            l = matrix_cross(m,n);
            pjk.X = POI.X-c.X;
            pjk.Y = POI.Y-c.Y;
            pjk.Z = POI.Z-c.Z;
            PN = matrix_dot(pjk,n);

            %1âÒñ⁄
            a.X = POI.X-N1.X;
            a.Y = POI.Y-N1.Y;
            a.Z = POI.Z-N1.Z;
            b.X = POI.X-N2.X;
            b.Y = POI.Y-N2.Y;
            b.Z = POI.Z-N2.Z;
            s.X = N2.X-N1.X;
            s.Y = N2.Y-N1.Y;
            s.Z = N2.Z-N1.Z;
            Al = matrix_dot(n,matrix_cross(s,a));

            PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
            PB = PA-Al.*matrix_dot(s,m);
            num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
            denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
            phiV = atan2(num,denom);
            A = matrix_norm(a);
            B = matrix_norm(b);
            S = matrix_norm(s);
            GL = 1./S.*log(((A+B+S)./(A+B-S+sqrt(eps))));
            srcV = Al.*GL-PN.*phiV;
            VortexA = VortexA-phiV;
            VortexB = VortexB-srcV;
            %2âÒñ⁄
            a.X = POI.X-N2.X;
            a.Y = POI.Y-N2.Y;
            a.Z = POI.Z-N2.Z;
            b.X = POI.X-N3.X;
            b.Y = POI.Y-N3.Y;
            b.Z = POI.Z-N3.Z;
            s.X = N3.X-N2.X;
            s.Y = N3.Y-N2.Y;
            s.Z = N3.Z-N2.Z;

            Al = matrix_dot(n,matrix_cross(s,a));

            PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
            PB = PA-Al.*matrix_dot(s,m);
            num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
            denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
            phiV = atan2(num,denom);
            A = matrix_norm(a);
            B = matrix_norm(b);
            S = matrix_norm(s);
            GL = 1./S.*log(((A+B+S)./(A+B-S+sqrt(eps))));
            srcV = Al.*GL-PN.*phiV;
            VortexA = VortexA-phiV;
            VortexB = VortexB-srcV;
            %3âÒñ⁄
            a.X = POI.X-N3.X;
            a.Y = POI.Y-N3.Y;
            a.Z = POI.Z-N3.Z;
            b.X = POI.X-N1.X;
            b.Y = POI.Y-N1.Y;
            b.Z = POI.Z-N1.Z;
            s.X = N1.X-N3.X;
            s.Y = N1.Y-N3.Y;
            s.Z = N1.Z-N3.Z;

            Al = matrix_dot(n,matrix_cross(s,a));

            PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
            PB = PA-Al.*matrix_dot(s,m);
            num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
            denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
            phiV = atan2(num,denom);
            A = matrix_norm(a);
            B = matrix_norm(b);
            S = matrix_norm(s);
            GL = 1./S.*log(((A+B+S)./(A+B-S+sqrt(eps))));
            srcV = Al.*GL-PN.*phiV;
            VortexA = VortexA-phiV;
            VortexB = VortexB-srcV;
        end
        infA = VortexA./4./pi; 
        infB = VortexB./4./pi;
        eyeRbuff = eye(frame.nbPanel);
        eyeR=eyeRbuff(:,cnt);
        infA(find(eyeR)) = -0.5;
        for j = 1:size(frame.tri,1)
            if frame.isWake(j,1) == 1
                interpCoeff = zeros(1,size(frame.wakelines,2));
                interpID = zeros(2,size(frame.wakelines,2));
                for i = 1:size(frame.wakelines,2)
                    upLag = 1;
                    loLag = 1;
                    if frame.surfID(j)==frame.wakelines{i}.ID
                        iwakey = (frame.verts(frame.wakelines{i}.nodeID(1),2)+frame.verts(frame.wakelines{i}.nodeID(2),2))/2;
                        for k = 1:size(frame.wakelines,2)
                            if frame.surfID(j)==frame.wakelines{k}.ID && not(i==k)
                                    kwakey = (frame.verts(frame.wakelines{k}.nodeID(1),2)+frame.verts(frame.wakelines{k}.nodeID(2),2))/2;
                                    upLag = upLag*(frame.center(j,2)-kwakey);
                                    loLag = loLag*(iwakey-kwakey);
                            end
                        end
                        interpCoeff(i) = upLag/loLag;
                        interpID(:,i) = [frame.wakelines{i}.upperID;frame.wakelines{i}.lowerID];
                    end
                end
                isIncluded = 0;
                for i = 1:numel(cnt)
                   if or(any(cnt(i)==interpID(1,:)),any(cnt(i)==interpID(2,:)))
                       isIncluded = 1;break;
                   end
                end
                if isIncluded==1
                    Nw1.X = frame.verts(frame.tri(j,1),1);
                    Nw1.Y = frame.verts(frame.tri(j,1),2);
                    Nw1.Z = frame.verts(frame.tri(j,1),3);
                    Nw2.X = frame.verts(frame.tri(j,2),1);
                    Nw2.Y = frame.verts(frame.tri(j,2),2);
                    Nw2.Z = frame.verts(frame.tri(j,2),3);
                    Nw3.X = frame.verts(frame.tri(j,3),1);
                    Nw3.Y = frame.verts(frame.tri(j,3),2);
                    Nw3.Z = frame.verts(frame.tri(j,3),3);
                    n.X = frame.n(j,1);
                    n.Y = frame.n(j,2);
                    n.Z = frame.n(j,3);
                    c.X = frame.center(j,1);
                    c.Y = frame.center(j,2);
                    c.Z = frame.center(j,3);
                    n12.X = (Nw1.X+Nw2.X)./2;
                    n12.Y = (Nw1.Y+Nw2.Y)./2;
                    n12.Z = (Nw1.Z+Nw2.Z)./2;
                    m = getUnitVector(c,n12);
                    l = matrix_cross(m,n);
                    for i = 1:frame.nbPanel
                        if frame.isBody(i,1)==1
                            pjk.X = frame.center(i,1)-c.X;
                            pjk.Y = frame.center(i,2)-c.Y;
                            pjk.Z = frame.center(i,3)-c.Z;
                            PN = matrix_dot(pjk,n);
                           %1âÒñ⁄
                            a.X = frame.center(i,1)-Nw1.X;
                            a.Y = frame.center(i,2)-Nw1.Y;
                            a.Z = frame.center(i,3)-Nw1.Z;
                            b.X = frame.center(i,1)-Nw2.X;
                            b.Y = frame.center(i,2)-Nw2.Y;
                            b.Z = frame.center(i,3)-Nw2.Z;
                            s.X = Nw2.X-Nw1.X;
                            s.Y = Nw2.Y-Nw1.Y;
                            s.Z = Nw2.Z-Nw1.Z;
                            Al = matrix_dot(n,matrix_cross(s,a));

                            PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
                            PB = PA-Al.*matrix_dot(s,m);
                            num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
                            denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
                            phiV = atan2(num,denom);
                            VortexA = phiV;
                            %2âÒñ⁄
                            a.X = frame.center(i,1)-Nw2.X;
                            a.Y = frame.center(i,2)-Nw2.Y;
                            a.Z = frame.center(i,3)-Nw2.Z;
                            b.X = frame.center(i,1)-Nw3.X;
                            b.Y = frame.center(i,2)-Nw3.Y;
                            b.Z = frame.center(i,3)-Nw3.Z;
                            s.X = Nw3.X-Nw2.X;
                            s.Y = Nw3.Y-Nw2.Y;
                            s.Z = Nw3.Z-Nw2.Z;
                            Al = matrix_dot(n,matrix_cross(s,a));

                            PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
                            PB = PA-Al.*matrix_dot(s,m);
                            num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
                            denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
                            phiV = atan2(num,denom);
                            VortexA = VortexA+phiV;
                            %3âÒñ⁄
                            a.X = frame.center(i,1)-Nw3.X;
                            a.Y = frame.center(i,2)-Nw3.Y;
                            a.Z = frame.center(i,3)-Nw3.Z;
                            b.X = frame.center(i,1)-Nw1.X;
                            b.Y = frame.center(i,2)-Nw1.Y;
                            b.Z = frame.center(i,3)-Nw1.Z;
                            s.X = Nw1.X-Nw3.X;
                            s.Y = Nw1.Y-Nw3.Y;
                            s.Z = Nw1.Z-Nw3.Z;
                            Al = matrix_dot(n,matrix_cross(s,a));

                            PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
                            PB = PA-Al.*matrix_dot(s,m);
                            num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
                            denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
                            phiV = atan2(num,denom);
                            VortexA = VortexA+phiV;
                            influence = VortexA./pi./4;
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            A = zeros(1,frame.nbPanel);B = zeros(1,frame.nbPanel);
                            A(1,interpID(1,find(interpID(1,:)))) =  + influence*interpCoeff(find(interpCoeff));
                            B(1,interpID(2,find(interpID(2,:)))) =  - influence*interpCoeff(find(interpCoeff));
                            infA(i,:) = infA(i,:)+A(1,cnt);
                            infA(i,:) = infA(i,:)+B(1,cnt);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        end
                    end
                   if frame.isSymm == 1
                        Nw1.X = frame.verts(frame.tri(j,1),1);
                        Nw1.Y = -frame.verts(frame.tri(j,1),2);
                        Nw1.Z = frame.verts(frame.tri(j,1),3);
                        Nw2.X = frame.verts(frame.tri(j,2),1);
                        Nw2.Y = -frame.verts(frame.tri(j,2),2);
                        Nw2.Z = frame.verts(frame.tri(j,2),3);
                        Nw3.X = frame.verts(frame.tri(j,3),1);
                        Nw3.Y = -frame.verts(frame.tri(j,3),2);
                        Nw3.Z = frame.verts(frame.tri(j,3),3);
                        n.X = frame.n(j,1);
                        n.Y = -frame.n(j,2);
                        n.Z = frame.n(j,3);
                        c.X = frame.center(j,1);
                        c.Y = -frame.center(j,2);
                        c.Z = frame.center(j,3);
                        n12.X = (Nw1.X+Nw2.X)./2;
                        n12.Y = (Nw1.Y+Nw2.Y)./2;
                        n12.Z = (Nw1.Z+Nw2.Z)./2;
                        m = getUnitVector(c,n12);
                        l = matrix_cross(m,n);
                        for i = 1:frame.nbPanel
                            if frame.isBody(i,1)==1
                                pjk.X = frame.center(i,1)-c.X;
                                pjk.Y = frame.center(i,2)-c.Y;
                                pjk.Z = frame.center(i,3)-c.Z;
                                PN = matrix_dot(pjk,n);
                                %1âÒñ⁄
                                a.X = frame.center(i,1)-Nw1.X;
                                a.Y = frame.center(i,2)-Nw1.Y;
                                a.Z = frame.center(i,3)-Nw1.Z;
                                b.X = frame.center(i,1)-Nw2.X;
                                b.Y = frame.center(i,2)-Nw2.Y;
                                b.Z = frame.center(i,3)-Nw2.Z;
                                s.X = Nw2.X-Nw1.X;
                                s.Y = Nw2.Y-Nw1.Y;
                                s.Z = Nw2.Z-Nw1.Z;
                                Al = matrix_dot(n,matrix_cross(s,a));

                                PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
                                PB = PA-Al.*matrix_dot(s,m);
                                num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
                                denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
                                phiV = atan2(num,denom);
                                VortexA = phiV;
                                %2âÒñ⁄
                                a.X = frame.center(i,1)-Nw2.X;
                                a.Y = frame.center(i,2)-Nw2.Y;
                                a.Z = frame.center(i,3)-Nw2.Z;
                                b.X = frame.center(i,1)-Nw3.X;
                                b.Y = frame.center(i,2)-Nw3.Y;
                                b.Z = frame.center(i,3)-Nw3.Z;
                                s.X = Nw3.X-Nw2.X;
                                s.Y = Nw3.Y-Nw2.Y;
                                s.Z = Nw3.Z-Nw2.Z;
                                Al = matrix_dot(n,matrix_cross(s,a));

                                PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
                                PB = PA-Al.*matrix_dot(s,m);
                                num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
                                denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
                                phiV = atan2(num,denom);
                                VortexA = VortexA+phiV;
                                %3âÒñ⁄
                                a.X = frame.center(i,1)-Nw3.X;
                                a.Y = frame.center(i,2)-Nw3.Y;
                                a.Z = frame.center(i,3)-Nw3.Z;
                                b.X = frame.center(i,1)-Nw1.X;
                                b.Y = frame.center(i,2)-Nw1.Y;
                                b.Z = frame.center(i,3)-Nw1.Z;
                                s.X = Nw1.X-Nw3.X;
                                s.Y = Nw1.Y-Nw3.Y;
                                s.Z = Nw1.Z-Nw3.Z;
                                Al = matrix_dot(n,matrix_cross(s,a));

                                PA = matrix_dot(a,matrix_cross(l,matrix_cross(a,s)));
                                PB = PA-Al.*matrix_dot(s,m);
                                num = matrix_dot(s,m).*PN.*(matrix_norm(b).*PA-matrix_norm(a).*PB);
                                denom = PA.*PB+PN.^2.*matrix_norm(a).*matrix_norm(b).*(matrix_dot(s,m)).^2;
                                phiV = atan2(num,denom);
                                VortexA = VortexA+phiV;
                                influence = VortexA./pi./4;
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                A = zeros(1,frame.nbPanel);B = zeros(1,frame.nbPanel);
                                A(1,interpID(1,find(interpID(1,:)))) =  + influence*interpCoeff(find(interpCoeff));
                                B(1,interpID(2,find(interpID(2,:)))) =  - influence*interpCoeff(find(interpCoeff));
                                infA(i,:) = infA(i,:)+A(1,cnt);
                                infA(i,:) = infA(i,:)+B(1,cnt);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            end
                        end
                   end
               end
            end
        end
   else
       infA = [];
       infB = [];             
   end
end


function m = getUnitVector(p1,p2)
m.X = p2.X-p1.X;
m.Y = p2.Y-p1.Y;
m.Z = p2.Z-p1.Z;
no = sqrt(m.X.^2+m.Y.^2+m.Z.^2);
index = find(no);
m.X(index) = m.X(index) ./ no(index);
m.Y(index) = m.Y(index) ./ no(index);
m.Z(index) = m.Z(index) ./ no(index);
end

function no = matrix_norm(m)
no = sqrt(m.X.^2+m.Y.^2+m.Z.^2);
end

function md = matrix_dot(A,B)
md = A.X.*B.X+A.Y.*B.Y+A.Z.*B.Z;
end

function mc = matrix_cross(A,B)
mc.X = A.Y.*B.Z-A.Z.*B.Y;
mc.Y = A.Z.*B.X-A.X.*B.Z;
mc.Z = A.X.*B.Y-A.Y.*B.X;
end



function [darea, dvolume, n] = vartex_iter(P0,P1,P2)
    %ç¿ïWÇó^Ç¶ÇƒSTLÇÃåJÇËï‘Çµïîï™Çí«â¡èëÇ´çûÇ›Ç∑ÇÈ
    x10 = P1(1)-P0(1);y10 = P1(2)-P0(2);z10 = P1(3)-P0(3);
    x20 = P2(1)-P0(1);y20 = P2(2)-P0(2);z20 = P2(3)-P0(3);
    nx=y10*z20-z10*y20;
    ny=z10*x20-x10*z20;
    nz=x10*y20-y10*x20;
    N = norm([nx,ny,nz]);
    darea = sqrt(nx*nx+ny*ny+nz*nz)/2;
    dvolume = (nx*P0(1)+ny*P0(2)+nz*P0(3))/6;
    n = [nx,ny,nz]./N;
end
