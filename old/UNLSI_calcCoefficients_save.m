function ansf = UNLSI_calcCoefficients(u,ansf,frame,flow)
ansf.u = u;
ansf.Cp = zeros(size(frame.tri,1),1);
if ansf.Mach < 1
    uall = zeros(frame.nPanel,1);
    uall(frame.isBody == 1,1) = u;
    potential = -uall + sum(ansf.Vinf.*frame.center,2);
    ansf.Vx = frame.VxMat*(potential);
    ansf.Vy = frame.VyMat*(potential);
    ansf.Vz = frame.VzMat*(potential);
    ansf.Cp(frame.isBody == 1,1) = (1-sum([ansf.Vx(frame.isBody == 1,1),ansf.Vy(frame.isBody == 1,1),ansf.Vz(frame.isBody == 1,1)].^2,2))./sqrt(1-ansf.Mach^2);
    ansf.Cp(frame.isBase == 1,1) = (-0.139-0.419.*(ansf.Mach-0.161).^2);
else
    ansf.Cp(frame.isBody==1,1) = u;
    ansf.Cp(frame.isBase==1,1) = (-ansf.Mach.^(-2)+0.57.*ansf.Mach.^(-4));
end

ansf.Cfe = zeros(size(ansf.Cp));
ansf.Tw = zeros(size(ansf.Cp));
if isfield(flow,'ppCfe')
   if ansf.Mach<1
       ansf.Cfe(frame.isBody==1,1) = ppval(flow.ppCfe,ansf.Cp(frame.isBody==1,1));
       ansf.Tw(frame.isBody==1,1) = ppval(flow.ppTw,ansf.Cp(frame.isBody==1,1));
   else
       ansf.Cfe(frame.isBody==1,1) = ppval(flow.ppCfe,ansf.delta);
       ansf.Tw(frame.isBody==1,1) = ppval(flow.ppTw,ansf.delta);
   end
end
ansf.Cp(isnan(ansf.Cp)==1) = 0;
ansf.Cfe(isnan(ansf.Cfe)==1) = 0;
frame.n(isnan(frame.n(:,1))==1,1) = 0;
frame.n(isnan(frame.n(:,2))==1,2) = 0;
frame.n(isnan(frame.n(:,3))==1,3) = 0;
ansf.s(isnan(ansf.s(:,1))==1,1) = 0;
ansf.s(isnan(ansf.s(:,2))==1,2) = 0;
ansf.s(isnan(ansf.s(:,3))==1,3) = 0;
dCA = (-ansf.Cp.*frame.n(:,1)+ansf.Cfe.*ansf.s(:,1)).*frame.area./frame.REFS(2,1);
dCY = (-ansf.Cp.*frame.n(:,2)+ansf.Cfe.*ansf.s(:,2)).*frame.area./frame.REFS(2,1);
dCN = (-ansf.Cp.*frame.n(:,3)-ansf.Cfe.*ansf.s(:,3)).*frame.area./frame.REFS(2,1);
dCMX = -dCY.*(frame.center(:,3)-frame.REFS(1,3))./frame.REFS(2,2)+dCN.*(frame.center(:,2)-frame.REFS(1,2))./frame.REFS(2,2);
dCMY = -dCN.*(frame.center(:,1)-frame.REFS(1,1))./frame.REFS(2,3)+dCA.*(frame.center(:,3)-frame.REFS(1,3))./frame.REFS(2,2);
dCMZ = dCY.*(frame.center(:,1)-frame.REFS(1,1))./frame.REFS(2,3)-dCA.*(frame.center(:,2)-frame.REFS(1,2))./frame.REFS(2,2);

if frame.isSymm == 1
    CN = sum(dCN)*2;
    CA = sum(dCA)*2;
    CY = sum(dCY)*0;
    ansf.CMX = sum(dCMX)*0; %left turn: positive
    ansf.CMY = sum(dCMY)*2; %pitch down: negative
    ansf.CMZ = sum(dCMZ)*0; %left yow: positive
else
    CN = sum(dCN);
    CA = sum(dCA);
    CY = sum(dCY);
    ansf.CMX = sum(dCMX);
    ansf.CMY = sum(dCMY);
    ansf.CMZ = sum(dCMZ);
end
%}
Cwind = ansf.FlowRot\[CA;CY;CN];
ansf.CY = Cwind(2);
if ansf.Mach<1
    if isempty(flow.LLTwakeNo)
        ansf.CD = Cwind(1);
        ansf.CL = Cwind(3);
    else
        for i = 1:size(frame.wakelines,2)
            if any(flow.LLTwakeNo == frame.wakelines{i}.ID)
                y(i) = (frame.verts(frame.wakelines{i}.nodeID(1),2)+frame.verts(frame.wakelines{i}.nodeID(2),2))/2;
                ID(i) = frame.wakelines{i}.ID;
                muw(i) = ansf.u(frame.wakelines{i}.upperID)-ansf.u(frame.wakelines{i}.lowerID);
            end
        end
        [~,index] = sort(y);y = y(index);muw = muw(index);ID = ID(index);
        [~,ansf.CD,AR] = sparseLLT(y,-muw,frame.REFS(2,1),frame.REFS(2,2));
        dCA = (ansf.Cfe.*ansf.s(:,1)).*frame.area./frame.REFS(2,1);
        dCY = (ansf.Cfe.*ansf.s(:,2)).*frame.area./frame.REFS(2,1);
        dCN = (-ansf.Cfe.*ansf.s(:,3)).*frame.area./frame.REFS(2,1);
        if frame.isSymm == 1
            CN = sum(dCN)*2;
            CA = sum(dCA)*2;
            CY = sum(dCY)*0;
        else
            CN = sum(dCN);
            CA = sum(dCA);
            CY = sum(dCY);
        end
        Cfewind = ansf.FlowRot\[CA;CY;CN];
        %ansf.CL =  ansf.CL/sqrt(1-ansf.Mach^2);
        ansf.CL = Cwind(3);
        ansf.CD =  ansf.CD/sqrt(1-ansf.Mach^2)+Cfewind(1);
        ansf.efficiency = ansf.CL^2/(pi*AR*ansf.CD);
    end
else
    ansf.CD = Cwind(1);
    ansf.CL = Cwind(3);
end
    
ansf.L_D = ansf.CL./ansf.CD;
end

function [CL,CDi,AR] = sparseLLT(y,muw,S,b)
    theta = acos(y./(b/2));
    N = size(theta,2)*2;
    optMat = zeros(N+size(theta,2));
    for i = 1:N
        optMat(i,i) = pi/2*(2*i-1)^4;
        optMat(i,N+1:N+size(theta,2)) = sin((2*i-1).*theta);
    end
    for j = 1:size(theta,2)
        optMat(N+j,1:N) = sin((2*(1:N)'-1)*theta(j));
    end
    optRHS(N+1:N+size(theta,2),1) = muw'./(2*b);
    Aftbuff = optMat\optRHS;
    Aft = Aftbuff(1:N,1);
    AR = b^2/S;
    CL = Aft(1)*pi*AR;
    CDi = pi*AR*sum((1:2:2*N-1)'.*Aft.^2);
end