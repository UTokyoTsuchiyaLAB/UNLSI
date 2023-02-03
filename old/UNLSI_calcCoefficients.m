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

if ansf.Mach<1 && frame.isSymm == 1 && isfield(frame,'yVertsID')
    dCA = (-ansf.Cp.*frame.n(:,1)+ansf.Cfe.*ansf.s(:,1)).*frame.area./frame.REFS(2,1);
    dCY = (-ansf.Cp.*frame.n(:,2)+ansf.Cfe.*ansf.s(:,2)).*frame.area./frame.REFS(2,1);
    dCN = (-ansf.Cp.*frame.n(:,3)+ansf.Cfe.*ansf.s(:,3)).*frame.area./frame.REFS(2,1);
    dCM = cross(frame.center-repmat(frame.REFS(1,1:3),[size(frame.center,1),1]),[dCA,dCY,dCN]);
    dCMX = dCM(:,1)./frame.REFS(2,2);
    dCMY = dCM(:,2)./frame.REFS(2,3);
    dCMZ = dCM(:,3)./frame.REFS(2,2);
    ansf.CMX = sum(dCMX)*0; %left turn: positive
    ansf.CMY = sum(dCMY)*2; %pitch down: negative
    ansf.CMZ = sum(dCMZ)*0; %left yow: positive
    CN = sum(dCN)*2;
    CA = sum(dCA)*2;
    CY = sum(dCY)*0;
    Cwind = ansf.FlowRot\[CA;CY;CN];
    ansf.CY = Cwind(2);
    ansf.CL = Cwind(3);
    
    Cp = ansf.Cp;
    for i = 1:numel(frame.wakeCalcID)
        Cp(frame.surfID == frame.wakeCalcID(i),1) = 0;
    end
    dCA = (-Cp.*frame.n(:,1)+ansf.Cfe.*ansf.s(:,1)).*frame.area./frame.REFS(2,1);
    dCY = (-Cp.*frame.n(:,2)+ansf.Cfe.*ansf.s(:,2)).*frame.area./frame.REFS(2,1);
    dCN = (-Cp.*frame.n(:,3)+ansf.Cfe.*ansf.s(:,3)).*frame.area./frame.REFS(2,1);
    CN = sum(dCN)*2;
    CA = sum(dCA)*2;
    CY = sum(dCY)*0;
    Cwind = ansf.FlowRot\[CA;CY;CN];
    ansf.CD = Cwind(1);
    ansf.CL_LLT = Cwind(3);
    AR = frame.REFS(2,2)^2/frame.REFS(2,1);
    for i = 1:numel(frame.yVertsID)
        muw = (ansf.u(frame.muBodyID{i}(1,:),1)-ansf.u(frame.muBodyID{i}(2,:),1))';
        [CLi,CDi] = sparseLLT(frame.LLTMat{i},-muw,frame.REFS(2,1),frame.bLLT(i),frame.isSymm);
        ansf.CD = ansf.CD + CDi./sqrt(1-ansf.Mach^2);
        ansf.CL_LLT = ansf.CL_LLT+CLi./sqrt(1-ansf.Mach^2);
        ansf.eff(i) = CLi^2/(pi*AR*CDi);
    end
    
    %{
    muw = -muw;
    b = frame.bLLT(i);optMat = frame.LLTMat{i};
    y = ((frame.verts(frame.yVertsID{i}(1,:),2)+frame.verts(frame.yVertsID{i}(2,:),2))./2)';
    N = numel(muw)*2;
    optRHS(N+1:N+numel(muw),1) = muw'./(2*b);
    Aftbuff = optMat\optRHS;
    Aft = Aftbuff(1:N,1);
    figure(3);clf;set(gcf,'Position',[100 100 840 630]);set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial');hold on;grid on
    s = linspace(0,pi/2,100);
    ys = b/2.*cos(s);
    i = 1;
    g = Aft(i).*sin((2*i-1).*s).*2*b;
    for i = 2:numel(Aft)
        g =g+ Aft(i).*sin((2*i-1).*s).*2*b;
    end
    plot(y,muw.*2,'ko','LineWidth',2);plot(ys,g.*2,'k-','LineWidth',2);
    legend({'CL(on wake)','Interpolation of CL by sparseLLT'},'FontSize',15,'Location','southeast')
    xlabel('y(m)');ylabel('CL')
    %}
elseif ansf.Mach<1 && frame.isSymm == 0 && isfield(frame,'yVertsID')
    dCA = (-ansf.Cp.*frame.n(:,1)+ansf.Cfe.*ansf.s(:,1)).*frame.area./frame.REFS(2,1);
    dCY = (-ansf.Cp.*frame.n(:,2)+ansf.Cfe.*ansf.s(:,2)).*frame.area./frame.REFS(2,1);
    dCN = (-ansf.Cp.*frame.n(:,3)+ansf.Cfe.*ansf.s(:,3)).*frame.area./frame.REFS(2,1);
    dCM = cross(frame.center-repmat(frame.REFS(1,1:3),[size(frame.center,1),1]),[dCA,dCY,dCN]);
    dCMX = dCM(:,1)./frame.REFS(2,2);
    dCMY = dCM(:,2)./frame.REFS(2,3);
    dCMZ = dCM(:,3)./frame.REFS(2,2);
    CN = sum(dCN);
    CA = sum(dCA);
    CY = sum(dCY);
    ansf.CMX = sum(dCMX);
    ansf.CMY = sum(dCMY);
    ansf.CMZ = sum(dCMZ);
    Cwind = ansf.FlowRot\[CA;CY;CN];
    ansf.CL = Cwind(3);
    ansf.CY = Cwind(2);
    
    Cp = ansf.Cp;
    for i = 1:numel(frame.wakeCalcID)
        Cp(frame.surfID == frame.wakeCalcID(i),1) = 0;
    end
    dCA = (-Cp.*frame.n(:,1)+ansf.Cfe.*ansf.s(:,1)).*frame.area./frame.REFS(2,1);
    dCY = (-Cp.*frame.n(:,2)+ansf.Cfe.*ansf.s(:,2)).*frame.area./frame.REFS(2,1);
    dCN = (-Cp.*frame.n(:,3)+ansf.Cfe.*ansf.s(:,3)).*frame.area./frame.REFS(2,1);
    CN = sum(dCN);
    CA = sum(dCA);
    CY = sum(dCY);
    Cwind = ansf.FlowRot\[CA;CY;CN];
    ansf.CD = Cwind(1);
    ansf.CL_LLT = Cwind(3);
    AR = frame.REFS(2,2)^2/frame.REFS(2,1);
    for i = 1:numel(frame.yVertsID)
        muw = (ansf.u(frame.muBodyID{i}(1,:),1)-ansf.u(frame.muBodyID{i}(2,:),1))';
        [CLi,CDi] = sparseLLT(frame.LLTMat{i},-muw,frame.REFS(2,1),frame.bLLT(i),frame.isSymm);
        ansf.CD = ansf.CD + CDi./sqrt(1-ansf.Mach^2);
        ansf.CL_LLT = ansf.CL_LLT+CLi./sqrt(1-ansf.Mach^2);
        ansf.eff(i) = CLi^2/(pi*AR*CDi);
    end
    
else
    dCA = (-ansf.Cp.*frame.n(:,1)+ansf.Cfe.*ansf.s(:,1)).*frame.area./frame.REFS(2,1);
    dCY = (-ansf.Cp.*frame.n(:,2)+ansf.Cfe.*ansf.s(:,2)).*frame.area./frame.REFS(2,1);
    dCN = (-ansf.Cp.*frame.n(:,3)+ansf.Cfe.*ansf.s(:,3)).*frame.area./frame.REFS(2,1);
    dCM = cross(frame.center-repmat(frame.REFS(1,1:3),[size(frame.center,1),1]),[dCA,dCY,dCN]);
    dCMX = dCM(:,1)./frame.REFS(2,2);
    dCMY = dCM(:,2)./frame.REFS(2,3);
    dCMZ = dCM(:,3)./frame.REFS(2,2);
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
    Cwind = ansf.FlowRot\[CA;CY;CN];
    ansf.CD = Cwind(1);
    ansf.CL = Cwind(3);
    ansf.CY = Cwind(2);
    
end
ansf.L_D = ansf.CL./ansf.CD;
end

function [CL,CDi,AR] = sparseLLT(optMat,muw,S,b,isSymm)
    N = numel(muw)*2;
    optRHS(N+1:N+numel(muw),1) = muw'./(2*b);
    Aftbuff = optMat\optRHS;
    Aft = Aftbuff(1:N,1);
    %{
    figure(8);clf;hold on
    s = linspace(0,pi,100);
    for i = 1:size(frame.yVertsID,2)
        y = ((frame.verts(frame.yVertsID{i}(1,:),2)+frame.verts(frame.yVertsID{i}(2,:),2))./2)';
    end
    ys = b/2.*cos(s);
    i = 1;
    g = Aft(i).*sin(i.*s);
    for i = 2:numel(Aft)
        g =g+ Aft(i).*sin(i.*s);
    end
    plot(ys,g);plot(y,muw./(2*b),'o');
    %}
    AR = b^2/S;
    CL = Aft(1)*pi*AR;
    if isSymm == 1
        CDi = pi*AR*sum((1:2:2*N-1)'.*Aft.^2);
    else
        CDi = pi*AR*sum((1:N)'.*Aft.^2);
    end
end
