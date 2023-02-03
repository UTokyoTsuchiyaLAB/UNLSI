function flow = UNLSI_makeFlow(Mach,Lmean,Pinf,Tinf)
flow.Mach = Mach;
flow.kappa = 1.4;


method = 'OldTangentCone';
if flow.Mach>1
    %Mach1ˆÈã‚È‚çsplinePP‚ðƒpƒlƒ‹Šp“xdelta‚É‘Î‚µ‚Äì‚Á‚Ä‚¨‚­
    delta = linspace(-pi,pi,500);
    Cp = zeros(size(delta));
    for j =1:size(delta,2)
        if delta(j) >= 0
            if strcmp(method,'TangentConeEdwards')
                Mas = (0.87*flow.Mach-0.554)*sin(delta(j))+0.53;
                Cp(j) = (2 * sin(delta(j))^2)/(1 - 1/4*((Mas^2+5)/(6*Mas^2)));
            elseif strcmp(method,'OldTangentCone')
                Mas = 1.090909*flow.Mach*sin(delta(j))+exp(-1.090909*flow.Mach*sin(delta(j)));
                Cp(j) = (2 * sin(delta(j))^2)/(1 - 1/4*((Mas^2+5)/(6*Mas^2)));
            elseif strcmp(method,'TangentWedge')
                if delta(j)>45.585*pi/180
                    Cp(j)=((1.2*flow.Mach*sin(delta(j))+exp(-0.6*flow.Mach*sin(delta(j))))^2-1.0)/(0.6*flow.Mach^2);
                    %R=1/flow.Mach^2+(flow.kappa+1)/2*delta(j)/sqrt(flow.Mach^2-1);
                elseif delta(j)<0.035
                    Cp(j)=flow.kappa*flow.Mach^2*delta(j)/sqrt(flow.Mach^4-1.0);
                else
                    b = -((flow.Mach^2+2)/flow.Mach^2)-flow.kappa*sin(delta(j))^2;
                    c = (2*flow.Mach^2+1)/flow.Mach^4+((flow.kappa+1)^2/4+(flow.kappa-1)/flow.Mach^2)*sin(delta(j))^2;
                    d = -cos(delta(j))^2/flow.Mach^4;
                    q = (b^2-3*c)/9;
                    rr = (b*(2*b^2-9*c)+27*d)/54;
                    disc = q^3-rr^2;
                    if disc<0
                        Cp(j)=((1.2*flow.Mach*sin(delta(j))+exp(-0.6*flow.Mach*sin(delta(j))))^2-1.0)/(0.6*flow.Mach^2);
                    else
                        r = roots([1,b,c,d]);
                        %ts = asin(sqrt(r(2)));
                        R=r(2);
                        Cp(j) = 4*(flow.Mach^2*R-1)/((flow.kappa+1)*flow.Mach^2);
                    end
                end
            end
            Pe_Pinf= Cp(j)*(flow.kappa*flow.Mach^2)/2+1;
            Te_Tinf = Pe_Pinf^(1/(flow.kappa/(flow.kappa-1)));
            Me = sqrt((2+(flow.kappa+1)*flow.Mach^2)/((flow.kappa+1)*Te_Tinf)-2/(flow.kappa+1));
        else
            Me = bisection(@(M2)pmsolve(flow.Mach,M2,-delta(j),flow.kappa),flow.Mach,300);
            Te_Tinf = (2+(flow.kappa+1)*flow.Mach^2)/(2+(flow.kappa+1)*Me^2);
            Pe_Pinf=(Te_Tinf)^(flow.kappa/(flow.kappa-1));
            Cp(j) = 2/(flow.kappa*flow.Mach^2)*(Pe_Pinf-1);
        end
        if not(Lmean==0)
            flow.Pinf = Pinf;
            flow.Tinf = Tinf;
            [Cfe(j),Tw(j)] = VonDriest(Lmean,Me,Pe_Pinf,Te_Tinf,flow); 
        end
    end
    
    flow.ppCp = interp1(delta,Cp,'spline','pp');
    if not(Lmean==0)
        flow.ppCfe = interp1(delta,Cfe,'spline','pp');
        flow.ppTw = interp1(delta,Tw,'spline','pp');
    end
else
    Cp = linspace(-10,1,100);
    for j = 1:size(Cp,2)
        Pe_Pinf= Cp(j)*(flow.kappa*flow.Mach^2)/2+1;
        Te_Tinf = Pe_Pinf^(1/(flow.kappa/(flow.kappa-1)));
        Me = sqrt((2+(flow.kappa+1)*flow.Mach^2)/((flow.kappa+1)*Te_Tinf)-2/(flow.kappa+1));
        if not(Lmean==0)
            flow.Pinf = Pinf;
            flow.Tinf = Tinf;
            [Cfe(j),Tw(j)] = VonDriest(Lmean,Me,Pe_Pinf,Te_Tinf,flow); 
        end
    end
    if not(Lmean==0)
        flow.ppCfe = interp1(Cp,Cfe,'spline','pp');
        flow.ppTw = interp1(Cp,Tw,'spline','pp');
    end
end
end

function [p info] = bisection(fun,a,b)

% provide the equation you want to solve with R.H.S = 0 form. 
% Write the L.H.S by using inline function
% Give initial guesses.
% Solves it by method of bisection.
% A very simple code. But may come handy
i = 1;
if fun(a)*fun(b)>0 
    p = (a + b)/2;
    info = -1;
else
    p = (a + b)/2;
    err = abs(fun(p));
    while(err > sqrt(eps)) 
       if fun(a)*fun(p)<0 
           b = p;
       else
           a = p;          
       end
        p = (a + b)/2; 
       err = abs(fun(p));
       i = i+1;
    end
    info = 1;
end
end

function res = pmsolve(M1,M2,nu,kappa)
M1nu = sqrt((kappa+1)/(kappa-1))*atan(sqrt((kappa-1)/(kappa+1)*(M1^2-1)))-atan(sqrt(M1^2-1));
M2nu = M1nu+nu;
res = sqrt((kappa+1)/(kappa-1))*atan(sqrt((kappa-1)/(kappa+1)*(M2^2-1)))-atan(sqrt(M2^2-1))-M2nu;
end

function M1nu = Prandtl_Meyer(M1,kappa)
    M1nu = sqrt((kappa+1)/(kappa-1))*atan(sqrt((kappa-1)/(kappa+1)*(M1^2-1)))-atan(sqrt(M1^2-1));
end
    
function [Cfe,Tw] = VonDriest(Lx,Me,Pe_Pinf,Te_Tinf,flow)
Te = flow.Tinf*Te_Tinf;
Pe = flow.Pinf*Pe_Pinf;

Taw =Te.*(1+(flow.kappa-1)./2.*Me.^2);
Tw = Taw;
mu0 = 1.458e-6;
S = 110.4;
R = 287;

sose = sqrt(flow.kappa.*R.*Te);
Ve = sose.*Me;
rhoe = Pe./(R.*Te);

%panel.mue = mu0.*(panel.Te./T0).*(T0+S)./(panel.Te+S);
%panel.muw = mu0.*(panel.Tw./T0).*(T0+S)./(panel.Tw+S);
mue = mu0*Te^(1.5)/(Te+S);
muw = mu0*Tw^(1.5)/(Tw+S);

Rexe = rhoe.*Ve.*Lx./mue;

a = sqrt((flow.kappa-1)./2.*Me.^2.*Te./Tw);
b = Taw./Tw-1;
A = (2.*a.^2-b)./sqrt(b.^2+4.*a.^2);
B = b./sqrt(b.^2+4.*a.^2);
options = optimset('Display','off');
Cfe = real(fsolve(@(Cfe)VDsolve(Cfe,A,B,Taw,Te,Rexe,mue,muw),sqrt(eps),options)); 
end

function res = VDsolve(Cfe,A,B,Taw,Te,Rexe,mue,muw)
    res = asin(A)+asin(B)-(sqrt(Cfe*(Taw/Te-1)))*(4.15*log10(Rexe*Cfe*mue/muw)+1.7);
end
