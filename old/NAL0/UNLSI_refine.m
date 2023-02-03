%
clear 
p(1) = 9.1;p(2)= 6.0;p(3) = 2.6;p(4) = 49.5;p(5) = 50.5;p(6) = 22.2;p(7) = 6.67;p(8) = 23.4;p(9) = 0.38;p(10)=11.2;p(11) = 15.6;p(12) = 4.4;p(13)=60*pi/180;p(14)=1/100*0.115;
wgsname = 'NAL0.wgs';
vehicle = NAL0Design(p,wgsname,'UNLSI');

[verts,tri,surfID,surfThn,surfE,surfRho] = NAL02tri(vehicle,p(14));
isSymm = 1;
REFS = [0 0 0 0;1300*p(14),46.8,100,1].*p(14);


frame = UNLSI_makeFrame(verts,tri,surfID);
frame = UNLSI_makeEquations(verts,REFS,isSymm,frame);

%solve flow
Mach = 0.3;
alpha = 15;
flow = UNLSI_makeFlow(Mach,p(14)*100,100000,300);
beta = 0;additionaldata.omega = [0,0,0];additionaldata.omegaCenter = [0,0,0];
[ansf] = UNLSI_calcResidual(frame,flow,alpha,beta,additionaldata);
ansf = UNLSI_calcCoefficients(ansf.u,ansf,frame,flow);

%below:structure 
frame = UNLSI_initFem(frame);
frame = UNLSI_makeFemLHS(frame,surfThn,surfE,surfRho);
fem = UNLSI_makeFemRHS(ansf,frame);
anss = UNLSI_calcFemResidual(frame,fem,50000);

figure(1);clf;
frameCpView(frame.verts+anss.displacement*100,frame,ansf.Cp,1,1,[1,1,1],[30,30],[-0.1,0.1])
frame2stl(frame,'NAL0.stl',1,0,0);
%===========================
%以下検証用
%{
Mach = [0.5,0.7,0.9,1.5,2,3,4,7.1];
alpha = [0,5,10,15];
CLvalid(:,1) = [-0.037 0.165 0.431 0.730]';
CDvalid(:,1) = [0.018 0.028 0.081 0.194]';
CLvalid(:,2) = [-0.039 0.177 0.460 0.768]';
CDvalid(:,2) = [0.022 0.031 0.087 0.207]';
CLvalid(:,3) = [-0.026 0.214 0.526 0.827]';
CDvalid(:,3) = [0.026 0.038 0.107 0.229]';
CLvalid(:,4) = [0.004 0.211 0.426 0.634]';
CDvalid(:,4) = [0.036 0.056 0.111 0.202]';
CLvalid(:,5) = [0.000 0.181 0.355 0.528]';
CDvalid(:,5) = [0.030 0.048 0.096 0.169]';
CLvalid(:,6) = [0.000 0.133 0.263 0.402]';
CDvalid(:,6) = [0.018 0.036 0.070 0.134]';
CLvalid(:,7) = [-0.008 0.108 0.231 0.365]';
CDvalid(:,7) = [0.018 0.032 0.062 0.118]';
CLvalid(:,8) = [-0.004 0.069 0.159 0.267]';
CDvalid(:,8) = [0.018 0.023 0.047 0.088]';
%

for i = 1:8
    if i <=3
        figure(2*i-1);set(gcf,'Position',[100 100 840 630])
        clf;set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial');hold on;grid on
        plot(alpha,CL(:,i),'k-','linewidth',2);
        plot(alpha,CLvalid(:,i),'k-.','linewidth',2);
        legend({'UNLSI(Cp integral)','wind tunnel'},'FontSize',15,'Location','southeast');xlabel('AoA(deg)');ylabel('CL')
        print(strcat('NAL0CL_M',mat2str(2*i-1),'.png'),'-dpng')
        figure(2*i);set(gcf,'Position',[100 100 840 630])
        clf;set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial');hold on;grid on
        plot(alpha,CD(:,i),'k-','linewidth',2);
        plot(alpha,CDvalid(:,i),'k-.','linewidth',2);
        legend({'UNLSI(Cp integral)','wind tunnel'},'FontSize',15,'Location','northwest');xlabel('AoA(deg)');ylabel('CD')
        print(strcat('NAL0CD_M',mat2str(2*i),'.png'),'-dpng')
    else
        figure(2*i-1);set(gcf,'Position',[100 100 840 630])
        clf;set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial');hold on;grid on
        plot(alpha,CL(:,i),'k-','linewidth',2);
        plot(alpha,CLvalid(:,i),'k-.','linewidth',2);
        legend({'UNLSI(Cp integral)','wind tunnel'},'FontSize',15,'Location','southeast');xlabel('AoA(deg)');ylabel('CL')
        print(strcat('NAL0CL_M',mat2str(2*i-1),'.png'),'-dpng')
        figure(2*i);set(gcf,'Position',[100 100 840 630])
        clf;set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial');hold on;grid on
        plot(alpha,CD(:,i),'k-','linewidth',2);
        plot(alpha,CDvalid(:,i),'k-.','linewidth',2);
        legend({'UNLSI(Cp integral)','wind tunnel'},'FontSize',15,'Location','northwest');xlabel('AoA(deg)');ylabel('CD')
        print(strcat('NAL0CD_M',mat2str(2*i),'.png'),'-dpng')
    end
end
%}