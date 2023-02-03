clear all
%{
clear all
[verts,tri,surfID]=wgs2tri('HIMICO.wgs','HIMICO.tri',3,1);
REFS = [0.8775,0,0,0;0.1929,0.3700,1.5,1];isSymm = 0;

frame = UNLSI_makeFrame(verts,tri,surfID);

frame = UNLSI_makeEquations(verts,REFS,isSymm,frame);
Mach = 5.0;
flow = UNLSI_makeFlow(Mach);
save HIMICO.mat
%}
load HIMICO.mat
vari = [-180:5:180];

for i = 1:numel(vari)
    disp(vari(i))
    alpha =vari(i);beta = 0;omega = [0,0,0];omegaCenter = [0.8775,0,0];
    [ansf] = UNLSI_calcResidual(frame,flow,alpha,beta,omega,omegaCenter);
    ansf = UNLSI_calcCoefficients(ansf.u,ansf,frame,flow);
    CL(1,i) = ansf.CL;
    CD(1,i) = ansf.CD;
    CY(1,i) = ansf.CD;
    CMX(1,i) =  ansf.CMX;
    CMY(1,i) = -ansf.CMY;
    CMZ(1,i) =  ansf.CMZ;
    disp([vari(1:i);CMZ])
end
%{
%below:structure
frame = UNLSI_initFem(frame);
frame = UNLSI_makeFemLHS(frame,surfThn,surfE,surfRho);
fem = UNLSI_makeFemRHS(ansf,frame);
[anss,delta_p] = UNLSI_calcFemResidual(frame,fem,5000);
frameCpView(frame.verts,frame,ansf.Cp,1,9,[-0.3,0.3])
%}
figure(1);clf;plot(vari,CL);legend('UNLSI')
figure(2);clf;plot(vari,CD);legend('UNLSI')
figure(3);clf;plot(vari,CMY);legend('UNLSI')
figure(4);clf;trimesh(tri(surfID<=100,:),verts(:,1),verts(:,2),verts(:,3));axis equal;colormap([0,0,0])
%frameCpView(frame.verts,frame,ansf.u,1,9)
%frameCpView(frame.verts,frame,ansf.Cp,1,10)
