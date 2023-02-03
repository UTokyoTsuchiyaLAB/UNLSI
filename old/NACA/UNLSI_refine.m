%
clear all
%

isSymm = 0;
frame.triname = 'NACA0012.tri';
wgsname = 'NAC0012.wgs';

fp=fopen('naca0012.dat','r');
buff=fgetl(fp);
chordbuff=(fscanf(fp,'%f',[2,350]))';
fclose(fp);
[a, nLE] = min(chordbuff(:,1));
n_wing = 7;
sb = 3;
chordbuff = flipud(chordbuff);
chordbuff(:,2) = chordbuff(:,2);
Ybuff = linspace(0,sb,n_wing);
%Ybuff = sin(sbuff).*sb;
panelfoil.Y = repmat(Ybuff,[size(chordbuff,1),1]);
panelfoil.X = repmat(chordbuff(:,1),[1,n_wing]);
panelfoil.Z = repmat(chordbuff(:,2),[1,n_wing]);
panelfoil1.Y = repmat(Ybuff,[size(chordbuff(1:nLE,1),1),1]);
panelfoil1.X = repmat(chordbuff(1:nLE,1),[1,n_wing]);
panelfoil1.Z = repmat(chordbuff(1:nLE,2),[1,n_wing]);
panelfoil2.Y = repmat(Ybuff,[size(chordbuff(nLE:end,1),1),1]);
panelfoil2.X = repmat(chordbuff(nLE:end,1),[1,n_wing]);
panelfoil2.Z = repmat(chordbuff(nLE:end,2),[1,n_wing]);
panelfoil.X(end,:) = panelfoil.X(1,:);panelfoil.Y(end,:) = panelfoil.Y(1,:);panelfoil.Z(end,:) = panelfoil.Z(1,:);
for i = 1:n_wing
    foilwall{i}.X = [chordbuff(1:nLE,1),flipud(chordbuff(nLE:end,1))];
    foilwall{i}.Z = [chordbuff(1:nLE,2),flipud(chordbuff(nLE:end,2))];
    foilwall{i}.Y = [Ybuff(i).*ones(size(chordbuff(1:nLE,1))),Ybuff(i).*ones(size(chordbuff(nLE:end,1)))];
end
wake.X = [panelfoil.X(end,:);panelfoil.X(end,:)+100];
wake.Y = [panelfoil.Y(end,:);panelfoil.Y(end,:)];
wake.Z = [panelfoil.Z(end,:);panelfoil.Z(end,:)];
buff = wgs_make(wgsname,1,'wing1',1,panelfoil1.X,panelfoil1.Y,panelfoil1.Z,[0 0 0],[1 0],1);
buff = wgs_make(wgsname,1,'wing2',0,panelfoil2.X,panelfoil2.Y,panelfoil2.Z,[0 0 0],[1 0],1);
buff = wgs_make(wgsname,1,'wingwall',0,foilwall{end}.X,foilwall{end}.Y,foilwall{end}.Z,[0 0 0],[1 0],1);
buff = wgs_make(wgsname,18,'wing_wake',0,wake.X,wake.Y,wake.Z,[0 0 0],[1 0],1);
for i = 1:n_wing-1
    buff = wgs_make(wgsname,1,'wingwall',0,foilwall{i}.X,foilwall{i}.Y,foilwall{i}.Z,[0 0 0],[1 0],1);
end
[verts,tri,surfID]=wgs2tri(wgsname,frame.triname,n_wing-1,1,not(isSymm));
%}
%tri(surfID==2,:) = [tri(surfID==2,3),tri(surfID==2,2),tri(surfID==2,1)];

%
figure(4);clf;hold on;
trimesh(tri(surfID <=100 ,:),verts(:,1),verts(:,2),verts(:,3))
%trimesh(tri(surfID == 1 |surfID == 3 | surfID >300,:),verts(:,1),-verts(:,2),verts(:,3))
colormap([0 0 0]);axis equal;axis off;
%{
%surfID(surfID == 2) = 101;
wgsname = 'NAC0012p.wgs';
buff = wgs_make(wgsname,1,'wing',1,panelfoil.X,panelfoil.Y,panelfoil.Z,[0 0 0],[1 0],1);
buff = wgs_make(wgsname,1,'wingwall',0,foilwall{end}.X,foilwall{end}.Y,foilwall{end}.Z,[0 0 0],[1 0],1);
buff = wgs_make(wgsname,18,'wing_wake',0,wake.X,wake.Y,wake.Z,[0 0 0],[0 0],1);
[~,~,~,surfaceX,surfaceY,surfaceZ]=wgs2tri(wgsname,frame.triname,0,1,0);
REFS = [0 0 0 0;sb*2,sb*2,1,1];isSymm = 1;
hypername = 'NAC0012h.wgs';
buff = wgs_make(hypername,1,'wing',1,panelfoil.X,panelfoil.Y,panelfoil.Z,[0 0 0],[1 0],1);
buff = wgs_make(hypername,1,'wingwall',0,foilwall{end}.X,foilwall{end}.Y,foilwall{end}.Z,[0 0 0],[1 0],1);
%}

%[tri, verts, surfID] = readtri('BWB.tri');
%
frame = UNLSI_makeFrame(verts,tri,surfID,[1,201;1,251]);
frame.isSymm = isSymm;
%frame2stl(frame,'naca.stl',1,0,0);
REFS = [0 0 0 0;sb*2,sb*2,1,1];
frame = UNLSI_makeEquations(verts,REFS,isSymm,frame);
Mach = 0.0;
flow = UNLSI_makeFlow(Mach,0,100000,300);
%}
%
alphai = 0:15;
clear CL CD CMY
for i = 1:numel(alphai)
    disp(alphai(i))
    alpha =alphai(i);beta = 0;additionaldata.omega = [0,0,0];additionaldata.omegaCenter = [0,0,0];
    [ansf] = UNLSI_calcResidual(frame,flow,alpha,beta,additionaldata);
    ansf = UNLSI_calcCoefficients(ansf.u,ansf,frame,flow);
    CL(1,i) = ansf.CL;
    CD(1,i) = ansf.CD;
    CMY(1,i) = ansf.CMY;
    if i == 6
        figure(5);clf;
        frameCpView(frame.verts,frame,ansf.Cp,1,5,[1,1,1],[90,90],[-0.5,0.5])
    end
    %{
    if Mach <1
        %[info, now_work] = PANAIRstart(1,wgsname,Mach,alpha,REFS,1,0,isSymm);
        %[CpData_de0{1}.PC] = PANAIRCp(1);
        %[status, presult] = PANAIRend(1,wgsname,Mach,alpha,REFS,0,0,now_work);

        %CL(2,i) = presult.CL;
        %CD(2,i) = presult.CD;
        %CMY(2,i) = presult.CMY;

        %{
        if i == 6
            figure(4);clf;set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial');hold on;grid on

            ndata = size(CpData_de0{1}.PC{1}.X(:,1),1);
            plot([frame.center(1:2:ndata-1,1);frame.center(547:2:547+ndata-1,1)],[ansf.Cp(1:2:ndata-1,1);ansf.Cp(547:2:547+ndata-1,1)],'k-','linewidth',2);
            plot(CpData_de0{1}.PC{1}.X(:,1),CpData_de0{1}.PC{1}.Cp{1}(:,1),'k-.','linewidth',2)
            ax = gca;
            ax.YDir = 'reverse';xlabel('X');ylabel('Cp')
            legend({'UNLSI','PANAIR'},'FontSize',15,'Location','northeast');
            figure(5);clf;
            frameCpView(frame.verts,frame,ansf.Cp,1,5,[2,1,1],[ 90,90])
            caxis([-1.5,1]);
            frameCpView(frame.verts,frame,ansf.Cp,1,5,[2,1,2],[-90,-90])
            caxis([-1.5,1]);
            figure(6);clf;set(gcf,'Position',[100 100 840 630]);set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial')
            for j = 1
                ppCp{j} = RbfppMake([CpData_de0{1}.PC{j}.X(:) CpData_de0{1}.PC{j}.Y(:) CpData_de0{1}.PC{j}.Z(:)],CpData_de0{1}.PC{j}.Cp{1}(:),1,0.0001);
                [CpEdge{j}] = execRbfInterp(1,ppCp{j},[surfaceX{j}(:),surfaceY{j}(:),surfaceZ{j}(:)]);
                CpEdge{j} = reshape(CpEdge{j},size(surfaceX{j}));
                subplot(2,1,1);hold on;view([ 90 90]);caxis([-1.5,1]);colormap jet;colorbar
                surf(surfaceX{j},surfaceY{j},surfaceZ{j},CpEdge{j})
                subplot(2,1,2);hold on;view([-90 -90]);caxis([-1.5,1]);colormap jet;colorbar
                surf(surfaceX{j},surfaceY{j},surfaceZ{j},CpEdge{j})
            end
        end
    %}
    else
        [hresult,CpData_de0{1}.PC] = hyper(hypername,Mach,alphai(i),0,REFS,[17 17],[3,3],1);
        CL(2,i) = hresult.CL;
        CD(2,i) = hresult.CD;
        CMY(2,i) = hresult.CMY;
        if i == 6
            figure(4);clf;set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial');hold on;grid on

            ndata = size(CpData_de0{1}.PC{1}.X(:,1),1);
            plot([frame.center(1:2:ndata-1,1);frame.center(547:2:547+ndata-1,1)],[ansf.Cp(1:2:ndata-1,1);ansf.Cp(547:2:547+ndata-1,1)],'k-','linewidth',2);
            plot(CpData_de0{1}.PC{1}.X(:,1),CpData_de0{1}.PC{1}.Cp{1}(:,1),'k-.','linewidth',2)
            ax = gca;
            ax.YDir = 'reverse';xlabel('X');ylabel('Cp')
            legend({'UNLSI','PANAIR'},'FontSize',15,'Location','southeast');
            figure(5);clf;
            frameCpView(frame.verts,frame,ansf.Cp,1,5,[2,1,1],[ 90,90])
            caxis([0,2]);
            frameCpView(frame.verts,frame,ansf.Cp,1,5,[2,1,2],[-90,-90])
            caxis([0,2]);
            figure(6);clf;set(gcf,'Position',[100 100 840 630]);set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial')
            for j = 1
                ppCp{j} = RbfppMake([CpData_de0{1}.PC{j}.X(:) CpData_de0{1}.PC{j}.Y(:) CpData_de0{1}.PC{j}.Z(:)],CpData_de0{1}.PC{j}.Cp{1}(:),1,0.0001);
                [CpEdge{j}] = execRbfInterp(1,ppCp{j},[surfaceX{j}(:),surfaceY{j}(:),surfaceZ{j}(:)]);
                CpEdge{j} = reshape(CpEdge{j},size(surfaceX{j}));
                subplot(2,1,1);hold on;view([ 90 90]);caxis([0,2]);colormap jet;colorbar
                surf(surfaceX{j},surfaceY{j},surfaceZ{j},CpEdge{j})
                surf(surfaceX{j},-surfaceY{j},surfaceZ{j},CpEdge{j})
                subplot(2,1,2);hold on;view([-90 -90]);caxis([0,2]);colormap jet;colorbar
                surf(surfaceX{j},surfaceY{j},surfaceZ{j},CpEdge{j})
                surf(surfaceX{j},-surfaceY{j},surfaceZ{j},CpEdge{j})
            end
        end
    end
    %}
    disp([alphai(1:i);CL])

end
%}
%{

alphai = 5;
clear CL CD CMY
for i = 1:numel(alphai)
    disp(alphai(i))
    alpha =alphai(i);beta = 0;additionaldata.omega = [0,0,0];additionaldata.omegaCenter = [0,0,0];
    [ansf] = UNLSI_calcResidual(frame,flow,alpha,beta,additionaldata);
    ansf = UNLSI_calcCoefficients(ansf.u,ansf,frame,flow);
end
%below:structure
E_material =  73500000000;rho_material = 2770;stress_yield=470*10^6; %A2000Œn
surfThn = ones(size(surfID)).*0.1./1000;
surfE = ones(size(surfID)).*E_material;
surfRho = ones(size(surfID)).*rho_material;
frame = UNLSI_initFem(frame);
frame = UNLSI_makeFemLHS(frame,surfThn,surfE,surfRho);
fem = UNLSI_makeFemRHS(ansf,frame);
anss = UNLSI_calcFemResidual(frame,fem,0.5*1.225*10^2);
figure(1);clf;
frameCpView(frame.verts+anss.displacement.*50,frame,anss.vonMises./stress_yield,1,1,[1,1,1],[0,0])

%—ƒŒ^‚Ì’f–Ê“ñŽŸƒ‚[ƒƒ“ƒg‚ðŒvŽZ‚·‚é
figure(2);clf;set(gcf,'Position',[100 100 840 630]);set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial');hold on;grid on
Hardness = airfoil_hardness(E_material,0.1/1000,1,chordbuff,0.5);
index = 40*(0:n_wing-1)+1;
dispdata = [frame.verts(index,:),anss.delta_Press(index,:)];
plot(dispdata(:,2),dispdata(:,6).*1000,'ko-','LineWidth',2);

w = 0.5*1.225*10^2.*ansf.CL;
ya = linspace(0,sb,100);
anlydisp = w/24/Hardness(2).*(3.*sb^4-4.*sb^3.*(sb-ya)+(sb-ya).^4);
plot(ya,anlydisp.*1000,'k--','LineWidth',2');legend({'UNLSI(Leading Edge)','Solution of Cantilever'},'FontSize',15,'Location','southeast');xlabel('y(m)');ylabel('Deflection(mm)')

%}
if Mach<1
    CLa = 2*pi/(1+2/sb/2);
    CLapx = CLa.*alphai.*pi./180;
    CDapx = CLapx.^2/pi/0.863/(sb*2);
    figure(1);
    clf;set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial');hold on;grid on
    plot(alphai,CL(1,:),'k-','linewidth',2);
    %plot(alphai,CL(2,:),'k-.','linewidth',2);
    plot(alphai,CLapx,'k:','linewidth',2);
    legend({'UNLSI(sparseLLT)','statistical'},'FontSize',15,'Location','southeast');xlabel('AoA(deg)');ylabel('CL')
    
    figure(2);
    clf;set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial');hold on;grid on
    plot(alphai,CD(1,:),'k-','linewidth',2);
    %plot(alphai,CD(2,:),'k-.','linewidth',2);
    plot(alphai,CDapx,'k:','linewidth',2);
    legend({'UNLSI(sparseLLT)','statistical'},'FontSize',15,'Location','northwest');xlabel('AoA(deg)');ylabel('CD')

    figure(3);
    clf;set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial');hold on;grid on
    plot(alphai,CMY(1,:),'k-','linewidth',2);
    %plot(alphai,CMY(2,:),'k-.','linewidth',2);
    legend({'UNLSI(sparse LLT)'},'FontSize',15);xlabel('AoA(deg)');ylabel('CMY')
else
    figure(1);
    clf;set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial');hold on;grid on
    plot(alphai,CL(1,:),'k-','linewidth',2);
    plot(alphai,CL(2,:),'k-.','linewidth',2);
    legend({'UNLSI(Cp integral)','HYPER'},'FontSize',15,'Location','southeast');xlabel('AoA(deg)');ylabel('CL')
    
    figure(2);
    clf;set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial');hold on;grid on
    plot(alphai,CD(1,:),'k-','linewidth',2);
    plot(alphai,CD(2,:),'k-.','linewidth',2);
    legend({'UNLSI(Cp integral)','HYPER',},'FontSize',15,'Location','northwest');xlabel('AoA(deg)');ylabel('CD')

    figure(3);
    clf;set(gca, 'LineWidth',1.0, 'FontSize',13, 'Fontname','Arial');hold on;grid on
    plot(alphai,CMY(1,:),'k-','linewidth',2);
    plot(alphai,CMY(2,:),'k-.','linewidth',2);
    legend({'UNLSI(Cp integral)','HYPER'},'FontSize',15);xlabel('AoA(deg)');ylabel('CMY')
end

%frameCpView(frame.verts,frame,ansf.Cp,1,10)
%}
%}
