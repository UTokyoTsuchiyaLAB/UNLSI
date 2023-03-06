clear 
orgVal = modifyDesFile("org.des","org.des");
[~, orgSurf] = readvspgeom( "orgSurf.vspgeom", 0);
ndim = numel(orgVal);
designScale = [2,1,10,0.5,5,0.02,0.02];
pert = 0.01;
for i = 1:ndim
    sampleDes = orgVal;
    sampleDes(i) = orgVal(i) + pert*(designScale(i));
    modifyDesFile("org.des","mod.des",sampleDes);
    CommandVSP = strcat("vsp wing.vsp3 -script makeModSurf.vspscript"); %手順2のコマンド環境に合わせて各自書き換え
    [~,~] = system(CommandVSP);
    [~, modSurf] = readvspgeom( "modSurf.vspgeom", 0);
    dmodSurf = modSurf'-orgSurf';
    sampleSurff = dmodSurf(:);
    sampleDes = orgVal;
    sampleDes(i) = orgVal(i) - pert*(designScale(i));
    modifyDesFile("org.des","mod.des",sampleDes);
    CommandVSP = strcat("vsp wing.vsp3 -script makeModSurf.vspscript"); %手順2のコマンド環境に合わせて各自書き換え
    [~,~] = system(CommandVSP);
    [~, modSurf] = readvspgeom( "modSurf.vspgeom", 0);
    dmodSurf = modSurf'-orgSurf';
    sampleSurfr = dmodSurf(:);
    surfMod(:,i) = (sampleSurff-sampleSurfr)./2./(pert*(designScale(i)));
end

surfVal = orgVal+designScale.*0.001;
verts = orgSurf'+reshape(surfMod*(surfVal-orgVal)',size(orgSurf'));
modifyDesFile("org.des","mod.des",surfVal);
CommandVSP = strcat("vsp wing.vsp3 -script makeModSurf.vspscript"); %手順2のコマンド環境に合わせて各自書き換え
[~,~] = system(CommandVSP);
[con, testSurf] = readvspgeom( "modSurf.vspgeom", 0);
figure(1);clf;axis equal;
trisurf(con',testSurf(1,:)',testSurf(2,:)',testSurf(3,:)')
hold on
trisurf(con',verts(:,1),verts(:,2),verts(:,3),sum((testSurf'-verts).^2,2));
axis equal;colorbar;