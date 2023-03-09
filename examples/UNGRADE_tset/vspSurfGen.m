function [vspSurf, vspCon,SREF,BREF,CREF,XYZREF,argin_x, x] = vspSurfGen(x,vspTitlename,desFilename)
            [~,x] = modifyDesFile(desFilename,"mod.des",x);
            CommandVSP = strcat("vsp ",vspTitlename,".vsp3 -script makeModSurf.vspscript -des mod.des"); %手順2のコマンド環境に合わせて各自書き換え
            [~,~] = system(CommandVSP);
            [vspCon, vspSurf] = readvspgeom( "modSurf.vspgeom", 0);
            vspSurf = vspSurf';
            vspCon = vspCon';
            SREF = XYZProj2Mat(3);
            BREF = 20;
            CREF = 4;
            XYZREF = [0,0,0];
            argin_x = 0;
end