function [vspSurf, vspCon,SREF,BREF,CREF,XYZREF,argin_x, x] = vspGeomGen(x,vspTitlename,desFilename)
            [~,x] = modifyDesFile(desFilename,"mod.des",x);
            if ispc == 1
                CommandVSP = strcat("vsp ",vspTitlename,".vsp3 -script makeModGeom.vspscript -des mod.des"); %手順2のコマンド環境に合わせて各自書き換え
            else
                CommandVSP = strcat("./vsp ",vspTitlename,".vsp3 -script makeModGeom.vspscript -des mod.des"); %手順2のコマンド環境に合わせて各自書き換え
            end
            [~,~] = system(CommandVSP);
            [vspCon, vspSurf] = readvspgeom( "modGeom.vspgeom", 0);
            vspSurf = vspSurf';
            vspCon = vspCon';
            comp_areas = XYZProj2Mat(3);
            SREF = comp_areas(2);
            BREF = 1.24;
            CREF = 1;
            XYZREF = [0,0,0];
            argin_x = 0;
end