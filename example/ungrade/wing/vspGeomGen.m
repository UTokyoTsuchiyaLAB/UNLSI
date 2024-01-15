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
            %comp_areas = XYZProj2Mat(3);
            SREF = (x(2)+x(4)+2*(x(1)+x(3)+x(5)+x(6)))*2;
            BREF = 20;
            CREF = 4;
            XYZREF = [0,0,0];
            argin_x = 0;
end