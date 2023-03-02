function [vspSurf, vspCon] = vspSurfGen(x,vspTitlename,desFilename)
            modifyDesFile(desFilename,"mod.des",x);
            CommandVSP = strcat("vsp ",vspTitlename,".vsp3 -script makeModSurf.vspscript -des mod.des"); %手順2のコマンド環境に合わせて各自書き換え
            [~,~] = system(CommandVSP);
            [vspCon, vspSurf] = readvspgeom( "modSurf.vspgeom", 0);
            vspSurf = vspSurf';
            vspCon = vspCon';
end