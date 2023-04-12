function [vspMesh, vspCon,surfID,wakeLineID, x] = vspMeshGen(x,vspTitlename,desFilename)
            [~,x] = modifyDesFile(desFilename,"org.des",x);
            if ispc == 1
                CommandVSP = strcat("vsp ",vspTitlename,".vsp3 -script makeOrgSurfandMesh.vspscript -des org.des"); %手順2のコマンド環境に合わせて各自書き換え
            else
                CommandVSP = strcat("./vsp ",vspTitlename,".vsp3 -script makeOrgSurfandMesh.vspscript -des org.des"); %手順2のコマンド環境に合わせて各自書き換え
            end
            [~,~] = system(CommandVSP);
            [con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "orgMesh.vspgeom", 0);
            vspMesh = p';
            vspCon = con';
            surfID = id';
            wakeLineID = wedata;
end