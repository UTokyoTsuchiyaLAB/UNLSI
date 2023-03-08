function [vspMesh, vspCon,surfID,wakeLineID, x] = vspMeshGen(x,vspTitlename,desFilename)
            [~,x] = modifyDesFile("org.des","org.des",x);
            CommandVSP = strcat("vsp wing.vsp3 -script makeOrgSurfandMesh.vspscript -des org.des"); %手順2のコマンド環境に合わせて各自書き換え
            [~,~] = system(CommandVSP);
            [con, p, uv1, uv2, uv3, wedata, id] = readvspgeom( "orgMesh.vspgeom", 0);
            vspMesh = p';
            vspCon = con';
            surfID = id';
            wakeLineID = wedata;
end