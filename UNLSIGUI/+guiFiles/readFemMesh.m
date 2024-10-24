function [con,verts,ID] = readFemMesh(filename)

    %if strcmpi(filetype,"stl")
    %    [TR,~,~,ID] = stlread(filename);
    %    con = TR.ConnectivityList;
    %    verts = TR.Points;
    %elseif strcmpi(filetype,"msh")
        fid = fopen(filename,'r');
        while 1
            tline = fgetl(fid);
            if ~ischar(tline); fclose(fid); break; end
            if strncmpi(tline,'$Elements',9)
               [con,ID]= getIDandCon( fid );
            elseif strncmpi(tline,'$Nodes',6)
               verts=getVerts( fid );
            end
        end
    %end

end

function verts = getVerts(fid)
    tline = fgetl(fid);
    nRows = sscanf(tline,'%d');
    vertsBuff= fscanf(fid,'%f',[4,nRows])';
    verts=vertsBuff(:,2:4);
end

function [con,ID] = getIDandCon( fid )
tline = fgetl(fid);
nRows = sscanf(tline,'%d');
for i = 1:nRows
    tline = fgetl(fid);
    n = sscanf(tline, '%d')';
    con(i,:) = n(n(3) + 4:end);
    if n(3) > 0
        tags = n(4:3+n(3));
        if length(tags) >= 1
            ID(i,1) = tags(1);
        end
    end
end
end