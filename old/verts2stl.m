function verts2stl(verts,tri,surfID,stlname,HR,isWake,isMirror)
    verts = verts.*HR;
    if not(isWake == 1)
        tri(surfID>200,:)=[];
    end
    fp=fopen(stlname,'wt');
	fprintf(fp,'solid %s\n',stlname);
	fclose(fp);
    
    for i = 1:size(tri,1)
        P0 = [verts(tri(i,1),1),verts(tri(i,1),2),verts(tri(i,1),3)];
        P1 = [verts(tri(i,2),1),verts(tri(i,2),2),verts(tri(i,2),3)];
        P2 = [verts(tri(i,3),1),verts(tri(i,3),2),verts(tri(i,3),3)];
        vartex_iter(P0,P1,P2,stlname);
    end
    if isMirror == 1
        for i = 1:size(tri,1)
            P0 = [verts(tri(i,1),1),-verts(tri(i,1),2),verts(tri(i,1),3)];
            P1 = [verts(tri(i,2),1),-verts(tri(i,2),2),verts(tri(i,2),3)];
            P2 = [verts(tri(i,3),1),-verts(tri(i,3),2),verts(tri(i,3),3)];
            vartex_iter(P1,P0,P2,stlname);
        end
    end
	fp=fopen(stlname,'at');
	fprintf(fp,'endsolid %s',stlname);
	fclose(fp);
end

function [darea dvolume] = vartex_iter(P0,P1,P2,filename)
    %���W��^����STL�̌J��Ԃ�������ǉ��������݂���
    x10 = P1(1)-P0(1);y10 = P1(2)-P0(2);z10 = P1(3)-P0(3);
    x20 = P2(1)-P0(1);y20 = P2(2)-P0(2);z20 = P2(3)-P0(3);
    nx=y10*z20-z10*y20;
    ny=z10*x20-x10*z20;
    nz=x10*y20-y10*x20;
    N = norm([nx,ny,nz]);
    darea = sqrt(nx*nx+ny*ny+nz*nz)/2;
    dvolume = (nx*P0(1)+ny*P0(2)+nz*P0(3))/6;

    fp=fopen(filename,'at');
    fprintf(fp,'facet normal %f %f %f\n',nx/N,ny/N,nz/N);
    fprintf(fp,' outer loop\n');
    fprintf(fp,'  vertex %f %f %f\n',P0(1),P0(2),P0(3));
    fprintf(fp,'  vertex %f %f %f\n',P1(1),P1(2),P1(3));
    fprintf(fp,'  vertex %f %f %f\n',P2(1),P2(2),P2(3));
    fprintf(fp,' endloop\n');
    fprintf(fp,'endfacet\n');
    fclose(fp);
end