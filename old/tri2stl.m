function tri2stl(triname,stlname)

    fp = fopen(triname,'r');
    buff = fscanf(fp,'%d %d',2);
    nVerts = buff(1);
    nTri = buff(2);

    verts = fscanf(fp,'%f %f %f',[3,nVerts])';
    tri = fscanf(fp,'%d %d %d',[3,nTri])';
    surfID = fscanf(fp,'%d',[1,nTri])';
    x = verts(:,1);
    y = verts(:,2);
    z = verts(:,3);
    fclose all;
    
    fp=fopen(stlname,'wt');
	fprintf(fp,'solid %s\n',stlname);
	fclose(fp);
    
    area = 0;
	volume = 0;
    for i = 1:size(tri,1)
        P0 = [verts(tri(i,1),1),verts(tri(i,1),2),verts(tri(i,1),3)];
        P1 = [verts(tri(i,2),1),verts(tri(i,2),2),verts(tri(i,2),3)];
        P2 = [verts(tri(i,3),1),verts(tri(i,3),2),verts(tri(i,3),3)];
        [darea dvolume] = vartex_iter(P0,P1,P2,stlname);
        area = area + darea;
        volume = volume + dvolume;
	end
	fp=fopen(stlname,'at');
	fprintf(fp,'endsolid %s',stlname);
	fclose(fp);
end

function [darea dvolume] = vartex_iter(P0,P1,P2,filename)
    %À•W‚ğ—^‚¦‚ÄSTL‚ÌŒJ‚è•Ô‚µ•”•ª‚ğ’Ç‰Á‘‚«‚İ‚·‚é
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