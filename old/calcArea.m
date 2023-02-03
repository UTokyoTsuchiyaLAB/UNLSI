function [area,pntX] = calcArea(verts,tri)
area = zeros(size(tri,1),1);pntX = zeros(size(tri,1),1);
for i = 1:size(tri,1)
    [area(i,1)] = vartex_iter(verts(tri(i,1),:),verts(tri(i,2),:),verts(tri(i,3),:));
    pntX(i,1) = mean(verts(tri(i,:),1));
end
    
end
function [darea, dvolume, n] = vartex_iter(P0,P1,P2)
    %À•W‚ğ—^‚¦‚ÄSTL‚ÌŒJ‚è•Ô‚µ•”•ª‚ğ’Ç‰Á‘‚«‚İ‚·‚é
    x10 = P1(1)-P0(1);y10 = P1(2)-P0(2);z10 = P1(3)-P0(3);
    x20 = P2(1)-P0(1);y20 = P2(2)-P0(2);z20 = P2(3)-P0(3);
    nx=y10*z20-z10*y20;
    ny=z10*x20-x10*z20;
    nz=x10*y20-y10*x20;
    N = norm([nx,ny,nz]);
    darea = sqrt(nx*nx+ny*ny+nz*nz)/2;
    dvolume = (nx*P0(1)+ny*P0(2)+nz*P0(3))/6;
    n = [nx,ny,nz]./N;
end
