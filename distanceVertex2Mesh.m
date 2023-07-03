function [pVError, nError] = distanceVertex2Mesh(mesh, vertex)
% DISTANCEVERTEX2MESH - calculate the distance between vertices and a mesh
%
% Syntax: [pVError, nError] = distanceVertex2Mesh(mesh, vertice)
%
% Inputs:
%   mesh -  the mesh which is used as a reference for the distance
%           calculation. 'mesh' needs to be a structure with two fields
%           called 'vertices' and 'faces', where 'vertices' is a n x 3
%           matrix defining n vertices in 3D space, and faces is a m x 3
%           matrix defining m faces with 3 vertice ids each.
%   vertice(s) -    vertices is a q x 3 matrix defining q vertices in 3D
%                   space
%
% Outputs:
%   pVError -   q x 1 array containing the shortest distance for each of
%               the q input vertices to the surface of the 'mesh'
%   nError -    average normalized error: the error is normalized for a
%               densely sampled mesh of a unit sphere of radius 1 and
%               center 0,0,0
%
%
% Example:
% [x,y,z] = sphere(20); % create a unit sphere of 441 samples
% ball = surf2patch(x,y,z,'triangles');
%                       % create triangulation of vertices
%                       % ball is a structure with faces and vertices
% vertices = 1.1 .* [x(:),y(:),z(:)];
%                       % create list of vertices, where all vertices
%                       % are shifted by 0.1 outwards
% [pVError, nError] = distanceVertex2Mesh(ball, vertices);
%                       % pVError is a list of 441 x 1 with 0.1 per entry
%                       % nError is a single value of 0.1
%
% Other m-files required: none
% Subfunctions:
%   distance3DP2P (calculate distance of point to vertex)
%   distance3DP2E (calculate distance of point to edge)
%   distance3DP2F (calculate distance of point to face)
% MAT-files required: none
%
% See also: surf2patch
%
% Author: Christopher Haccius
% Telecommunications Lab, Saarland University, Germany
% email: haccius@nt.uni-saarland.de
% March 2015; Last revision: 26-March-2015

warning ('off','MATLAB:rankDeficientMatrix'); % turn off warnings for 
            % defficient ranks occuring in point-to-face distance
            % calculation

%vertices = mesh.vertices;
%faces = mesh.faces;
%morita Modified
vertices = mesh.Points;
faces = mesh.ConnectivityList;

[numV,dim] = size(vertices);
[numF,pts] = size(faces);

[tV,dimT] = size(vertex);

if(pts~=3)
    error('Only Triangulations allowed (Faces do not have 3 Vertices)!');
elseif (dim~=dimT || dim~=3)
    error('Mesh and Vertices must be in 3D space!');
end

% initialie minimal distance to infinty
d_min = Inf(tV,1);

% first check: find closest vertex
for c1 = 1:tV % iterate over all test vertices
    for c2 = 1:numV % iterate over all mesh vertices
        v1 = vertex(c1,:);
        v2 = vertices(c2,:);
        d = distance3DP2V(v2,v1); % Euclidean distance
        if d < d_min(c1)
            d_min(c1) = d;
        end
    end
end

% second check: find closest edge
for c1 = 1:tV % iterate over all test vertices
    for c2 = 1:numF % iterate over all faces
        for c3 = 1:2 % iterate over all edges of the face
            for c4 = c3+1:3
                v1 = vertex(c1,:);
                v2 = vertices(faces(c2,c3),:);
                v3 = vertices(faces(c2,c4),:);
                % check if edge is possible
                if ( all(min(v2,v3) < (v1 + d_min(c1))) && ...
                     all(max(v2,v3) > (v1 - d_min(c1))) )
                    d = distance3DP2E(v1,v2,v3);
                    if d < d_min(c1) % d is shorter than previous shortest
                        d_min(c1) = d;
                    end
                end
            end
        end
    end
end

% third check: find closest face
for c1 = 1:tV % iterate over all test vertices
    for c2 = 1:numF % iterate over all faces
        v1 = vertex(c1,:);
        v2 = vertices(faces(c2,1),:);
        v3 = vertices(faces(c2,2),:);
        v4 = vertices(faces(c2,3),:);
        % check if face is possible
        if ( all(min([v2;v3;v4]) < (v1 + d_min(c1))) && ...
             all(max([v2;v3;v4]) > (v1 - d_min(c1))) )
            d = distance3DP2F(v1,v2,v3,v4);
            if d < d_min(c1)
                d_min(c1) = d;
            end

        end
    end
end

pVError = d_min; % output error per vertex is d_min
minS = min(vertices); % get size of mesh
maxS = max(vertices);
s = 1/sqrt(3) * norm(0.5 * (maxS - minS)); % calculate size of mesh
nError = sum(pVError) / (s * tV); % average and normalize error
end % end of function


% Point-to-Point Distance in 3D Space
function dist = distance3DP2V(v1,v2)
    dist = norm(v1-v2);% Euclidean distance
end

% Point-to-LineSegment Distance in 3D Space, Line defined by 2 Points (2,3)
function dist = distance3DP2E(v1,v2,v3)
    d = norm(cross((v3-v2),(v2-v1)))/norm(v3 - v2);
    % check if intersection is on edge
    s = - (v2-v1)*(v3-v2)' / (norm(v3-v2))^2;
    if (s>=0 && s<=1)
        dist = d;
    else
        dist = inf;
    end
end

% Point-to-Face Distance in 3D Space, Face defined by 3 Points (2,3,4)
function dist = distance3DP2F(v1,v2,v3,v4)
    n = cross((v4-v2),(v3-v2)) / norm(cross((v4-v2),(v3-v2)));
    d = abs(n * (v1 - v2)');
    % check if intersection is on face
    n = cross((v4-v2),(v3-v2)) / norm(cross((v4-v2),(v3-v2)));
    f1 = v1 + d * n;
    f2 = v1 - d * n;
    m = [v3-v2;v4-v2]';
    try 
        r1 = m\(f1-v2)';
    catch
        r1 = [inf;inf];
    end
    try
        r2 = m\(f2-v2)';
    catch
        r2 = [inf;inf];
    end
    if ((sum(r1)<=1 && sum(r1)>=0 && all(r1 >=0) && all(r1 <=1)) || ...
            (sum(r2)<=1 && sum(r2)>=0 && all(r2 >=0) && all(r2 <=1)))
        dist = d;
    else
        dist = inf;
    end
end