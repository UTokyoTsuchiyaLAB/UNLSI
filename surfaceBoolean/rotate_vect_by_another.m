%this function gives a rotation matrix to rotate one vector to another

function M = rotate_vect_by_another(u,v)
% create rotation matrix to align v with u
    %Input
    % u and v, vectors in R3
    
    %Output
    % M : Rotation matrix to align one vector with another
    
    u = u/norm(u); %the vector to be rotated
    v = v/norm(v);
    if ~isequal(round(u,4),v) &&  ~all(-round(u,4)==(v)) %equivalent vectors that are not oposite to each other
        ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
        RU = @(A,B) eye(3) + ssc(cross(A,B)) + ssc(cross(A,B))^2*(1-dot(A,B))/(norm(cross(A,B))^2);
        M = RU(v,u);
    elseif all(-round(u,4)==(v))
        %normals are 180 degrees apart
        M = -1*eye(3);
    else
        M = eye(3); 
    end
    
end

