function [pp] = RbfppMake(xd,fd,rbfMode,r0,colmax)
    if rbfMode ==4
        phi = @(r,r0)phi4(r,r0);
    elseif rbfMode == 1
        phi = @(r,r0)phi1(r,r0);
    elseif rbfMode == 2
        phi = @(r,r0)phi2(r,r0);
    else
        error('未実装')
    end

    if nargin > 4
        colmax = colmax(:);
        scaleShift = zeros(size(colmax,1),1);
        scaleWeight = colmax-scaleShift;
    else
        for i = 1:size(xd,2)
            scaleShift(i,1) = min(xd(:,i));
            colmax(i,1) = max(xd(:,i));
            scaleWeight(i,1) = colmax(i,1)-scaleShift(i,1);
            if scaleWeight(i,1) == 0
                scaleWeight(i,1) = 1;
            end 
        end
    end
    %xdのスケールの変更
    for i = 1:size(xd,2)
        xd(:,i) = (xd(:,i)-scaleShift(i,1))./scaleWeight(i,1);
    end
    
    X = xd';
    H = sum(xd.^2,2);
    H = repmat(H,[1,size(X,2)]);
    r = sqrt(H'-2.*X'*X+H);
    a = phi(r,r0);
    invR = pinv(a);
    w = invR*fd;
    pp.w = w;
    pp.rbfMode = rbfMode;
    pp.nSample = size(xd,1);
    pp.nDesign = size(xd,2);
    pp.val_samp = xd;
    pp.res_samp = fd;
    pp.R0 = r0;
    pp.scaleShift = scaleShift;
    pp.scaleWeight = scaleWeight;
    
end

function  phi = phi1(r,r0)
    phi = sqrt(r.*r+r0.*r0);
end

function  phi = phi2(r,r0)
    phi = 1./sqrt(r.*r+r0.*r0);
end

function phi = phi4(r,r0)
    phi = exp(-0.5.*r.^2/r0.^2);
end

