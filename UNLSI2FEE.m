function fee = UNLSI2FEE(fee,massPropName,propCalcFlag,Mach,alphaDyn,betaDyn,UREF,feeDirectory)
    
    [mass,cgPoint,Inatia] = readMassPropResult(massPropName);
    fee = fee.flowCondition(1,Mach); %flowNoのマッハ数を指定 
    fee = fee.setRotationCenter(cgPoint); %回転中心の設定
    fee = fee.setCf(1,UREF*fee.CREF/15.01*10^6,fee.CREF,0.052*(10^-5),0); %摩擦係数パラメータの指定
    fee = fee.calcDynCoef(1,alphaDyn,betaDyn,[],propCalcFlag,1);
    DYNCOEF = fee.DYNCOEF; 
    REFS = [fee.SREF,fee.BREF,fee.CREF];
    alpharange1 = -10:2:20;
    betarange1 = -10:2:10;
    index = [3,1,11,12,6,18:20];
    sign = [1,1,1,1,1,-1,1,-1];
    [x1,x2] = ndgrid(alpharange1,betarange1);
    fee = fee.solveFlow(1,x1(:),x2(:),[],propCalcFlag);
    coefdata = fee.AERODATA{1}(:,index).*repmat(sign,[numel(x1),1]);
    alpharange2 = -10:1:20;
    betarange2 = -10:1:10;
    [X1,X2] = ndgrid(alpharange2,betarange2);
    for i = 1:6
        coefrbf{i} = RbfppMake(coefdata(:,1:2),coefdata(:,i+2),1,0.01);
    
        for a = 1:numel(alpharange2)
            for b = 1:numel(betarange2)
                coef{i}(a,b) = execRbfInterp(1,coefrbf{i},[X1(a,b),X2(a,b)]);
            end
        end
        ppCoef{i} = griddedInterpolant(X1,X2,coef{i},'linear','linear');
    end
    
    save aeroRBF.mat ppCoef DYNCOEF REFS mass Inatia UREF -mat;
    if nargin>7
        copyfile("aeroRBF.mat",fullfile(feeDirectory,"aeroRBF.mat"),"f");
    end
end
function [mass,cgPoint,Inatia] = readMassPropResult(massPropName)

    fp = fopen(massPropName,"r");
    while(1)
        data = fgetl(fp);
        if strncmpi(data,"Totals",6)
            massPropTxt = data;
            massPropVec = str2num(massPropTxt(7:end));
            break;
        end
        if data == -1
            break;
        end
    end
    fclose(fp);
    
    mass = massPropVec(1);
    cgPoint = massPropVec(2:4);
    Inatia(1,1) = massPropVec(5);
    Inatia(2,2) = massPropVec(6);
    Inatia(3,3) = massPropVec(7);
    Inatia(1,2) = massPropVec(8);
    Inatia(2,1) = massPropVec(8);
    Inatia(1,3) = massPropVec(9);
    Inatia(3,1) = massPropVec(9);
    Inatia(2,3) = massPropVec(10);
    Inatia(3,2) = massPropVec(10);


end

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

function fi = execRbfInterp(mode,pp,val_interp)

    nSamp = size(pp.val_samp,1);
    nInterp = size(val_interp,1);
    val_interp= (val_interp-repmat(pp.scaleShift(:)',[nInterp,1]))./repmat(pp.scaleWeight(:)',[nInterp,1]);
    Xi = sum(val_interp.^2,2);
    H1 = repmat(Xi,[1,nSamp]);
    Xs = sum(pp.val_samp.^2,2);
    H2 = repmat(Xs',[nInterp,1]);
    M = val_interp*pp.val_samp';
    r = sqrt(H1-2.*M+H2);
    fi = phi1(r,pp.R0)*pp.w(:);
end


