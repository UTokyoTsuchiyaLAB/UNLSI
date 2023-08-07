function [sp,roll,spiral,dr,dynCoefStruct] = DYNCOEF2Mode(DYNCOEF,rho,UREF,SREF,BREF,CREF,mass,Inatia)
        dynCoefStruct.Cyb = DYNCOEF(2,1);
        dynCoefStruct.Clb = DYNCOEF(4,1);
        dynCoefStruct.Cnb = DYNCOEF(6,1);
        dynCoefStruct.Cxa = DYNCOEF(1,2);
        dynCoefStruct.Cza = DYNCOEF(3,2);
        dynCoefStruct.Cma = DYNCOEF(5,2);
        dynCoefStruct.Cyp = DYNCOEF(2,3);
        dynCoefStruct.Clp = DYNCOEF(4,3);
        dynCoefStruct.Cnp = DYNCOEF(6,3);
        dynCoefStruct.Cxq = DYNCOEF(1,4);
        dynCoefStruct.Czq = DYNCOEF(3,4);
        dynCoefStruct.Cmq = DYNCOEF(5,4);
        dynCoefStruct.Cyr = DYNCOEF(2,5);
        dynCoefStruct.Clr = DYNCOEF(4,5);
        dynCoefStruct.Cnr = DYNCOEF(6,5);
        

	    %Mq = rho*UREF*SREF*CREF*CREF/(4*Inatia(2,2))*dynCoefStruct.Cmq;
	    Yb = rho*UREF*UREF*SREF/(2*mass)*dynCoefStruct.Cyb;
	    Lb = rho*UREF*UREF*SREF*BREF/(2*Inatia(1,1))*dynCoefStruct.Clb;
	    Nb = rho*UREF*UREF*SREF*BREF/(2*Inatia(3,3))*dynCoefStruct.Cnb;
	    %Yp = rho*UREF*SREF*BREF/(4*mass)*dynCoefStruct.Cyp;
	    Lp = rho*UREF*SREF*BREF*BREF/(4*Inatia(1,1))*dynCoefStruct.Clp;
	    Np = rho*UREF*SREF*BREF*BREF/(4*Inatia(3,3))*dynCoefStruct.Cnp;
	    %Yr = rho*UREF*SREF*BREF/(4*mass)*dynCoefStruct.Cyr;
	    Lr = rho*UREF*SREF*BREF*BREF/(4*Inatia(1,1))*dynCoefStruct.Clr;
	    Nr = rho*UREF*SREF*BREF*BREF/(4*Inatia(3,3))*dynCoefStruct.Cnr;
        %Xu = rho*UREF*SREF/(2*mass)*(dynCoefStruct.Cxu);
        %Zu = rho*UREF*SREF/(2*mass)*(dynCoefStruct.CZu-2*dynCoefStruct.CL);
        %Xa = rho*UREF^2*SREF/(2*mass)*dynCoefStruct.Cxa;
        Za = rho*UREF^2*SREF/(2*mass)*dynCoefStruct.Cza;
        %Zq = rho*UREF*SREF*CREF/(4*mass)*(dynCoefStruct.Czq);
        %Mu = rho*UREF*SREF*CREF/(2*Inatia(2,2))*dynCoefStruct.Cmu;
        Ma = rho*UREF^2*SREF*CREF/(2*Inatia(2,2))*dynCoefStruct.Cma;
        Mq = rho*UREF*SREF*CREF^2/(4*Inatia(2,2))*dynCoefStruct.Cmq;

        
        sp.omegan = sqrt(-Ma+Za/UREF*Mq);
        sp.eta =(-Za/UREF-Mq)/2/sp.omegan;
        
        if not(Yb == 0 && Lp == 0)
            roll.T = -1/(Lp);
            lambdas_u = (Lb*Nr-Nb*Lr);
            lambdas_l = (Lp*Nr-Np*Lr)*Yb+(Lb*Nr-Nb*Lr)*UREF-Lb*9.8;
            spiral.T = lambdas_l/lambdas_u;
    
            dr.omegan = sqrt(Nb-(Np/Lp)*Lb);
            dr.eta = -(Nr-(Np/Lp)*Lr+(Np/Lp^2)*Lb)/2/sp.omegan;
        else
            roll = [];
            spiral = [];
            dr = [];
        end
        
end
