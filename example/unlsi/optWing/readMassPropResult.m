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
