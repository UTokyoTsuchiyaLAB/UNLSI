function orgVal = modifyDesFile(filename,outputname,modVal)
    fid = fopen(filename,"r");
    nrow = fscanf(fid,"%d",1);
    for i = 1:nrow
        dataName{i} = fscanf(fid,"%s:%s:%s:%s:",1);
        orgVal(i) = fscanf(fid,"%f",1);
    end
    fclose(fid);
    if nargin == 2
        modVal = orgVal;
    end
    fid = fopen(outputname,"w");
    fprintf(fid,"%d\n",nrow);
    for i = 1:nrow
        fprintf(fid,"%s %f\n",dataName{i},modVal(i));
    end
    fclose(fid);
end